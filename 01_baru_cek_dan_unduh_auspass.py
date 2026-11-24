#!/usr/bin/env python3
"""
01_cek_dan_unduh_auspass.py
Cek ketersediaan continuous waveform (AusPass / FDSN) lalu unduh per hari.

Instalasi (Python 3.10-3.12 disarankan):
    pip install obspy numpy scipy matplotlib pandas requests

Tiga langkah (jalankan berurutan, di komputer yang punya internet):
    1) discover      -> cari stasiun yang epoch-nya menutupi seluruh jendela waktu
    2) availability  -> hitung % cakupan data per stasiun per hari (FDSN availability)
    3) download      -> unduh MiniSEED mentah per hari + StationXML (level=response)
                        + katalog event (QuakeML) + log CSV (sampling rate, npts, gap)

Contoh (Opsi A: sekuens Woods Point, VIC, 21-28 Sep 2021 UTC):
    python 01_cek_dan_unduh_auspass.py discover     --start 2021-09-21 --days 7 --network S1 --maxdist-km 300
    python 01_cek_dan_unduh_auspass.py availability --start 2021-09-21 --days 7
    python 01_cek_dan_unduh_auspass.py download     --start 2021-09-21 --days 7 --stations S1.AUGRF,S1.AUMTS

Catatan:
- Data mentah TIDAK diubah; preprocessing disimpan terpisah (skrip berikutnya).
- Katalog event di AusPass adalah salinan katalog Geoscience Australia yang diperbarui
  manual sekitar sebulan sekali dan tidak memuat phase pick.
- Skrip hanya menulis ke --outdir (default ./data). Menjalankan ulang dengan parameter
  yang sama menimpa file dengan isi yang sama (reproducible).
"""
import argparse
import csv
import hashlib
import sys
import time
from pathlib import Path

import requests
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

BASE_URL = "https://auspass.edu.au"
DAY = 86400.0


# ----------------------------------------------------------------- util
def window(args):
    t0 = UTCDateTime(args.start)
    return t0, t0 + args.days * DAY


def ref_event(client, t0, t1, args):
    """Titik acuan jarak: koordinat manual, atau event terbesar di jendela waktu."""
    if args.ref_lat is not None and args.ref_lon is not None:
        return args.ref_lat, args.ref_lon, "koordinat manual"
    cat = client.get_events(starttime=t0, endtime=t1, minmagnitude=args.ref_minmag,
                            orderby="magnitude", limit=5)
    if not len(cat):
        sys.exit("Tidak ada event acuan di katalog. Beri --ref-lat dan --ref-lon.")
    ev = cat[0]
    o = ev.preferred_origin() or ev.origins[0]
    m = ev.preferred_magnitude() or (ev.magnitudes[0] if ev.magnitudes else None)
    desc = f"{o.time} M{m.mag if m else '?'} (katalog {args.server})"
    return o.latitude, o.longitude, desc


def read_candidates(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# ------------------------------------------------------------ 1) discover
def cmd_discover(args):
    client = Client(args.server)
    t0, t1 = window(args)
    lat, lon, desc = ref_event(client, t0, t1, args)
    print(f"Acuan: lat={lat:.3f} lon={lon:.3f} | {desc}")
    # startbefore/endafter: epoch kanal harus menutupi SELURUH jendela waktu
    inv = client.get_stations(network=args.network, channel=args.channels,
                              startbefore=t0, endafter=t1, level="channel",
                              latitude=lat, longitude=lon,
                              maxradius=kilometers2degrees(args.maxdist_km))
    rows = []
    for net in inv:
        for sta in net:
            for cha in sta:
                d_km = gps2dist_azimuth(lat, lon, sta.latitude, sta.longitude)[0] / 1000.0
                rows.append(dict(net=net.code, sta=sta.code, loc=cha.location_code or "",
                                 cha=cha.code, lat=sta.latitude, lon=sta.longitude,
                                 elev=sta.elevation, sample_rate=cha.sample_rate,
                                 epoch_start=cha.start_date, epoch_end=cha.end_date,
                                 dist_km=round(d_km, 1)))
    rows.sort(key=lambda r: r["dist_km"])
    out = Path(args.outdir) / "meta"
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "candidate_stations.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["info"])
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} kanal kandidat -> {out/'candidate_stations.csv'}")
    for r in rows[:25]:
        print(f"{r['net']}.{r['sta']}.{r['loc']}.{r['cha']:4s} {r['dist_km']:7.1f} km "
              f"{r['sample_rate']} Hz  {r['lat']:.3f},{r['lon']:.3f}")


# ------------------------------------------------------- 2) availability
def fetch_availability(net, sta, loc, cha, t0, t1):
    r = requests.get(f"{BASE_URL}/fdsnws/availability/1/query",
                     params=dict(network=net, station=sta, location=loc or "*",
                                 channel=cha, starttime=t0.isoformat(),
                                 endtime=t1.isoformat(), format="text",
                                 merge="samplerate,quality,overlap", mergegaps=1.0),
                     timeout=120)
    if r.status_code == 204:
        print(f"   [info] HTTP 204 (server: tidak ada data) untuk {net}.{sta}.{cha}")
        return []
    r.raise_for_status()
    out = []
    for line in r.text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split()
        # location boleh kosong (7 kolom) atau "--" (8 kolom): ambil dari BELAKANG
        # kolom terakhir: ... SampleRate Earliest Latest
        try:
            out.append(dict(loc="", cha=cha, sr=float(p[-3]),
                            t0=UTCDateTime(p[-2]), t1=UTCDateTime(p[-1])))
        except Exception:
            print(f"   [peringatan] baris tak terbaca: {line!r}")
    if not out:
        print(f"   [info] HTTP {r.status_code}, 0 baris terbaca. Contoh respons: {r.text[:200]!r}")
    return out


def cmd_availability(args):
    t0, t1 = window(args)
    cands = read_candidates(Path(args.outdir) / "meta" / "candidate_stations.csv")
    cands = [c for c in cands if c["cha"].endswith("Z")][: args.max_stations]
    seen, summary = set(), []
    for c in cands:
        key = (c["net"], c["sta"], c["loc"], c["cha"])
        if key in seen:
            continue
        seen.add(key)
        segs = fetch_availability(c["net"], c["sta"], c["loc"], c["cha"], t0, t1)
        per_day = []
        for d in range(args.days):
            ds, de = t0 + d * DAY, t0 + (d + 1) * DAY
            cov = sum(max(0.0, min(s["t1"], de) - max(s["t0"], ds)) for s in segs)
            per_day.append(round(100 * cov / DAY, 2))
        total = round(sum(per_day) / args.days, 2)
        summary.append(dict(net=c["net"], sta=c["sta"], loc=c["loc"], cha=c["cha"],
                            dist_km=c["dist_km"], n_segments=len(segs),
                            coverage_pct=total,
                            **{f"day{d+1}_pct": v for d, v in enumerate(per_day)}))
        print(f"{key}: {total:6.2f}% | {per_day} | {len(segs)} segmen")
        time.sleep(0.3)
    out = Path(args.outdir) / "meta" / "availability_summary.csv"
    if summary:
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(sorted(summary, key=lambda r: (-r["coverage_pct"], r["dist_km"])))
        print(f"Ringkasan -> {out}")
    print("Pilih stasiun dengan coverage_pct mendekati 100 pada SEMUA hari.")


# --------------------------------------------------------- 3) download
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cmd_download(args):
    client = Client(args.server)
    t0, t1 = window(args)
    raw = Path(args.outdir) / "raw"
    meta = Path(args.outdir) / "meta"
    raw.mkdir(parents=True, exist_ok=True)
    meta.mkdir(parents=True, exist_ok=True)
    log_rows = []

    for code in args.stations.split(","):
        net, sta = code.strip().split(".")
        # metadata respons instrumen (StationXML level=response)
        inv = client.get_stations(network=net, station=sta, channel=args.channels,
                                  startbefore=t0, endafter=t1, level="response")
        inv.write(str(meta / f"{net}.{sta}.response.xml"), format="STATIONXML")

        for d in range(args.days):
            ds, de = t0 + d * DAY, t0 + (d + 1) * DAY
            status, st = "OK", None
            for attempt in range(1, 4):
                try:
                    st = client.get_waveforms(net, sta, args.location, args.channels, ds, de)
                    break
                except Exception as e:  # 204/timeout/dll.
                    status = f"GAGAL: {type(e).__name__}: {str(e)[:80]}"
                    time.sleep(3 * attempt)
            if st is None or len(st) == 0:
                log_rows.append(dict(file="", net=net, sta=sta, loc="", cha="",
                                     start=ds, end=de, sampling_rate="", npts=0,
                                     n_gaps="", gap_total_s="", status=status, sha256=""))
                print(f"{net}.{sta} {ds.date}: {status}")
                continue
            # simpan per kanal (satu file per net.sta.loc.cha.hari)
            for cha in sorted({t.stats.channel for t in st}):
                sub = st.select(channel=cha)
                loc = sub[0].stats.location or "--"
                fn = raw / f"{net}.{sta}.{loc}.{cha}.{ds.strftime('%Y-%m-%d')}.mseed"
                sub.write(str(fn), format="MSEED")
                g_sub = sub.get_gaps()
                log_rows.append(dict(
                    file=fn.name, net=net, sta=sta, loc=loc, cha=cha,
                    start=min(t.stats.starttime for t in sub),
                    end=max(t.stats.endtime for t in sub),
                    sampling_rate=sub[0].stats.sampling_rate,
                    npts=sum(t.stats.npts for t in sub),
                    n_gaps=len([g for g in g_sub if g[6] > 0]),
                    gap_total_s=round(sum(g[6] for g in g_sub if g[6] > 0), 2),
                    status=status, sha256=sha256(fn)))
                print(f"{fn.name}: {log_rows[-1]['npts']} sampel, "
                      f"{log_rows[-1]['n_gaps']} gap")

    # katalog event pada jendela waktu (pembanding verifikasi nanti)
    try:
        lat, lon, _ = ref_event(client, t0, t1, args)
        cat = client.get_events(starttime=t0, endtime=t1, minmagnitude=args.cat_minmag,
                                latitude=lat, longitude=lon,
                                maxradius=kilometers2degrees(args.maxdist_km))
        cat.write(str(meta / "catalog_window.xml"), format="QUAKEML")
        print(f"Katalog: {len(cat)} event -> {meta/'catalog_window.xml'}")
    except Exception as e:
        print(f"Katalog gagal diunduh: {e}")

    with open(Path(args.outdir) / "download_log.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["file", "net", "sta", "loc", "cha", "start", "end",
                                          "sampling_rate", "npts", "n_gaps", "gap_total_s",
                                          "status", "sha256"])
        w.writeheader()
        w.writerows(log_rows)
    print("Log -> download_log.csv. PERIKSA isinya: selesai skrip != data lengkap.")


# ----------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("step", choices=["discover", "availability", "download"])
    p.add_argument("--server", default="AUSPASS", help="nama client ObsPy (AUSPASS/IRIS/GEONET)")
    p.add_argument("--start", required=True, help="awal jendela, UTC, mis. 2021-09-21")
    p.add_argument("--days", type=int, default=7)
    p.add_argument("--network", default="S1")
    p.add_argument("--channels", default="BH?,HH?,EH?")
    p.add_argument("--location", default="*")
    p.add_argument("--stations", default="", help="net.sta,net.sta (untuk download)")
    p.add_argument("--ref-lat", type=float)
    p.add_argument("--ref-lon", type=float)
    p.add_argument("--ref-minmag", type=float, default=5.0)
    p.add_argument("--cat-minmag", type=float, default=1.0)
    p.add_argument("--maxdist-km", type=float, default=300.0)
    p.add_argument("--max-stations", type=int, default=25)
    p.add_argument("--outdir", default="data")
    args = p.parse_args()
    {"discover": cmd_discover, "availability": cmd_availability,
     "download": cmd_download}[args.step](args)


if __name__ == "__main__":
    main()
