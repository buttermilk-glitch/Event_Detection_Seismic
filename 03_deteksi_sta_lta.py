#!/usr/bin/env python3
"""
03_deteksi_sta_lta.py
Deteksi kandidat event (STA/LTA) pada data proses, uji sensitivitas parameter,
dan pencocokan dengan katalog (waktu tiba P/S teoretis, model iasp91).

Input  : data/processed/NET.STA.*.CHA.*.proc.mseed
         data/meta/catalog_window.xml   (dari skrip 01)
         koordinat stasiun: dari data/meta/NET.STA.response.xml atau --sta-lat/--sta-lon
Output : data/results/parameter_sweep.csv
         data/results/candidates.csv         (semua trigger + status)
         data/results/catalog_vs_detection.csv
         data/figures/event_*.png, timeline_events.png

Contoh:
    python 03_deteksi_sta_lta.py --net S1 --sta AUGRF --cha BHZ \
        --sta-len 1.0 --lta-len 30 --on 4.0 --off 1.5 --utc-offset 10

PENTING tentang status:
  Skrip hanya memberi status "PROBABLE" (cocok katalog + 1 stasiun) atau
  "UNVERIFIED" (tidak cocok). "CONFIRMED" menuntut bukti tambahan (stasiun
  lain / picking P-S konsisten) yang Anda periksa sendiri; skrip tidak
  menetapkannya otomatis.

Keterbatasan: deteksi per file harian; LTA butuh pemanasan (lta-len detik awal
tiap hari diabaikan), sehingga event dalam lta-len detik pertama tiap hari
bisa terlewat. Dokumentasikan di slide.
"""
import argparse
import csv
import glob
import itertools
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from obspy import Stream, UTCDateTime, read, read_events, read_inventory
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
from obspy.signal.trigger import classic_sta_lta, trigger_onset
from obspy.taup import TauPyModel
from scipy.fft import next_fast_len
from scipy.signal import hilbert


def load_stream(out, net, sta, cha):
    files = sorted(glob.glob(str(out / "processed" / f"{net}.{sta}.*.{cha}.*.proc.mseed")))
    if not files:
        raise SystemExit("Tidak ada file processed. Jalankan skrip 02 dulu.")
    return files


def char_function(data, sr, kind):
    if kind == "energy":
        return data.astype(np.float64) ** 2
    n = len(data)
    nf = next_fast_len(n)
    return np.abs(hilbert(data.astype(np.float64), N=nf))[:n]  # envelope (Earle & Shearer 1994)


def detect(tr, sta_s, lta_s, on, off, cf_kind, min_dur, min_sep):
    sr = tr.stats.sampling_rate
    cf = char_function(tr.data, sr, cf_kind)
    nsta, nlta = int(sta_s * sr), int(lta_s * sr)
    ratio = classic_sta_lta(cf, nsta, nlta)
    ratio[: nlta] = 0.0  # abaikan masa pemanasan LTA
    trig = trigger_onset(ratio, on, off)
    t0 = tr.stats.starttime
    ev = []
    for a, b in trig:
        dur = (b - a) / sr
        if dur < min_dur:
            continue
        ev.append([t0 + a / sr, t0 + b / sr, float(ratio[a:b].max())])
    ev.sort(key=lambda e: e[0])
    merged = []
    for e in ev:  # gabungkan trigger yang terlalu dekat (minimum separation)
        if merged and e[0] - merged[-1][1] < min_sep:
            merged[-1][1] = max(merged[-1][1], e[1])
            merged[-1][2] = max(merged[-1][2], e[2])
        else:
            merged.append(e)
    return merged


def station_coords(out, net, sta, args):
    if args.sta_lat is not None:
        return args.sta_lat, args.sta_lon
    xml = out / "meta" / f"{net}.{sta}.response.xml"
    inv = read_inventory(str(xml))
    s = inv.select(network=net, station=sta)[0][0]
    return s.latitude, s.longitude


def catalog_table(out, slat, slon, model, args):
    cat = read_events(str(out / "meta" / "catalog_window.xml"))
    rows = []
    for ev in cat:
        o = ev.preferred_origin() or ev.origins[0]
        m = ev.preferred_magnitude() or (ev.magnitudes[0] if ev.magnitudes else None)
        d_m = gps2dist_azimuth(o.latitude, o.longitude, slat, slon)[0]
        dist_deg = kilometers2degrees(d_m / 1000.0)
        depth = (o.depth or 10000.0) / 1000.0
        tp = model.get_travel_times(depth, dist_deg, ["p", "P", "Pn"])
        ts = model.get_travel_times(depth, dist_deg, ["s", "S", "Sn"])
        rows.append(dict(origin=o.time, lat=o.latitude, lon=o.longitude, depth_km=round(depth, 1),
                         mag=m.mag if m else None, mag_type=m.magnitude_type if m else "",
                         dist_km=round(d_m / 1000.0, 1),
                         p_arr=o.time + min(a.time for a in tp) if tp else None,
                         s_arr=o.time + min(a.time for a in ts) if ts else None))
    return rows


def snr_and_freq(tr, t_on, noise_s=60.0, sig_s=20.0):
    noise = tr.slice(t_on - noise_s - 5, t_on - 5).data
    sig = tr.slice(t_on, t_on + sig_s).data
    if len(noise) < 10 or len(sig) < 10:
        return np.nan, np.nan
    snr = float(np.max(np.abs(sig)) / (np.sqrt(np.mean(noise.astype(float) ** 2)) + 1e-30))
    spec = np.abs(np.fft.rfft(sig * np.hanning(len(sig))))
    f = np.fft.rfftfreq(len(sig), 1 / tr.stats.sampling_rate)
    return snr, float(f[np.argmax(spec[1:]) + 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="data")
    ap.add_argument("--net", required=True)
    ap.add_argument("--sta", required=True)
    ap.add_argument("--cha", default="BHZ")
    ap.add_argument("--sta-lat", type=float)
    ap.add_argument("--sta-lon", type=float)
    ap.add_argument("--cf", choices=["envelope", "energy"], default="envelope")
    ap.add_argument("--sta-len", type=float, default=1.0)
    ap.add_argument("--lta-len", type=float, default=30.0)
    ap.add_argument("--on", type=float, default=4.0)
    ap.add_argument("--off", type=float, default=1.5)
    ap.add_argument("--min-dur", type=float, default=2.0)
    ap.add_argument("--min-sep", type=float, default=30.0)
    ap.add_argument("--match-tol", type=float, default=15.0,
                    help="toleransi (detik) trigger vs waktu tiba P/S teoretis")
    ap.add_argument("--utc-offset", type=float, default=0.0,
                    help="selisih jam lokal terhadap UTC (untuk menandai jam kerja)")
    ap.add_argument("--max-plots", type=int, default=12)
    args = ap.parse_args()

    out = Path(args.outdir)
    for d in ("results", "figures"):
        (out / d).mkdir(parents=True, exist_ok=True)
    files = load_stream(out, args.net, args.sta, args.cha)
    traces = []
    for fn in files:
        traces.extend(read(fn))

    # ---------- 1) uji sensitivitas parameter
    sweep = []
    for sta_s, lta_s, on in itertools.product([0.5, 1.0, 2.0], [10, 30, 60], [3.0, 4.0, 5.0]):
        n = 0
        for tr in traces:
            n += len(detect(tr, sta_s, lta_s, on, max(1.5, on * 0.4), args.cf,
                            args.min_dur, args.min_sep))
        sweep.append(dict(sta_s=sta_s, lta_s=lta_s, on=on, off=round(max(1.5, on * 0.4), 2),
                          n_triggers=n))
        print(f"STA={sta_s} LTA={lta_s} on={on}: {n} trigger")
    with open(out / "results" / "parameter_sweep.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sweep[0].keys()))
        w.writeheader()
        w.writerows(sweep)

    # ---------- 2) deteksi dengan parameter terpilih
    cands = []
    for tr in traces:
        for t_on, t_off, peak in detect(tr, args.sta_len, args.lta_len, args.on, args.off,
                                        args.cf, args.min_dur, args.min_sep):
            snr, fdom = snr_and_freq(tr, t_on)
            lh = (t_on.hour + t_on.minute / 60 + args.utc_offset) % 24
            cands.append(dict(t_on=t_on, t_off=t_off, dur_s=round(t_off - t_on, 1),
                              peak_ratio=round(peak, 2), snr=round(snr, 1),
                              f_dom_hz=round(fdom, 2), local_hour=round(lh, 1),
                              work_hours=(8 <= lh <= 17)))
    print(f"{len(cands)} kandidat dengan parameter terpilih")

    # ---------- 3) pencocokan katalog
    slat, slon = station_coords(out, args.net, args.sta, args)
    model = TauPyModel("iasp91")
    catrows = catalog_table(out, slat, slon, model, args)
    for c in cands:
        c.update(match_origin="", match_mag="", match_dist_km="", match_phase="",
                 dt_s="", status="UNVERIFIED")
        best = None
        for e in catrows:
            for ph, key in (("P", "p_arr"), ("S", "s_arr")):
                if e[key] is None:
                    continue
                dt = c["t_on"] - e[key]
                if abs(dt) <= args.match_tol and (best is None or abs(dt) < abs(best[0])):
                    best = (dt, e, ph)
        if best:
            dt, e, ph = best
            c.update(match_origin=str(e["origin"]), match_mag=e["mag"],
                     match_dist_km=e["dist_km"], match_phase=ph, dt_s=round(dt, 1),
                     status="PROBABLE")
        elif c["work_hours"] and c["f_dom_hz"] > 5:
            c["status"] = "UNVERIFIED (jam kerja, frekuensi tinggi: curiga antropogenik)"

    with open(out / "results" / "candidates.csv", "w", newline="") as f:
        if cands:
            w = csv.DictWriter(f, fieldnames=list(cands[0].keys()))
            w.writeheader()
            w.writerows(cands)

    detected = set()
    cv = []
    for e in catrows:
        hit = [c for c in cands if c["match_origin"] == str(e["origin"])]
        if hit:
            detected.add(str(e["origin"]))
        cv.append(dict(origin=e["origin"], mag=e["mag"], mag_type=e["mag_type"],
                       dist_km=e["dist_km"], depth_km=e["depth_km"], p_arr=e["p_arr"],
                       s_arr=e["s_arr"], detected="ya" if hit else "tidak",
                       peak_ratio=max([h["peak_ratio"] for h in hit], default=""),
                       snr=max([h["snr"] for h in hit], default="")))
    with open(out / "results" / "catalog_vs_detection.csv", "w", newline="") as f:
        if cv:
            w = csv.DictWriter(f, fieldnames=list(cv[0].keys()))
            w.writeheader()
            w.writerows(cv)
    n_match = sum(1 for c in cands if c["status"] == "PROBABLE")
    print(f"Katalog: {len(catrows)} event | terdeteksi: {len(detected)} | "
          f"trigger PROBABLE: {n_match} | UNVERIFIED: {len(cands) - n_match}")

    # ---------- 4) plot per event (prioritaskan yang cocok katalog)
    full = Stream(traces)
    order = sorted(cands, key=lambda c: (c["status"] != "PROBABLE", -c["snr"]
                                         if np.isfinite(c["snr"]) else 0))
    for i, c in enumerate(order[: args.max_plots], 1):
        seg = full.slice(c["t_on"] - 30, c["t_on"] + 120).merge(fill_value=0)
        if not seg:
            continue
        tr = seg[0]
        fig, ax = plt.subplots(figsize=(11, 3.5))
        ax.plot(tr.times("matplotlib"), tr.data, lw=0.5, color="k")
        ax.axvline(c["t_on"].matplotlib_date, color="r", ls="--", label="trigger on")
        for e in catrows:
            if str(e["origin"]) == c["match_origin"]:
                for lab, key, col in (("P teoretis", "p_arr", "b"), ("S teoretis", "s_arr", "g")):
                    ax.axvline(e[key].matplotlib_date, color=col, ls=":", label=lab)
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
        ax.set_title(f"Kandidat {i}: {c['t_on']} | rasio={c['peak_ratio']} SNR={c['snr']} | {c['status']}",
                     fontsize=9)
        ax.legend(fontsize=7, loc="upper right")
        fig.tight_layout()
        fig.savefig(out / "figures" / f"event_{i:02d}.png", dpi=120)
        plt.close(fig)

    # ---------- 5) timeline katalog vs deteksi
    fig, ax = plt.subplots(figsize=(13, 3))
    if cands:
        ax.scatter([c["t_on"].matplotlib_date for c in cands], [1] * len(cands), marker="|",
                   s=200, color="C1", label="trigger STA/LTA")
    if catrows:
        ax.scatter([e["origin"].matplotlib_date for e in catrows], [2] * len(catrows),
                   marker="v", color="C0", label="katalog (origin)")
    ax.set_yticks([1, 2])
    ax.set_yticklabels(["deteksi", "katalog"])
    ax.xaxis_date()
    fig.autofmt_xdate()
    ax.legend(fontsize=8)
    ax.set_title("Distribusi waktu: katalog vs deteksi")
    fig.tight_layout()
    fig.savefig(out / "figures" / "timeline_events.png", dpi=130)
    plt.close(fig)
    print("Selesai. Periksa SETIAP plot event secara visual sebelum mengklaim status.")


if __name__ == "__main__":
    main()
