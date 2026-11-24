#!/usr/bin/env python3
"""
04_sweep_recall.py
Uji sensitivitas parameter STA/LTA yang menghitung seberapa banyak event katalog
berhasil ditemukan (recall) dan berapa trigger yang tidak cocok katalog.

Memakai fungsi dari 03_deteksi_sta_lta.py (harus satu folder dengan file ini),
jadi definisi trigger identik dengan script 03.

Input : data/processed/NET.STA.*.CHA.*.proc.mseed, data/meta/catalog_window.xml
Output: hasil_STA/results/parameter_sweep_recall.csv

Catatan istilah:
  "tidak cocok katalog" BUKAN otomatis false trigger. Bisa noise, bisa juga
  event kecil yang tidak ada di katalog. Statusnya baru jelas setelah diperiksa
  (stasiun kedua, bentuk waveform).

Contoh:
    python 04_sweep_recall.py --net S1 --sta AUSSC --cha HHZ --sta-lat -38.0977 --sta-lon 147.0825
"""
import argparse
import csv
import importlib.util
import itertools
from pathlib import Path

import numpy as np
from obspy import read
from obspy.signal.trigger import classic_sta_lta, trigger_onset
from obspy.taup import TauPyModel


def load_03():
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location("deteksi03", here / "03_deteksi_sta_lta.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def triggers(cfs, sta_s, lta_s, on, off, min_dur, min_sep):
    """Logika sama dengan detect() di script 03, tetapi CF dihitung sekali saja."""
    allev = []
    for tr, cf in cfs:
        sr = tr.stats.sampling_rate
        nsta, nlta = int(sta_s * sr), int(lta_s * sr)
        ratio = classic_sta_lta(cf, nsta, nlta)
        ratio[:nlta] = 0.0
        t0 = tr.stats.starttime
        for a, b in trigger_onset(ratio, on, off):
            if (b - a) / sr < min_dur:
                continue
            allev.append([t0 + a / sr, t0 + b / sr])
    allev.sort(key=lambda e: e[0])
    merged = []
    for e in allev:
        if merged and e[0] - merged[-1][1] < min_sep:
            merged[-1][1] = max(merged[-1][1], e[1])
        else:
            merged.append(e)
    return merged


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="data")
    ap.add_argument("--net", required=True)
    ap.add_argument("--sta", required=True)
    ap.add_argument("--cha", default="HHZ")
    ap.add_argument("--sta-lat", type=float)
    ap.add_argument("--sta-lon", type=float)
    ap.add_argument("--cf", choices=["envelope", "energy"], default="envelope")
    ap.add_argument("--min-dur", type=float, default=2.0)
    ap.add_argument("--min-sep", type=float, default=30.0)
    ap.add_argument("--match-tol", type=float, default=15.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    d = load_03()
    data = Path(args.data)
    out_csv = Path(args.out) if args.out else \
        Path(f"hasil_{args.sta}") / "results" / "parameter_sweep_recall.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    files = d.load_stream(data, args.net, args.sta, args.cha)
    traces = []
    for fn in files:
        traces.extend(read(fn))
    days = sum(t.stats.npts / t.stats.sampling_rate for t in traces) / 86400.0
    print(f"{len(traces)} trace, total {days:.2f} hari data")

    print("Menghitung characteristic function (sekali saja)...")
    cfs = [(tr, d.char_function(tr.data, tr.stats.sampling_rate, args.cf)) for tr in traces]

    slat, slon = d.station_coords(data, args.net, args.sta, args)
    model = TauPyModel("iasp91")
    catrows = d.catalog_table(data, slat, slon, model, args)
    p_ts = np.array([e["p_arr"].timestamp if e["p_arr"] is not None else np.nan for e in catrows])
    s_ts = np.array([e["s_arr"].timestamp if e["s_arr"] is not None else np.nan for e in catrows])
    n_cat = len(catrows)
    print(f"Katalog: {n_cat} event")

    grid = [(s, l, o, round(max(1.5, o * 0.4), 2), "grid")
            for s, l, o in itertools.product([0.5, 1.0, 2.0], [10, 30, 60], [3.0, 4.0, 5.0])]
    grid.append((1.0, 30, 4.0, 1.5, "default-script-03"))

    rows = []
    for sta_s, lta_s, on, off, tag in grid:
        trg = triggers(cfs, sta_s, lta_s, on, off, args.min_dur, args.min_sep)
        n = len(trg)
        if n:
            t = np.array([a.timestamp for a, _ in trg])
            hit = (np.abs(t[:, None] - p_ts[None, :]) <= args.match_tol) | \
                  (np.abs(t[:, None] - s_ts[None, :]) <= args.match_tol)
            n_match = int(hit.any(axis=1).sum())
            n_rec = int(hit.any(axis=0).sum())
        else:
            n_match = n_rec = 0
        rows.append(dict(
            sta_s=sta_s, lta_s=lta_s, on=on, off=off, catatan=tag,
            n_trigger=n, trigger_per_hari=round(n / days, 1),
            n_trigger_cocok_katalog=n_match,
            event_katalog_ditemukan=n_rec, dari_total_katalog=n_cat,
            recall_pct=round(100 * n_rec / n_cat, 1) if n_cat else 0,
            presisi_pct=round(100 * n_match / n, 1) if n else 0,
            n_tidak_cocok_katalog=n - n_match))
        print(f"STA={sta_s:<4} LTA={lta_s:<3} on={on:<3} off={off:<4}: trigger={n:<5} "
              f"event katalog={n_rec}/{n_cat} presisi={rows[-1]['presisi_pct']}%  [{tag}]")

    rows.sort(key=lambda r: (-r["recall_pct"], r["n_trigger"]))
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nTersimpan: {out_csv}")


if __name__ == "__main__":
    main()