#!/usr/bin/env python3
"""
02_qc_dan_preprocessing.py
Quality control + preprocessing data harian hasil 01_cek_dan_unduh_auspass.py.

Input  : data/raw/NET.STA.LOC.CHA.YYYY-MM-DD.mseed   (TIDAK diubah)
         data/meta/NET.STA.response.xml               (opsional, untuk remove_response)
Output : data/processed/*.mseed   (float32, satuan m/s jika ada respons, selain itu counts)
         data/qc/qc_report.csv    (gap, sampling rate, clipping, RMS, % cakupan per file)
         data/figures/qc_*.png

Contoh:
    python 02_qc_dan_preprocessing.py --net S1 --sta AUGRF --cha BHZ --freqmin 1 --freqmax 10

Alasan singkat tiap langkah (rujukan: NMSOP-2 Ch.16; ObsPy docs):
  - split() pada gap: filter/detrend tidak boleh melintasi gap (artefak).
  - demean+linear detrend + taper: mencegah transien filter pada tepi segmen.
  - remove_response + pre_filt: mengubah counts -> m/s dan menstabilkan dekonvolusi
    (pre_filt meredam frekuensi di luar band respons instrumen).
  - bandpass zero-phase: meningkatkan SNR event lokal/regional (default 1-10 Hz);
    zero-phase tidak menggeser waktu tiba (penting untuk picking).
  - Tanpa resampling kecuali diminta (--resample), agar tidak mengubah isi spektral.
"""
import argparse
import csv
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import Stream, UTCDateTime, read, read_inventory

DAY = 86400.0


def qc_file(fn):
    """QC satu file harian mentah. Mengembalikan (dict, Stream)."""
    st = read(fn)
    gaps = st.get_gaps()
    real_gaps = [g for g in gaps if g[6] > 0]
    overlaps = [g for g in gaps if g[6] < 0]
    srs = sorted({float(t.stats.sampling_rate) for t in st})
    npts = sum(t.stats.npts for t in st)
    sr = srs[0]
    t0 = min(t.stats.starttime for t in st)
    t1 = max(t.stats.endtime for t in st)
    data = np.concatenate([t.data.astype(float) for t in st])
    amax = np.max(np.abs(data)) if data.size else 0
    n_at_max = int(np.sum(np.abs(data) == amax)) if amax > 0 else 0
    rec = dict(
        file=Path(fn).name, start=t0, end=t1, sampling_rates=";".join(map(str, srs)),
        npts=npts, coverage_pct=round(100 * npts / (DAY * sr), 3),
        n_gaps=len(real_gaps), gap_total_s=round(sum(g[6] for g in real_gaps), 2),
        n_overlaps=len(overlaps), rms_counts=round(float(np.sqrt(np.mean(data ** 2))), 2),
        max_abs_counts=amax, n_samples_at_max=n_at_max,
        n_nonfinite=int(np.sum(~np.isfinite(data))),
        flat_fraction=round(float(np.mean(np.diff(data) == 0)), 4) if data.size > 1 else 1,
    )
    flags = []
    if len(srs) > 1:
        flags.append("SAMPLING_RATE_TIDAK_KONSISTEN")
    if rec["coverage_pct"] < 95:
        flags.append("COVERAGE<95%")
    if n_at_max > 10 and n_at_max / max(npts, 1) > 1e-4:
        flags.append("KEMUNGKINAN_CLIPPING")
    if rec["flat_fraction"] > 0.2:
        flags.append("BANYAK_SAMPEL_DATAR")
    if rec["n_nonfinite"]:
        flags.append("NaN/Inf")
    rec["flags"] = "|".join(flags) or "OK"
    return rec, st


def preprocess(st, inv, args):
    st = st.copy()
    st = st.split()  # pecah di gap, hasil: trace tanpa masked array
    st = Stream([t for t in st if t.stats.npts > 20 * t.stats.sampling_rate])
    for tr in st:
        tr.data = tr.data.astype(np.float64)
    st.detrend("demean")
    st.detrend("linear")
    st.taper(max_percentage=0.02, max_length=10.0, type="hann")
    unit = "counts"
    if inv is not None and not args.no_response:
        nyq = st[0].stats.sampling_rate / 2.0
        pre = (0.02, 0.05, 0.85 * nyq, 0.95 * nyq)
        st.remove_response(inventory=inv, output="VEL", pre_filt=pre, zero_mean=True,
                           taper=False)
        unit = "m/s"
    if args.resample:
        st.resample(args.resample, no_filter=False)
    nyq = st[0].stats.sampling_rate / 2.0
    fmax = min(args.freqmax, 0.9 * nyq)
    st.filter("bandpass", freqmin=args.freqmin, freqmax=fmax, corners=4, zerophase=True)
    for tr in st:
        tr.data = tr.data.astype(np.float32)
        tr.stats.processing_unit = unit
    return st, unit


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="data")
    ap.add_argument("--net", required=True)
    ap.add_argument("--sta", required=True)
    ap.add_argument("--cha", default="BHZ")
    ap.add_argument("--freqmin", type=float, default=1.0)
    ap.add_argument("--freqmax", type=float, default=10.0)
    ap.add_argument("--resample", type=float, default=None)
    ap.add_argument("--no-response", action="store_true")
    ap.add_argument("--demo-start", default=None,
                    help="UTC mulai plot 10 menit mentah vs proses (default: tengah hari ke-1)")
    args = ap.parse_args()

    out = Path(args.outdir)
    (out / "processed").mkdir(parents=True, exist_ok=True)
    (out / "qc").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(out / "raw" / f"{args.net}.{args.sta}.*.{args.cha}.*.mseed")))
    if not files:
        raise SystemExit("Tidak ada file raw yang cocok. Cek --net/--sta/--cha dan folder data/raw.")

    inv = None
    xml = out / "meta" / f"{args.net}.{args.sta}.response.xml"
    if xml.exists():
        inv = read_inventory(str(xml))
    else:
        print("PERINGATAN: StationXML tidak ditemukan -> tanpa koreksi respons (satuan counts).")

    report, minmax = [], []
    first_raw = first_proc = None
    for fn in files:
        rec, st = qc_file(fn)
        proc, unit = preprocess(st, inv, args)
        rec["processed_unit"] = unit
        report.append(rec)
        o = out / "processed" / Path(fn).name.replace(".mseed", ".proc.mseed")
        proc.write(str(o), format="MSEED", encoding="FLOAT32")
        print(f"{rec['file']}: cov={rec['coverage_pct']}% gaps={rec['n_gaps']} "
              f"flags={rec['flags']} -> {o.name}")
        # min/max per menit untuk plot 7 hari
        for tr in proc:
            sr = tr.stats.sampling_rate
            n = int(60 * sr)
            m = len(tr.data) // n
            if m == 0:
                continue
            a = tr.data[: m * n].reshape(m, n)
            tt = tr.stats.starttime.matplotlib_date + (np.arange(m) * 60 + 30) / DAY
            minmax.append((tt, a.min(1), a.max(1)))
        if first_raw is None:
            first_raw, first_proc = st, proc

    with open(out / "qc" / "qc_report.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(report[0].keys()))
        w.writeheader()
        w.writerows(report)
    print(f"QC report -> {out/'qc'/'qc_report.csv'}")

    # ---- plot 1: raw vs processed (10 menit) + PSD sebelum/sesudah
    t_mid = UTCDateTime(args.demo_start) if args.demo_start else \
        min(t.stats.starttime for t in first_raw) + 12 * 3600
    raw10 = first_raw.copy().trim(t_mid, t_mid + 600)
    pro10 = first_proc.copy().trim(t_mid, t_mid + 600)
    if raw10 and pro10:
        fig, ax = plt.subplots(3, 1, figsize=(11, 8))
        r, p = raw10[0], pro10[0]
        ax[0].plot(r.times(), r.data, lw=0.4)
        ax[0].set_title(f"Mentah (counts) {r.id} {t_mid}")
        ax[1].plot(p.times(), p.data, lw=0.4, color="C1")
        ax[1].set_title(f"Setelah preprocessing ({p.stats.processing_unit}), "
                        f"bandpass {args.freqmin}-{args.freqmax} Hz")
        ax[1].set_xlabel("detik sejak awal jendela")
        ax[2].psd(r.data - np.mean(r.data), NFFT=2048, Fs=r.stats.sampling_rate,
                  label="mentah (counts)", color="C0")
        ax[2].legend()
        ax[2].set_title("PSD jendela 10 menit (mentah); bandingkan dengan hasil proses di ax lain")
        fig.tight_layout()
        fig.savefig(out / "figures" / "qc_raw_vs_processed.png", dpi=130)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.psd(r.data - np.mean(r.data), NFFT=2048, Fs=r.stats.sampling_rate,
               label="mentah (counts)")
        ax.psd(p.data, NFFT=2048, Fs=p.stats.sampling_rate, label=f"proses ({p.stats.processing_unit})")
        ax.legend()
        ax.set_title("Noise sebelum vs sesudah preprocessing (10 menit)")
        fig.tight_layout()
        fig.savefig(out / "figures" / "qc_noise_psd.png", dpi=130)
        plt.close(fig)

    # ---- plot 2: continuous N hari (min/max per menit)
    if minmax:
        fig, ax = plt.subplots(figsize=(14, 4))
        for tt, lo, hi in minmax:
            ax.fill_between(tt, lo, hi, color="k", lw=0)
        ax.xaxis_date()
        ax.set_title(f"Continuous waveform {args.net}.{args.sta}..{args.cha} (min/max per menit, proses)")
        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(out / "figures" / "qc_continuous_7d.png", dpi=130)
        plt.close(fig)
    print("Plot -> data/figures/")


if __name__ == "__main__":
    main()
