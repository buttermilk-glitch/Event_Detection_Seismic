#!/usr/bin/env python3
"""
Bandingkan trigger kandidat dari satu stasiun dengan waveform stasiun kedua.

Input:
  data/results/candidates.csv
  data/meta/catalog_window.xml
  data/meta/NET.STA.response.xml
  data/processed/NET.STA.*.CHA.*.proc.mseed

Output:
  data/verification/verification_summary.csv
  data/verification/review_shortlist.csv
  data/verification/events/*.png

Example:
  python 05_verifikasi_dua_stasiun.py --data data --net S1 ^
      --candidate-station AUSSC --compare-station AUKSC --cha HHZ

Garis P/S pada plot adalah waktu teoretis dari katalog dan model iasp91.
Trigger otomatis yang berdekatan bukan manual pick atau konfirmasi gempa.
"""
import argparse
import csv
import glob
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from obspy import read, read_events, read_inventory
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
from obspy.taup import TauPyModel


def load_detector():
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location(
        "deteksi_sta_lta", here / "03_deteksi_sta_lta.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("Tidak bisa memuat 03_deteksi_sta_lta.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def read_station_streams(data, net, station, channel):
    pattern = str(data / "processed" /
                  f"{net}.{station}.*.{channel}.*.proc.mseed")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"Tidak ada waveform processed untuk {net}.{station}.{channel}: {pattern}")
    return files


def detect_station(files, detector, args):
    found = []
    skipped = []
    min_samples = None
    for fn in files:
        for tr in read(fn):
            required = int(args.lta_len * tr.stats.sampling_rate)
            min_samples = required if min_samples is None else min_samples
            if len(tr.data) < required:
                skipped.append((tr.id, float(tr.stats.endtime - tr.stats.starttime)))
                continue
            for start, end, peak in detector.detect(
                    tr, args.sta_len, args.lta_len, args.on, args.off,
                    args.cf, args.min_dur, args.min_sep):
                found.append(dict(start=start, end=end, peak=peak))
    found.sort(key=lambda item: item["start"])
    return found, skipped, min_samples


def station_coordinates(data, net, station):
    path = data / "meta" / f"{net}.{station}.response.xml"
    if not path.exists():
        raise SystemExit(f"StationXML tidak ditemukan: {path}")
    inv = read_inventory(str(path))
    selected = inv.select(network=net, station=station)
    if not selected or not selected[0]:
        raise SystemExit(f"Koordinat stasiun {net}.{station} tidak ada di {path}")
    sta = selected[0][0]
    return sta.latitude, sta.longitude


def phase_times(event, station_lat, station_lon, model):
    origin = event.preferred_origin() or event.origins[0]
    magnitude = event.preferred_magnitude()
    if magnitude is None and event.magnitudes:
        magnitude = event.magnitudes[0]
    distance_m = gps2dist_azimuth(
        origin.latitude, origin.longitude, station_lat, station_lon)[0]
    distance_deg = kilometers2degrees(distance_m / 1000.0)
    depth_km = (origin.depth or 10000.0) / 1000.0
    p = model.get_travel_times(depth_km, distance_deg, ["p", "P", "Pn"])
    s = model.get_travel_times(depth_km, distance_deg, ["s", "S", "Sn"])
    return dict(
        origin=origin.time,
        magnitude=magnitude.mag if magnitude else "",
        distance_km=distance_m / 1000.0,
        p=origin.time + min(a.time for a in p) if p else None,
        s=origin.time + min(a.time for a in s) if s else None,
    )


def closest_match(triggers, phases, tolerance):
    matches = []
    for phase in ("P", "S"):
        predicted = phases[phase.lower()]
        if predicted is None:
            continue
        for trigger in triggers:
            delta = float(trigger["start"] - predicted)
            if abs(delta) <= tolerance:
                matches.append((abs(delta), phase, delta, trigger))
    if not matches:
        return None
    _, phase, delta, trigger = min(matches, key=lambda item: item[0])
    return dict(phase=phase, delta=delta, trigger=trigger)


def event_traces(files, start, end):
    selected = []
    for fn in files:
        headers = read(fn, headonly=True)
        if not any(tr.stats.endtime >= start and tr.stats.starttime <= end
                   for tr in headers):
            continue
        for tr in read(fn):
            if tr.stats.endtime >= start and tr.stats.starttime <= end:
                segment = tr.slice(start, end)
                if segment.stats.npts:
                    selected.append(segment)
    return selected


def draw_event_plot(index, origin_key, event_info, stations, files_by_station,
                    triggers_by_station, matches, outdir):
    phases = [event_info[station] for station in stations]
    predicted = [t for phase_set in phases for t in (phase_set["p"], phase_set["s"])
                 if t is not None]
    if not predicted:
        return
    start = min(predicted) - 10
    end = max(predicted) + 25
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
    for ax, station in zip(axes, stations):
        segments = event_traces(files_by_station[station], start, end)
        for i, tr in enumerate(segments):
            ax.plot(tr.times("matplotlib"), tr.data, linewidth=0.55,
                    color="black", label="Waveform processed (m/s)" if i == 0 else None)
        station_phase = event_info[station]
        for key, label, color in (
                ("p", "P teoretis", "tab:blue"),
                ("s", "S teoretis", "tab:green")):
            arrival = station_phase[key]
            if arrival is not None:
                ax.axvline(arrival.matplotlib_date, color=color, linestyle=":",
                           linewidth=1.3, label=label)
        for trigger in triggers_by_station[station]:
            if start <= trigger["start"] <= end:
                ax.axvline(trigger["start"].matplotlib_date, color="0.65",
                           linestyle="--", alpha=0.55, linewidth=0.8,
                           label="Trigger STA/LTA" if "Trigger STA/LTA" not in
                           ax.get_legend_handles_labels()[1] else None)
        match = matches[station]
        if match is not None:
            ax.axvline(match["trigger"]["start"].matplotlib_date, color="tab:red",
                       linestyle="-", linewidth=1.5,
                       label=f"Trigger cocok waktu {match['phase']} (otomatis)")
        if not segments:
            ax.text(0.5, 0.5, "Waveform tidak tersedia pada jendela ini",
                    ha="center", va="center", transform=ax.transAxes)
        ax.set_ylabel(f"{station}\n(m/s)")
        ax.grid(True, alpha=0.2)
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        if unique:
            ax.legend(unique.values(), unique.keys(), loc="upper right", fontsize=8)

    magnitude = event_info[stations[0]]["magnitude"]
    mag_text = f"M{magnitude:.2f}" if magnitude != "" else "magnitudo tidak tersedia"
    origin = event_info[stations[0]]["origin"]
    fig.suptitle(
        f"Kandidat {index:02d} | origin katalog {origin} | {mag_text}\n"
        "P/S = prediksi iasp91; trigger = STA/LTA otomatis; bukan pick manual")
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    axes[-1].set_xlabel("Waktu UTC")
    fig.autofmt_xdate()
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(outdir / f"event_{index:02d}_{origin_key[:10]}.png", dpi=140)
    plt.close(fig)


def write_rows(path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data")
    parser.add_argument("--net", default="S1")
    parser.add_argument("--candidate-station", default="AUSSC",
                        help="stasiun yang menghasilkan candidates.csv")
    parser.add_argument("--compare-station", default="AUKSC")
    parser.add_argument("--cha", default="HHZ")
    parser.add_argument("--sta-len", type=float, default=1.0)
    parser.add_argument("--lta-len", type=float, default=30.0)
    parser.add_argument("--on", type=float, default=4.0)
    parser.add_argument("--off", type=float, default=1.5)
    parser.add_argument("--cf", choices=["envelope", "energy"], default="envelope")
    parser.add_argument("--min-dur", type=float, default=2.0)
    parser.add_argument("--min-sep", type=float, default=30.0)
    parser.add_argument("--match-tol", type=float, default=15.0)
    parser.add_argument("--top", type=int, default=5)
    parser.add_argument("--outdir", default="data/verification")
    args = parser.parse_args()
    if args.sta_len <= 0 or args.lta_len <= args.sta_len:
        parser.error("--lta-len harus lebih besar dari --sta-len, dan keduanya positif")
    if args.top < 1:
        parser.error("--top harus minimal 1")

    data = Path(args.data)
    outdir = Path(args.outdir)
    plot_dir = outdir / "events"
    plot_dir.mkdir(parents=True, exist_ok=True)
    detector = load_detector()
    stations = [args.candidate_station, args.compare_station]

    candidate_path = data / "results" / "candidates.csv"
    catalog_path = data / "meta" / "catalog_window.xml"
    candidates = [row for row in read_csv(candidate_path)
                  if row.get("status") == "PROBABLE" and row.get("match_origin")]
    by_origin = {}
    for row in candidates:
        by_origin.setdefault(row["match_origin"], row)
    if not by_origin:
        raise SystemExit(f"Tidak ada kandidat PROBABLE di {candidate_path}")

    catalog = read_events(str(catalog_path))
    events = {}
    for event in catalog:
        origin = event.preferred_origin() or event.origins[0]
        events[str(origin.time)] = event
    missing = sorted(set(by_origin) - set(events))
    if missing:
        raise SystemExit(f"Origin kandidat tidak ditemukan di katalog: {missing}")

    files_by_station = {}
    triggers_by_station = {}
    for station in stations:
        files = read_station_streams(data, args.net, station, args.cha)
        triggers, skipped, min_samples = detect_station(files, detector, args)
        files_by_station[station] = files
        triggers_by_station[station] = triggers
        print(f"{args.net}.{station}.{args.cha}: {len(triggers)} trigger otomatis")
        for trace_id, seconds in skipped:
            print(f"  Dilewati: segmen {trace_id} hanya {seconds:.2f} dtk "
                  f"(LTA perlu {args.lta_len:g} dtk / {min_samples} sampel)")

    coordinates = {
        station: station_coordinates(data, args.net, station)
        for station in stations
    }
    model = TauPyModel("iasp91")
    summary = []
    plot_data = {}
    for origin_key, source_candidate in by_origin.items():
        event = events[origin_key]
        event_info = {
            station: phase_times(event, *coordinates[station], model)
            for station in stations
        }
        matches = {
            station: closest_match(
                triggers_by_station[station], event_info[station], args.match_tol)
            for station in stations
        }
        reference_match = matches[args.candidate_station]
        other_match = matches[args.compare_station]
        detected_at_both = reference_match is not None and other_match is not None
        row = dict(
            origin=origin_key,
            magnitude=event_info[args.candidate_station]["magnitude"],
            candidate_station=args.candidate_station,
            candidate_trigger=source_candidate.get("t_on", ""),
            candidate_snr=source_candidate.get("snr", ""),
            candidate_match_phase=source_candidate.get("match_phase", ""),
            candidate_dt_s=source_candidate.get("dt_s", ""),
            compare_station=args.compare_station,
            compare_trigger=(str(other_match["trigger"]["start"]) if other_match else ""),
            compare_phase=(other_match["phase"] if other_match else ""),
            compare_dt_s=(round(other_match["delta"], 2) if other_match else ""),
            trigger_match_both_stations="yes" if detected_at_both else "no",
            interpretation=("two-station automatic timing match; review waveform"
                            if detected_at_both else
                            "no automatic match at both stations; inspect waveform"),
        )
        for station in stations:
            row[f"{station}_distance_km"] = round(
                event_info[station]["distance_km"], 1)
            row[f"{station}_p_theoretical"] = (
                str(event_info[station]["p"]) if event_info[station]["p"] else "")
            row[f"{station}_s_theoretical"] = (
                str(event_info[station]["s"]) if event_info[station]["s"] else "")
        summary.append(row)
        plot_data[origin_key] = (event_info, matches)

    summary.sort(key=lambda row: (
        row["trigger_match_both_stations"] != "yes",
        -float(row["candidate_snr"] or 0)))
    write_rows(outdir / "verification_summary.csv", summary)

    shortlist = []
    for index, row in enumerate(summary[:args.top], 1):
        shortlist.append(dict(
            review_rank=index,
            origin=row["origin"],
            magnitude=row["magnitude"],
            trigger_match_both_stations=row["trigger_match_both_stations"],
            candidate_snr=row["candidate_snr"],
            reason="prioritas pemeriksaan visual; bukan status confirmed",
        ))
    write_rows(outdir / "review_shortlist.csv", shortlist)

    for index, row in enumerate(summary, 1):
        event_info, matches = plot_data[row["origin"]]
        draw_event_plot(
            index, row["origin"], event_info, stations, files_by_station,
            triggers_by_station, matches, plot_dir)
    print(f"Ringkasan: {outdir / 'verification_summary.csv'}")
    print(f"Prioritas review: {outdir / 'review_shortlist.csv'}")
    print(f"Plot per event: {plot_dir}")
    print("Trigger cocok di dua stasiun bukan manual pick atau konfirmasi otomatis.")


if __name__ == "__main__":
    main()
