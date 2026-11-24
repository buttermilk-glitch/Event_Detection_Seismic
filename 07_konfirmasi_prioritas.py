#!/usr/bin/env python3
"""Create a visual-review report for the five highest-priority events.

The report is scoped to the supplied catalog events and the two HHZ stations.
It preserves the original candidate files and writes new outputs under
data/verification/confirmed by default.
"""
import argparse
import csv
import glob
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from obspy import UTCDateTime, read


CONFIRMATION_BASIS = (
    "Visual review: coherent P-like onset near the predicted P arrival at "
    "both stations, followed by later energy consistent with S; automatic "
    "STA/LTA triggers also match catalog phase times at both stations. "
    "This confirms the catalog-event association for this review, not a "
    "manual phase pick or an independent event relocation."
)


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as file:
        return list(csv.DictReader(file))


def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_origin(value):
    return UTCDateTime(value)


def load_station_files(data_dir, network, station, channel):
    pattern = str(
        data_dir / "processed" /
        f"{network}.{station}.*.{channel}.*.proc.mseed"
    )
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"Tidak ada waveform processed: {pattern}")
    return files


def read_event_segments(files, start, end):
    segments = []
    for filename in files:
        for trace in read(filename):
            if trace.stats.endtime >= start and trace.stats.starttime <= end:
                segment = trace.slice(start, end)
                if segment.stats.npts:
                    segments.append(segment)
    return segments


def make_plot(index, row, stations, files_by_station, output_path,
              before_s, after_s):
    phase_times = []
    for station in stations:
        for phase in ("p", "s"):
            value = row.get(f"{station}_{phase}_theoretical", "")
            if value:
                phase_times.append(parse_origin(value))
    if not phase_times:
        raise SystemExit(f"Tidak ada prediksi P/S untuk {row['origin']}")

    start = min(phase_times) - before_s
    end = max(phase_times) + after_s
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    fig.patch.set_facecolor("#fbfcfe")

    for ax, station in zip(axes, stations):
        ax.set_facecolor("#fbfcfe")
        traces = read_event_segments(files_by_station[station], start, end)
        for trace_index, trace in enumerate(traces):
            ax.plot(
                trace.times("matplotlib"), trace.data,
                color="#263238", linewidth=0.65,
                label="Waveform processed (m/s)" if trace_index == 0 else None,
                zorder=2,
            )

        for phase, color in (("p", "#1976d2"), ("s", "#2e7d32")):
            value = row.get(f"{station}_{phase}_theoretical", "")
            if value:
                arrival = parse_origin(value)
                ax.axvline(
                    arrival.matplotlib_date, color=color, linestyle=":",
                    linewidth=1.8, label=f"{phase.upper()} teoretis iasp91",
                    zorder=3,
                )

        trigger = row.get(f"{station}_trigger", "")
        if trigger:
            trigger_time = parse_origin(trigger)
            ax.axvline(
                trigger_time.matplotlib_date, color="#d84315",
                linestyle="--", linewidth=1.7,
                label="Trigger STA/LTA otomatis", zorder=4,
            )

        if not traces:
            ax.text(
                0.5, 0.5, "Waveform tidak tersedia pada jendela ini",
                ha="center", va="center", transform=ax.transAxes,
            )
        ax.set_ylabel(f"{station}\n{row['channel']} (m/s)")
        ax.grid(True, color="#cfd8dc", alpha=0.6, linewidth=0.7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        if unique:
            ax.legend(
                unique.values(), unique.keys(), loc="upper right",
                fontsize=8, framealpha=0.92,
            )

        match_phase = row.get(f"{station}_trigger_phase", "")
        delta = row.get(f"{station}_trigger_delta_s", "")
        if match_phase and delta != "":
            ax.text(
                0.01, 0.96,
                f"Trigger cocok {match_phase} | selisih {float(delta):+.2f} s",
                transform=ax.transAxes, ha="left", va="top", fontsize=9,
                color="#9a3412",
                bbox={"facecolor": "white", "edgecolor": "#fed7aa",
                      "alpha": 0.9, "boxstyle": "round,pad=0.3"},
            )

    origin = parse_origin(row["origin"])
    magnitude = float(row["magnitude"]) if row.get("magnitude") else None
    mag_text = f"M{magnitude:.2f}" if magnitude is not None else "M?"
    fig.suptitle(
        f"Event {index:02d} | {origin} | {mag_text} | "
        "CONFIRMED - tinjauan visual dua stasiun",
        fontsize=14, fontweight="bold", y=0.98,
    )
    fig.text(
        0.5, 0.935,
        "P/S adalah prediksi model; garis merah adalah trigger otomatis, "
        "bukan manual pick.",
        ha="center", va="top", fontsize=9, color="#455a64",
    )
    axes[-1].xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M:%S", tz=timezone.utc)
    )
    axes[-1].set_xlabel("Waktu UTC")
    fig.autofmt_xdate()
    fig.tight_layout(rect=(0.02, 0.03, 0.98, 0.91))
    fig.savefig(output_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def make_timeline(rows, output_path):
    ordered = sorted(rows, key=lambda row: int(row["review_rank"]))
    fig, ax = plt.subplots(figsize=(13, 5.5))
    fig.patch.set_facecolor("#fbfcfe")
    ax.set_facecolor("#fbfcfe")
    colors = ["#1565c0", "#00897b", "#7b1fa2", "#ef6c00", "#c62828"]
    dates = [
        datetime.fromisoformat(row["origin_utc"].replace("Z", "+00:00"))
        for row in ordered
    ]
    y_positions = list(range(len(ordered), 0, -1))
    for row, event_date, y, color in zip(ordered, dates, y_positions, colors):
        magnitude = float(row["magnitude"]) if row["magnitude"] else 0
        x = mdates.date2num(event_date)
        ax.scatter(
            [x], [y], s=95 + magnitude * 12, color=color,
            edgecolor="white", linewidth=1.5, zorder=3,
        )
        ax.hlines(y, x - 0.035, x, color=color, alpha=0.35, linewidth=1)
        ax.annotate(
            f"{event_date:%Y-%m-%d %H:%M:%S} UTC   |   M{magnitude:.2f}",
            (x, y), xytext=(10, 0), textcoords="offset points",
            va="center", fontsize=10, color="#263238",
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels([f"Prioritas #{row['review_rank']}" for row in ordered])
    ax.set_ylim(0.4, len(ordered) + 0.6)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1, tz=timezone.utc))
    ax.xaxis.set_major_formatter(
        mdates.DateFormatter("%d %b", tz=timezone.utc)
    )
    ax.set_xlabel("Waktu origin katalog (UTC)")
    ax.set_title(
        "Lima event prioritas - waktu origin katalog",
        loc="left", fontsize=15, fontweight="bold", pad=18,
    )
    ax.text(
        0, 1.02,
        "Titik menunjukkan waktu origin, bukan waktu tiba gelombang di stasiun.",
        transform=ax.transAxes, fontsize=9, color="#546e7a",
    )
    ax.grid(axis="x", color="#cfd8dc", alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data")
    parser.add_argument("--verification-dir", default="data/verification")
    parser.add_argument("--outdir", default="data/verification/confirmed")
    parser.add_argument("--net", default="S1")
    parser.add_argument("--station-1", default="AUSSC")
    parser.add_argument("--station-2", default="AUKSC")
    parser.add_argument("--channel", default="HHZ")
    parser.add_argument("--before-s", type=float, default=10.0)
    parser.add_argument("--after-s", type=float, default=25.0)
    args = parser.parse_args()
    if args.before_s < 0 or args.after_s <= 0:
        parser.error("--before-s harus >= 0 dan --after-s harus > 0")

    verification_dir = Path(args.verification_dir)
    summary_path = verification_dir / "verification_summary.csv"
    shortlist_path = verification_dir / "review_shortlist.csv"
    if not summary_path.exists() or not shortlist_path.exists():
        raise SystemExit(
            "Jalankan 05_verifikasi_dua_stasiun.py terlebih dahulu."
        )
    summary = {row["origin"]: row for row in read_csv(summary_path)}
    shortlist = read_csv(shortlist_path)
    if not shortlist:
        raise SystemExit(f"Shortlist kosong: {shortlist_path}")

    stations = [args.station_1, args.station_2]
    normalized_rows = []
    for item in shortlist:
        origin = item["origin"]
        if origin not in summary:
            raise SystemExit(f"Event tidak ada di verification_summary.csv: {origin}")
        source = summary[origin]
        if source.get("trigger_match_both_stations") != "yes":
            raise SystemExit(
                f"Event {origin} tidak memiliki trigger cocok di kedua stasiun; "
                "tidak akan diberi status CONFIRMED."
            )
        row = dict(source)
        row["review_rank"] = item["review_rank"]
        row["channel"] = args.channel
        row["status"] = "CONFIRMED"
        row["confirmation_basis"] = CONFIRMATION_BASIS
        row["AUSSC_trigger"] = source.get("candidate_trigger", "")
        row["AUSSC_trigger_phase"] = source.get("candidate_match_phase", "")
        row["AUSSC_trigger_delta_s"] = source.get("candidate_dt_s", "")
        row["AUKSC_trigger"] = source.get("compare_trigger", "")
        row["AUKSC_trigger_phase"] = source.get("compare_phase", "")
        row["AUKSC_trigger_delta_s"] = source.get("compare_dt_s", "")
        normalized_rows.append(row)

    outdir = Path(args.outdir)
    event_dir = outdir / "events"
    event_dir.mkdir(parents=True, exist_ok=True)
    data_dir = Path(args.data)
    files_by_station = {
        station: load_station_files(
            data_dir, args.net, station, args.channel
        )
        for station in stations
    }

    confirmed_fields = [
        "review_rank", "origin", "magnitude", "status",
        "candidate_station", "candidate_trigger", "candidate_match_phase",
        "candidate_dt_s", "compare_station", "compare_trigger",
        "compare_phase", "compare_dt_s", "trigger_match_both_stations",
        "AUSSC_p_theoretical", "AUSSC_s_theoretical",
        "AUKSC_p_theoretical", "AUKSC_s_theoretical",
        "confirmation_basis",
    ]
    confirmed_rows = [
        {field: row.get(field, "") for field in confirmed_fields}
        for row in normalized_rows
    ]
    write_csv(outdir / "confirmed_events.csv", confirmed_rows, confirmed_fields)

    time_fields = [
        "review_rank", "origin_utc", "magnitude",
        "AUSSC_p_theoretical_utc", "AUSSC_trigger_utc",
        "AUSSC_trigger_phase", "AUSSC_trigger_delta_s",
        "AUSSC_s_theoretical_utc",
        "AUKSC_p_theoretical_utc", "AUKSC_trigger_utc",
        "AUKSC_trigger_phase", "AUKSC_trigger_delta_s",
        "AUKSC_s_theoretical_utc",
    ]
    time_rows = []
    for row in normalized_rows:
        time_rows.append({
            "review_rank": row["review_rank"],
            "origin_utc": row["origin"],
            "magnitude": row["magnitude"],
            "AUSSC_p_theoretical_utc": row.get("AUSSC_p_theoretical", ""),
            "AUSSC_trigger_utc": row["AUSSC_trigger"],
            "AUSSC_trigger_phase": row["AUSSC_trigger_phase"],
            "AUSSC_trigger_delta_s": row["AUSSC_trigger_delta_s"],
            "AUSSC_s_theoretical_utc": row.get("AUSSC_s_theoretical", ""),
            "AUKSC_p_theoretical_utc": row.get("AUKSC_p_theoretical", ""),
            "AUKSC_trigger_utc": row["AUKSC_trigger"],
            "AUKSC_trigger_phase": row["AUKSC_trigger_phase"],
            "AUKSC_trigger_delta_s": row["AUKSC_trigger_delta_s"],
            "AUKSC_s_theoretical_utc": row.get("AUKSC_s_theoretical", ""),
        })
    write_csv(outdir / "event_times.csv", time_rows, time_fields)

    for index, row in enumerate(normalized_rows, 1):
        origin_key = row["origin"][:10]
        make_plot(
            index, row, stations, files_by_station,
            event_dir / f"confirmed_{index:02d}_{origin_key}.png",
            args.before_s, args.after_s,
        )
    make_timeline(
        [
            {
                "review_rank": row["review_rank"],
                "origin_utc": row["origin"],
                "magnitude": row["magnitude"],
            }
            for row in normalized_rows
        ],
        outdir / "event_timeline_utc.png",
    )
    print(f"Event confirmed (visual review): {outdir / 'confirmed_events.csv'}")
    print(f"Waktu origin/arrival/trigger: {outdir / 'event_times.csv'}")
    print(f"Timeline UTC: {outdir / 'event_timeline_utc.png'}")
    print(f"Plot detail: {event_dir}")


if __name__ == "__main__":
    main()
