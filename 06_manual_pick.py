#!/usr/bin/env python3
"""
Interactive manual P/S picking for the two-station review shortlist.

Run from the project directory:
  python 06_manual_pick.py --data data --verification-dir data\\verification

Controls in the plot:
  P or S : choose the phase to mark
  click   : place that phase on the waveform in the selected station panel
  U       : undo the last change for the current event
  C       : clear picks for the current event
  N / B   : next / previous event
  Q       : save and quit

Theoretical arrivals and STA/LTA triggers are references only. The output
contains only the times manually clicked by the user.
"""
import argparse
import csv
import glob
from datetime import timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from obspy import UTCDateTime, read


PICK_FIELDS = [
    "origin_utc", "magnitude", "station", "channel",
    "p_pick_utc", "p_uncertainty_s", "s_pick_utc", "s_uncertainty_s",
    "p_theoretical_utc", "s_theoretical_utc", "automatic_trigger_utc",
    "picker", "notes",
]


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def read_event_segments(files, start, end):
    segments = []
    for filename in files:
        headers = read(filename, headonly=True)
        if not any(tr.stats.endtime >= start and tr.stats.starttime <= end
                   for tr in headers):
            continue
        for tr in read(filename):
            if tr.stats.endtime >= start and tr.stats.starttime <= end:
                segment = tr.slice(start, end)
                if segment.stats.npts:
                    segments.append(segment)
    return segments


def load_processed_files(data, net, station, channel):
    pattern = str(data / "processed" /
                  f"{net}.{station}.*.{channel}.*.proc.mseed")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"Tidak ada waveform processed: {pattern}")
    return files


class ManualPicker:
    def __init__(self, events, stations, files_by_station, output_path,
                 channel, uncertainty_s, before_s, after_s):
        self.events = events
        self.stations = stations
        self.files_by_station = files_by_station
        self.output_path = output_path
        self.channel = channel
        self.uncertainty_s = uncertainty_s
        self.before_s = before_s
        self.after_s = after_s
        self.index = 0
        self.mode = None
        self.undo_stack = []
        self.existing_rows = self._load_existing()
        self.rows = self._make_rows()
        self.fig, self.axes = plt.subplots(
            2, 1, figsize=(14, 8), sharex=True)
        self.help_text = self.fig.text(
            0.5, 0.015, "", ha="center", va="bottom", fontsize=9)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self.fig.canvas.mpl_connect("button_press_event", self._on_click)
        self.fig.canvas.mpl_connect("close_event", self._on_close)
        self._save()
        self._draw()

    def _load_existing(self):
        if not self.output_path.exists():
            return {}
        saved = {}
        with open(self.output_path, newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                key = (row.get("origin_utc", ""), row.get("station", ""))
                if all(key):
                    saved[key] = {field: row.get(field, "") for field in PICK_FIELDS}
        return saved

    def _make_rows(self):
        rows = {}
        for event in self.events:
            for station in self.stations:
                key = (event["origin"], station)
                saved = self.existing_rows.get(key, {})
                rows[key] = dict(
                    origin_utc=event["origin"],
                    magnitude=event["magnitude"],
                    station=station,
                    channel=self.channel,
                    p_pick_utc=saved.get("p_pick_utc", ""),
                    p_uncertainty_s=saved.get("p_uncertainty_s", ""),
                    s_pick_utc=saved.get("s_pick_utc", ""),
                    s_uncertainty_s=saved.get("s_uncertainty_s", ""),
                    p_theoretical_utc=event.get(f"{station}_p_theoretical", ""),
                    s_theoretical_utc=event.get(f"{station}_s_theoretical", ""),
                    automatic_trigger_utc=self._automatic_trigger(event, station),
                    picker=saved.get("picker", "manual"),
                    notes=saved.get("notes", ""),
                )
        return rows

    @staticmethod
    def _automatic_trigger(event, station):
        if station == event.get("candidate_station"):
            return event.get("candidate_trigger", "")
        if station == event.get("compare_station"):
            return event.get("compare_trigger", "")
        return ""

    def _current_event(self):
        return self.events[self.index]

    def _current_row(self, station):
        return self.rows[(self._current_event()["origin"], station)]

    def _save(self):
        all_rows = dict(self.existing_rows)
        all_rows.update(self.rows)
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=PICK_FIELDS)
            writer.writeheader()
            writer.writerows(sorted(
                all_rows.values(),
                key=lambda row: (row["origin_utc"], row["station"])))
        self.existing_rows = all_rows

    @staticmethod
    def _as_time(value):
        return UTCDateTime(value) if value else None

    def _draw(self):
        event = self._current_event()
        phase_times = []
        for station in self.stations:
            for phase in ("p", "s"):
                value = self._as_time(event.get(f"{station}_{phase}_theoretical", ""))
                if value:
                    phase_times.append(value)
        if not phase_times:
            raise SystemExit(f"Tidak ada prediksi P/S untuk {event['origin']}")
        start = min(phase_times) - self.before_s
        end = max(phase_times) + self.after_s

        for ax, station in zip(self.axes, self.stations):
            ax.clear()
            segments = read_event_segments(
                self.files_by_station[station], start, end)
            for i, tr in enumerate(segments):
                ax.plot(tr.times("matplotlib"), tr.data, color="black",
                        linewidth=0.65,
                        label="Waveform processed (m/s)" if i == 0 else None)

            row = self._current_row(station)
            for phase, color in (("p", "tab:blue"), ("s", "tab:green")):
                theoretical = self._as_time(row[f"{phase}_theoretical_utc"])
                if theoretical:
                    ax.axvline(
                        theoretical.matplotlib_date, color=color,
                        linestyle=":", linewidth=1.4,
                        label=f"{phase.upper()} teoretis")
                picked = self._as_time(row[f"{phase}_pick_utc"])
                if picked:
                    ax.axvline(
                        picked.matplotlib_date, color=color,
                        linestyle="-", linewidth=2.0,
                        label=f"{phase.upper()} manual")
            automatic = self._as_time(row["automatic_trigger_utc"])
            if automatic:
                ax.axvline(
                    automatic.matplotlib_date, color="tab:red",
                    linestyle="--", linewidth=1.3,
                    label="Trigger STA/LTA otomatis")
            if not segments:
                ax.text(
                    0.5, 0.5, "Waveform tidak tersedia pada jendela ini",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_ylabel(f"{station}\n{self.channel} (m/s)")
            ax.grid(True, alpha=0.25)
            handles, labels = ax.get_legend_handles_labels()
            unique = dict(zip(labels, handles))
            if unique:
                ax.legend(unique.values(), unique.keys(),
                          loc="upper right", fontsize=8)

        magnitude = event["magnitude"]
        mag_text = f"M{float(magnitude):.2f}" if magnitude else "M?"
        self.fig.suptitle(
            f"Manual pick {self.index + 1}/{len(self.events)} | "
            f"origin katalog {event['origin']} | {mag_text}\n"
            "Garis putus-putus = prediksi; garis merah = trigger otomatis; "
            "garis solid warna P/S = pick manual")
        self.axes[-1].xaxis.set_major_formatter(
            mdates.DateFormatter("%H:%M:%S"))
        self.axes[-1].set_xlabel("Waktu UTC")
        p_counts = [
            f"{station} P:{bool(self._current_row(station)['p_pick_utc'])} "
            f"S:{bool(self._current_row(station)['s_pick_utc'])}"
            for station in self.stations
        ]
        mode = self.mode or "belum dipilih"
        self.help_text.set_text(
            f"Mode: {mode} | {', '.join(p_counts)} | "
            "Tekan P/S lalu klik waveform; U undo; C hapus event; "
            "N berikutnya; B sebelumnya; Q simpan & keluar. "
            f"Ketidakpastian default: +/-{self.uncertainty_s:g} s")
        self.fig.autofmt_xdate()
        self.fig.tight_layout(rect=(0, 0.06, 1, 0.92))
        self.fig.canvas.draw_idle()

    def _on_click(self, event):
        if event.inaxes not in self.axes or event.xdata is None:
            return
        if self.mode not in ("P", "S"):
            self.help_text.set_text(
                "Pilih fase dulu dengan tombol P atau S, lalu klik waveform.")
            self.fig.canvas.draw_idle()
            return
        station = self.stations[list(self.axes).index(event.inaxes)]
        row = self._current_row(station)
        field = f"{self.mode.lower()}_pick_utc"
        uncertainty = f"{self.mode.lower()}_uncertainty_s"
        previous = (row[field], row[uncertainty])
        self.undo_stack.append((station, field, uncertainty, previous))
        pick_time = UTCDateTime(
            mdates.num2date(event.xdata, tz=timezone.utc).timestamp())
        row[field] = str(pick_time)
        row[uncertainty] = f"{self.uncertainty_s:g}"
        self._save()
        self._draw()

    def _on_key(self, event):
        key = (event.key or "").lower()
        if key in ("p", "s"):
            self.mode = key.upper()
            self._draw()
        elif key == "u":
            if self.undo_stack:
                station, field, uncertainty, previous = self.undo_stack.pop()
                row = self._current_row(station)
                row[field], row[uncertainty] = previous
                self._save()
                self._draw()
        elif key == "c":
            for station in self.stations:
                row = self._current_row(station)
                for phase in ("p", "s"):
                    field = f"{phase}_pick_utc"
                    uncertainty = f"{phase}_uncertainty_s"
                    previous = (row[field], row[uncertainty])
                    self.undo_stack.append((station, field, uncertainty, previous))
                    row[field], row[uncertainty] = "", ""
            self._save()
            self._draw()
        elif key in ("n", "right"):
            if self.index + 1 < len(self.events):
                self.index += 1
                self.mode = None
                self.undo_stack.clear()
                self._draw()
        elif key in ("b", "left"):
            if self.index > 0:
                self.index -= 1
                self.mode = None
                self.undo_stack.clear()
                self._draw()
        elif key == "q":
            self._save()
            plt.close(self.fig)

    def _on_close(self, _event):
        self._save()

    def run(self):
        print(f"CSV manual pick: {self.output_path}")
        print("Klik hanya setelah memilih P atau S. Tutup jendela atau tekan Q untuk selesai.")
        plt.show(block=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data")
    parser.add_argument("--verification-dir", default="data/verification")
    parser.add_argument("--net", default="S1")
    parser.add_argument("--station-1", default="AUSSC")
    parser.add_argument("--station-2", default="AUKSC")
    parser.add_argument("--channel", default="HHZ")
    parser.add_argument("--uncertainty-s", type=float, default=0.5,
                        help="ketidakpastian pick manual yang dicatat untuk setiap klik")
    parser.add_argument("--before-s", type=float, default=10.0)
    parser.add_argument("--after-s", type=float, default=25.0)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()
    if args.uncertainty_s < 0 or args.before_s < 0 or args.after_s <= 0:
        parser.error("Ketidakpastian/before harus >= 0 dan after harus > 0")

    verification_dir = Path(args.verification_dir)
    summary_path = verification_dir / "verification_summary.csv"
    shortlist_path = verification_dir / "review_shortlist.csv"
    if not summary_path.exists() or not shortlist_path.exists():
        raise SystemExit(
            "Jalankan 05_verifikasi_dua_stasiun.py terlebih dahulu "
            "untuk membuat verification_summary.csv dan review_shortlist.csv.")
    summary = {row["origin"]: row for row in read_csv(summary_path)}
    shortlist = read_csv(shortlist_path)
    events = []
    for item in shortlist:
        origin = item["origin"]
        if origin not in summary:
            raise SystemExit(f"Event shortlist tidak ada di {summary_path}: {origin}")
        events.append(summary[origin])
    if not events:
        raise SystemExit(f"Shortlist kosong: {shortlist_path}")

    stations = [args.station_1, args.station_2]
    files_by_station = {
        station: load_processed_files(
            Path(args.data), args.net, station, args.channel)
        for station in stations
    }
    output = Path(args.output) if args.output else verification_dir / "manual_picks.csv"
    picker = ManualPicker(
        events, stations, files_by_station, output, args.channel,
        args.uncertainty_s, args.before_s, args.after_s)
    picker.run()


if __name__ == "__main__":
    main()
