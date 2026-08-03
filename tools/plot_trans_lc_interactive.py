import argparse
from glob import glob

import matplotlib.pylab as plt
import numpy as np
from matplotlib.widgets import Slider


plt.rc("font", family="sans-serif")
plt.rc("font", serif="Helvetica Neue")
plt.rc("text", usetex="false")


DEFAULT_PERIOD_DAYS = 4.055
STATE_LABELS = {
    0: "out",
    1: "partial",
    2: "full",
}
QUANTITY_LABELS = {
    "flux": "Normalised flux",
    "depth": "Transit depth",
    "depth_opaque": "Opaque-body transit depth",
    "depth_atm": "Atmosphere-only transit depth",
    "east": "East-limb atmosphere transit depth",
    "west": "West-limb atmosphere transit depth",
}
QUANTITY_KEYS = ("flux", "depth", "depth_opaque", "depth_atm", "east", "west")
SPECTRUM_KEYS = QUANTITY_KEYS + ("all-components",)


def centered_phase(phase):
    return ((phase + 0.5) % 1.0) - 0.5


def read_transit_file(fname):
    """Read one Transit_*.txt light-curve output file."""
    header = np.loadtxt(fname, max_rows=1)
    data = np.atleast_2d(np.loadtxt(fname, skiprows=1))

    if header.size < 10:
        raise ValueError(f"{fname}: expected 10 header values, got {header.size}")
    if data.shape[1] < 5:
        raise ValueError(
            f"{fname}: expected wl, depth, depth_atm, depth_atm_east, depth_atm_west"
        )
    if data.shape[1] >= 6:
        depth_opaque = data[:, 5]
    else:
        depth_opaque = data[:, 1] - data[:, 2]

    meta = {
        "fname": fname,
        "n_wl": int(header[0]),
        "rp": header[1],
        "r_top": header[2],
        "viewphi": header[3],
        "viewtheta": header[4],
        "phase": header[5],
        "state": int(header[6]),
        "xstar": header[7],
        "ystar": header[8],
        "zstar": header[9],
    }

    return meta, data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], depth_opaque


def load_transit_run(pattern, period_days):
    rows = [read_transit_file(fname) for fname in sorted(glob(pattern))]
    if not rows:
        raise FileNotFoundError(f"No transit files matched {pattern!r}")

    ref_wl = rows[0][1]
    for meta, wl, *_ in rows[1:]:
        if len(wl) != len(ref_wl) or not np.allclose(wl, ref_wl):
            raise ValueError(f"{meta['fname']}: wavelength grid differs from first file")

    rows.sort(key=lambda item: centered_phase(item[0]["phase"]))

    phase = np.array([centered_phase(row[0]["phase"]) for row in rows])
    depth = np.array([row[2] for row in rows])
    return {
        "meta": [row[0] for row in rows],
        "wl": ref_wl,
        "phase": phase,
        "time_hours": phase * period_days * 24.0,
        "state": np.array([row[0]["state"] for row in rows]),
        "depth": depth,
        "depth_atm": np.array([row[3] for row in rows]),
        "east": np.array([row[4] for row in rows]),
        "west": np.array([row[5] for row in rows]),
        "depth_opaque": np.array([row[6] for row in rows]),
        "flux": 1.0 - depth,
    }


def wavelength_slice(n_wl, center_idx, half_width):
    center_idx = int(np.clip(round(center_idx), 0, n_wl - 1))
    half_width = int(np.clip(round(half_width), 0, n_wl - 1))
    start = max(0, center_idx - half_width)
    stop = min(n_wl, center_idx + half_width + 1)
    return center_idx, half_width, slice(start, stop)


def band_label(wl, band_slice):
    band = wl[band_slice]
    if len(band) == 1:
        return f"{band[0]:.6g} um"
    return f"{band[0]:.6g}-{band[-1]:.6g} um ({len(band)} bins)"


class TransitExplorer:
    def __init__(self, run, quantity, spectrum_quantity):
        self.run = run
        self.quantity = quantity
        self.spectrum_quantity = spectrum_quantity
        self.selected_phase_idx = len(run["time_hours"]) // 2

        self.fig, (self.ax_lc, self.ax_spec) = plt.subplots(2, 1, figsize=(10, 8))
        self.fig.subplots_adjust(left=0.10, right=0.96, top=0.92, bottom=0.20, hspace=0.35)

        self.lc_line, = self.ax_lc.plot([], [], "ko-", ms=4, lw=1.2)
        self.selected_point, = self.ax_lc.plot([], [], "o", ms=9, mfc="none", mec="firebrick", mew=1.8)

        wl_slider_ax = self.fig.add_axes([0.16, 0.105, 0.76, 0.03])
        width_slider_ax = self.fig.add_axes([0.16, 0.055, 0.76, 0.03])
        n_wl = len(run["wl"])
        slider_max = max(n_wl - 1, 1)
        self.wl_slider = Slider(
            wl_slider_ax,
            "Wavelength bin",
            0,
            slider_max,
            valinit=0,
            valstep=1,
        )
        self.width_slider = Slider(
            width_slider_ax,
            "Half-width bins",
            0,
            slider_max,
            valinit=0,
            valstep=1,
        )

        self.wl_slider.on_changed(self.update_light_curve)
        self.width_slider.on_changed(self.update_light_curve)
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)

        self.update_light_curve(None)

    def current_band(self):
        return wavelength_slice(
            len(self.run["wl"]),
            self.wl_slider.val,
            self.width_slider.val,
        )

    def update_light_curve(self, _event):
        _center_idx, _half_width, band_slice = self.current_band()
        y = np.mean(self.run[self.quantity][:, band_slice], axis=1)
        self.lc_line.set_data(self.run["time_hours"], y)
        self.selected_point.set_data(
            [self.run["time_hours"][self.selected_phase_idx]],
            [y[self.selected_phase_idx]],
        )

        self.ax_lc.relim()
        self.ax_lc.autoscale_view()
        self.ax_lc.set_xlabel("Time from mid-transit [hours]")
        self.ax_lc.set_ylabel(QUANTITY_LABELS[self.quantity])
        self.ax_lc.set_title(f"{QUANTITY_LABELS[self.quantity]} at {band_label(self.run['wl'], band_slice)}")
        self.update_spectrum()
        self.fig.canvas.draw_idle()

    def update_spectrum(self):
        self.ax_spec.clear()
        wl = self.run["wl"]
        idx = self.selected_phase_idx
        _center_idx, _half_width, band_slice = self.current_band()

        if self.spectrum_quantity == "all-components":
            self.ax_spec.plot(wl, self.run["depth"][idx], label="Total depth", color="black")
            self.ax_spec.plot(wl, self.run["depth_opaque"][idx], label="Opaque body", color="darkorange")
            self.ax_spec.plot(wl, self.run["depth_atm"][idx], label="Atmosphere", color="0.45")
            self.ax_spec.plot(wl, self.run["east"][idx], label="East limb", color="firebrick")
            self.ax_spec.plot(wl, self.run["west"][idx], label="West limb", color="royalblue")
            self.ax_spec.set_ylabel("Transit depth")
            self.ax_spec.legend(fontsize=9)
        else:
            self.ax_spec.plot(wl, self.run[self.spectrum_quantity][idx], color="black")
            self.ax_spec.set_ylabel(QUANTITY_LABELS[self.spectrum_quantity])

        self.ax_spec.axvspan(wl[band_slice][0], wl[band_slice][-1], color="gold", alpha=0.18)
        if len(wl) > 1:
            self.ax_spec.set_xscale("log")
        self.ax_spec.set_xlabel(r"$\lambda$ [$\mu$m]")

        meta = self.run["meta"][idx]
        state = STATE_LABELS.get(self.run["state"][idx], f"state {self.run['state'][idx]}")
        self.ax_spec.set_title(
            f"{meta['fname']} | phase {self.run['phase'][idx]:+.6f} | "
            f"{self.run['time_hours'][idx]:+.3f} h | {state}"
        )
        self.ax_spec.relim()
        self.ax_spec.autoscale_view()

    def on_click(self, event):
        if event.inaxes is not self.ax_lc or event.xdata is None:
            return

        self.selected_phase_idx = int(np.argmin(np.abs(self.run["time_hours"] - event.xdata)))
        self.update_light_curve(None)

    def show(self):
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Interactively explore 3D_sph_trans_lc Transit_*.txt output files."
    )
    parser.add_argument(
        "--pattern",
        default="Transit_*.txt",
        help="Glob pattern for transit light-curve output files.",
    )
    parser.add_argument(
        "--period-days",
        type=float,
        default=DEFAULT_PERIOD_DAYS,
        help=f"Orbital period in days. Default is {DEFAULT_PERIOD_DAYS} for WASP-39b.",
    )
    parser.add_argument(
        "--quantity",
        choices=QUANTITY_KEYS,
        default="flux",
        help="Quantity shown in the light-curve panel.",
    )
    parser.add_argument(
        "--spectrum-quantity",
        choices=SPECTRUM_KEYS,
        default="depth",
        help="Quantity shown in the clicked-phase spectrum panel.",
    )
    args = parser.parse_args()

    run = load_transit_run(args.pattern, args.period_days)
    TransitExplorer(run, args.quantity, args.spectrum_quantity).show()


if __name__ == "__main__":
    main()
