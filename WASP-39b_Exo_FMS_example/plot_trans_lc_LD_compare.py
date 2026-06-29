import argparse
from glob import glob
from pathlib import Path

import matplotlib.pylab as plt
import numpy as np


plt.rc("font", family="sans-serif")
plt.rc("font", serif="Helvetica Neue")
plt.rc("text", usetex="false")


DEFAULT_PERIOD_DAYS = 4.0552941
STATE_LABELS = {
    0: "out",
    1: "partial",
    2: "full",
}


def centered_phase(phase):
    return ((phase + 0.5) % 1.0) - 0.5


def read_transit_file(fname):
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


def load_run(directory):
    pattern = str(Path(directory) / "Transit_*.txt")
    rows = [read_transit_file(fname) for fname in sorted(glob(pattern))]
    if not rows:
        raise FileNotFoundError(f"No transit files matched {pattern!r}")

    rows.sort(key=lambda item: centered_phase(item[0]["phase"]))

    wl = rows[0][1]
    phases = np.array([centered_phase(row[0]["phase"]) for row in rows])
    states = np.array([row[0]["state"] for row in rows])
    depth = np.array([row[2] for row in rows])
    depth_atm = np.array([row[3] for row in rows])
    east = np.array([row[4] for row in rows])
    west = np.array([row[5] for row in rows])
    depth_opaque = np.array([row[6] for row in rows])

    for row in rows[1:]:
        if len(row[1]) != len(wl) or not np.allclose(row[1], wl):
            raise ValueError(f"{row[0]['fname']}: wavelength grid differs from first file")

    return {
        "directory": directory,
        "wl": wl,
        "phase": phases,
        "state": states,
        "depth": depth,
        "depth_opaque": depth_opaque,
        "depth_atm": depth_atm,
        "east": east,
        "west": west,
        "flux": 1.0 - depth,
    }


def validate_pair(off_ld, on_ld):
    if len(off_ld["phase"]) != len(on_ld["phase"]):
        raise ValueError("LD-off and LD-on runs have different numbers of transit phases")
    if len(off_ld["wl"]) != len(on_ld["wl"]) or not np.allclose(off_ld["wl"], on_ld["wl"]):
        raise ValueError("LD-off and LD-on runs have different wavelength grids")
    if not np.allclose(off_ld["phase"], on_ld["phase"]):
        raise ValueError("LD-off and LD-on runs have different phase sampling")


def plot_compare(off_ld, on_ld, period_days):
    validate_pair(off_ld, on_ld)

    wl = off_ld["wl"]
    phase = off_ld["phase"]
    time_hours = phase * period_days * 24.0
    state = off_ld["state"]

    off_depth = off_ld["depth"]
    on_depth = on_ld["depth"]
    off_flux = off_ld["flux"]
    on_flux = on_ld["flux"]

    depth_diff = on_depth - off_depth
    flux_diff = on_flux - off_flux

    fig, ax = plt.subplots()
    off_band = np.mean(off_flux, axis=1)
    on_band = np.mean(on_flux, axis=1)
    ax.plot(time_hours, off_band, "ko-", label="LD off")
    ax.plot(time_hours, on_band, "o-", color="firebrick", label="LD on")
    for s in np.unique(state):
        idx = state == s
        state_name = STATE_LABELS.get(s, f"state {s}")
        ax.scatter(time_hours[idx], on_band[idx], s=70, facecolors="none", label=state_name)
    ax.set_xlabel("Time from mid-transit [hours]", fontsize=14)
    ax.set_ylabel("Normalised flux", fontsize=14)
    ax.legend(fontsize=9)
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

    fig, ax = plt.subplots()
    ax.plot(time_hours, np.mean(on_flux - off_flux, axis=1), "ko-")
    ax.axhline(0.0, color="0.5", lw=1)
    ax.set_xlabel("Time from mid-transit [hours]", fontsize=14)
    ax.set_ylabel("LD on - LD off normalised flux", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

    fig, ax = plt.subplots()
    sample_idx = np.linspace(0, len(wl) - 1, min(5, len(wl)), dtype=int)
    for idx in sample_idx:
        ax.plot(time_hours, off_flux[:, idx], ls="--", alpha=0.6, label=f"off {wl[idx]:.3g} um")
        ax.plot(time_hours, on_flux[:, idx], alpha=0.8, label=f"on {wl[idx]:.3g} um")
    ax.set_xlabel("Time from mid-transit [hours]", fontsize=14)
    ax.set_ylabel("Normalised flux", fontsize=14)
    ax.legend(fontsize=8, ncol=2)
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(wl, time_hours, flux_diff, shading="auto", cmap="RdBu_r")
    ax.set_xlabel(r"$\lambda$ [$\mu$m]", fontsize=14)
    ax.set_ylabel("Time from mid-transit [hours]", fontsize=14)
    ax.set_xscale("log")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("LD on - LD off normalised flux")
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(wl, time_hours, depth_diff, shading="auto", cmap="RdBu_r")
    ax.set_xlabel(r"$\lambda$ [$\mu$m]", fontsize=14)
    ax.set_ylabel("Time from mid-transit [hours]", fontsize=14)
    ax.set_xscale("log")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("LD on - LD off transit depth")
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Compare 3D_sph_trans_lc outputs with limb darkening on and off."
    )
    parser.add_argument(
        "--off-dir",
        default="off_LD",
        help="Directory containing no-limb-darkening Transit_*.txt files.",
    )
    parser.add_argument(
        "--on-dir",
        default="on_LD",
        help="Directory containing limb-darkened Transit_*.txt files.",
    )
    parser.add_argument(
        "--period-days",
        type=float,
        default=DEFAULT_PERIOD_DAYS,
        help="Orbital period in days. Default is 4.0552941 for WASP-39b.",
    )
    args = parser.parse_args()

    off_ld = load_run(args.off_dir)
    on_ld = load_run(args.on_dir)
    plot_compare(off_ld, on_ld, args.period_days)


if __name__ == "__main__":
    main()
