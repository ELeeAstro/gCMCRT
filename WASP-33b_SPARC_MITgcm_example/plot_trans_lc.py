import argparse
from glob import glob

import matplotlib.pylab as plt
import numpy as np


plt.rc("font", family="sans-serif")
plt.rc("font", serif="Helvetica Neue")
plt.rc("text", usetex="false")


RSUN = 6.957e10
DEFAULT_RS = 1.444 * RSUN
DEFAULT_PERIOD_DAYS = 1.21987
STATE_LABELS = {
    0: "out",
    1: "partial",
    2: "full",
}


def centered_phase(phase):
    return ((phase + 0.5) % 1.0) - 0.5


def read_transit_file(fname):
    """Read one Transit_*.txt file.

    Current output columns are finished, dimensionless depths:
        wl, depth, depth_atm, depth_atm_east, depth_atm_west, depth_opaque
    Older five-column files are still accepted; depth_opaque is inferred as
    depth - depth_atm.
    """
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


def load_transit_files(pattern):
    rows = []
    for fname in sorted(glob(pattern)):
        rows.append(read_transit_file(fname))

    if not rows:
        raise FileNotFoundError(f"No transit files matched {pattern!r}")

    ref_wl = rows[0][1]
    for meta, wl, *_ in rows[1:]:
        if len(wl) != len(ref_wl) or not np.allclose(wl, ref_wl):
            raise ValueError(f"{meta['fname']}: wavelength grid differs from first file")

    rows.sort(key=lambda item: centered_phase(item[0]["phase"]))
    return rows


def plot_transit_results(rows, period_days):
    phases = np.array([centered_phase(row[0]["phase"]) for row in rows])
    time_hours = phases * period_days * 24.0
    states = np.array([row[0]["state"] for row in rows])
    wl = rows[0][1]
    depth_total = np.array([row[2] for row in rows])
    depth_atm = np.array([row[3] for row in rows])
    depth_east = np.array([row[4] for row in rows])
    depth_west = np.array([row[5] for row in rows])
    depth_opaque = np.array([row[6] for row in rows])
    flux = 1.0 - depth_total

    fig, ax = plt.subplots()
    for i, phase in enumerate(phases):
        state = STATE_LABELS.get(states[i], f"state {states[i]}")
        alpha = 0.35 if states[i] == 0 else 0.85
        ax.plot(wl, depth_total[i], label=f"{phase:+.5f} ({state})", alpha=alpha)

    ax.set_xlabel(r"$\lambda$ [$\mu$m]", fontsize=14)
    ax.set_ylabel("Transit depth", fontsize=14)
    ax.set_xscale("log")
    ax.legend(fontsize=8, ncol=2)
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

    fig, ax = plt.subplots()
    ax.plot(time_hours, np.mean(depth_total, axis=1), "ko-", label="Total")
    ax.plot(time_hours, np.mean(depth_opaque, axis=1), "o-", color="darkorange", label="Opaque body")
    ax.plot(time_hours, np.mean(depth_atm, axis=1), "o-", color="0.4", label="Atmosphere")
    ax.plot(time_hours, np.mean(depth_east, axis=1), "o-", color="firebrick", label="East limb (atm)")
    ax.plot(time_hours, np.mean(depth_west, axis=1), "o-", color="royalblue", label="West limb (atm)")
    ax.set_xlabel("Time from mid-transit [hours]", fontsize=14)
    ax.set_ylabel("Band mean transit depth", fontsize=14)
    ax.legend(fontsize=9)
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

    fig, ax = plt.subplots()
    flux_band = np.mean(flux, axis=1)
    ax.plot(time_hours, flux_band, "ko-", label="Band mean")
    sample_idx = np.linspace(0, len(wl) - 1, min(4, len(wl)), dtype=int)
    for idx in sample_idx:
        ax.plot(time_hours, flux[:, idx], "-", alpha=0.45, label=f"{wl[idx]:.3g} um")
    ax.set_xlabel("Time from mid-transit [hours]", fontsize=14)
    ax.set_ylabel("Normalised flux", fontsize=14)
    ax.legend(fontsize=9)
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(wl, time_hours, flux, shading="auto", cmap="inferno_r")
    ax.set_xlabel(r"$\lambda$ [$\mu$m]", fontsize=14)
    ax.set_ylabel("Time from mid-transit [hours]", fontsize=14)
    ax.set_xscale("log")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Normalised flux")
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot 3D_sph_trans_lc Transit_*.txt output files."
    )
    parser.add_argument(
        "--pattern",
        default="Transit_*.txt",
        help="Glob pattern for transit light-curve output files.",
    )
    parser.add_argument(
        "--rs",
        type=float,
        default=DEFAULT_RS,
        help="Deprecated; transit depths are already normalized by the Fortran output.",
    )
    parser.add_argument(
        "--period-days",
        type=float,
        default=DEFAULT_PERIOD_DAYS,
        help="Orbital period in days. Default is 1.21987 for WASP-33b.",
    )
    args = parser.parse_args()

    rows = load_transit_files(args.pattern)
    plot_transit_results(rows, args.period_days)


if __name__ == "__main__":
    main()
