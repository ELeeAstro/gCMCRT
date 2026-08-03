from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
from scipy.stats import binned_statistic


RSUN = 6.957e10
RS = 0.895 * RSUN


plt.rc("font", family="sans-serif")
plt.rc("font", serif="Helvetica Neue")
plt.rc("text", usetex="false")


def read_emission_file(path):
    head = np.atleast_1d(np.loadtxt(path, max_rows=1))
    data = np.atleast_2d(np.loadtxt(path, skiprows=1))

    frac_occ = data[:, 2] if data.shape[1] >= 4 else data[:, 1]
    ltot = data[:, 3] if data.shape[1] >= 4 else data[:, 2]

    phase = float(head[5]) if len(head) > 5 else float(head[3])

    return {
        "wl": data[:, 0],
        "frac": data[:, 1],
        "frac_occ": frac_occ,
        "ltot": ltot,
        "rp": float(head[1]),
        "r_top": float(head[2]),
        "viewphi": float(head[3]) if len(head) > 3 else np.nan,
        "viewthet": float(head[4]) if len(head) > 4 else np.nan,
        "phase": phase,
        "occ_state": int(head[6]) if len(head) > 6 else 0,
        "xstar": float(head[7]) if len(head) > 7 else np.nan,
        "ystar": float(head[8]) if len(head) > 8 else np.nan,
        "zstar": float(head[9]) if len(head) > 9 else np.nan,
        "separation": float(head[10]) if len(head) > 10 else np.nan,
    }


def centered_phase(phase):
    if phase < 0.0:
        return phase
    return ((phase + 0.5) % 1.0) - 0.5


def wavelength_edges_for(wl):
    edges = np.loadtxt("wavelengths_R100_UV_edge.txt", skiprows=1)[:, 1]
    centers = np.loadtxt("wavelengths.wl", skiprows=1)[:, 1]
    start = int(np.argmin(np.abs(centers - wl[0])))
    return edges[start:start + len(wl) + 1]


def load_stellar_spectrum(wl):
    candidates = sorted(Path(".").glob("*stellar*spectrum*.txt"))
    if not candidates:
        raise FileNotFoundError(
            "No stellar spectrum found. Add a file matching '*stellar*spectrum*.txt' "
            "with wavelength [um] and stellar flux [erg/s/cm2/um] columns."
        )
    data = np.loadtxt(candidates[0], skiprows=1)
    wl_s = data[:, 0]
    fs = data[:, 1] * 1e4
    edges = wavelength_edges_for(wl)
    fs_binned, _, _ = binned_statistic(wl_s, fs, bins=edges)
    return fs_binned


em_files = sorted(Path(".").glob("Em_*.txt"))
if not em_files:
    raise FileNotFoundError("No Em_*.txt files found")

first = read_emission_file(em_files[0])
n_ph = len(em_files)
n_wl = len(first["wl"])
fs_binned = load_stellar_spectrum(first["wl"])

phase = np.zeros(n_ph)
occ_state = np.zeros(n_ph, dtype=int)
mean_fpfs = np.zeros(n_ph)
mean_fpfs_occ = np.zeros(n_ph)

for n, fname in enumerate(em_files):
    print(fname)
    em = read_emission_file(fname)
    fp = (em["frac"] * em["ltot"]) / em["r_top"]**2
    fp_occ = (em["frac_occ"] * em["ltot"]) / em["r_top"]**2
    fpfs = (em["rp"] / RS) ** 2 * (fp / fs_binned)
    fpfs_occ = (em["rp"] / RS) ** 2 * (fp_occ / fs_binned)
    phase[n] = centered_phase(em["phase"])
    occ_state[n] = em["occ_state"]
    mean_fpfs[n] = np.mean(fpfs)
    mean_fpfs_occ[n] = np.mean(fpfs_occ)

order = np.argsort(phase)
phase = phase[order]
mean_fpfs = mean_fpfs[order]
mean_fpfs_occ = mean_fpfs_occ[order]
occ_state = occ_state[order]

fig, ax = plt.subplots()
ax.plot(phase, mean_fpfs, "o-", label="Unocculted")
ax.scatter(phase, mean_fpfs_occ, c=occ_state, label="Occulted")
ax.set_xlabel("Orbital phase", fontsize=14)
ax.set_ylabel(r"Band mean F$_{\rm p}$/F$_{\star}$", fontsize=14)
ax.legend()
ax.tick_params(axis="both", which="major", labelsize=12)
fig.tight_layout()

plt.show()
