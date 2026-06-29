from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
from scipy.stats import binned_statistic


plt.rc("font", family="sans-serif")
plt.rc("font", serif="Helvetica Neue")
plt.rc("text", usetex="false")

RSUN = 6.957e10
RS = 0.932 * RSUN


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


def brightness_temperature(wl_um, flux):
    c0 = 2.99792458e10
    h = 6.626176e-27
    kb = 1.3080662e-16
    wl_cm = wl_um * 1e-4
    return (h * c0) / (kb * wl_cm) / np.log(
        1.0 + ((2.0 * h * c0**2) / (flux / np.pi * wl_cm**5))
    )


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


def wavelength_edges_for(wl):
    edges = np.loadtxt("wavelengths_R100_UV_edge.txt", skiprows=1)[:, 1]
    centers = np.loadtxt("wavelengths.wl", skiprows=1)[:, 1]
    start = int(np.argmin(np.abs(centers - wl[0])))
    return edges[start:start + len(wl) + 1]


em_files = sorted(Path(".").glob("Em_*.txt"))
if not em_files:
    raise FileNotFoundError("No Em_*.txt files found")

first = read_emission_file(em_files[0])
wl = first["wl"]
n_ph = len(em_files)
n_wl = len(wl)

fp = np.zeros((n_ph, n_wl))
fp_occ = np.zeros((n_ph, n_wl))
bt = np.zeros((n_ph, n_wl))
viewphi = np.zeros(n_ph)
phase = np.zeros(n_ph)
occ_state = np.zeros(n_ph, dtype=int)

for n, fname in enumerate(em_files):
    print(fname)
    em = read_emission_file(fname)
    viewphi[n] = em["viewphi"]
    phase[n] = em["phase"]
    occ_state[n] = em["occ_state"]
    fp[n, :] = (em["frac"] * em["ltot"]) / em["r_top"]**2
    fp_occ[n, :] = (em["frac_occ"] * em["ltot"]) / em["r_top"]**2
    bt[n, :] = brightness_temperature(wl, fp[n, :])

fig, ax = plt.subplots()
for n in range(n_ph):
    ax.plot(wl, fp[n, :], label=f"Phase {phase[n]:.4f} (occ {occ_state[n]})")
ax.set_xlabel(r"$\lambda$ [$\mu$m]", fontsize=14)
ax.set_ylabel(r"F$_{\rm p}$ [erg s$^{-1}$ cm$^{-2}$ cm$^{-1}$]", fontsize=14)
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend(fontsize=8)
ax.tick_params(axis="both", which="major", labelsize=12)
fig.tight_layout()

fig, ax = plt.subplots()
for n in range(n_ph):
    ax.plot(wl, bt[n, :], label=f"Phase {phase[n]:.4f} (occ {occ_state[n]})")
ax.set_xlabel(r"$\lambda$ [$\mu$m]", fontsize=14)
ax.set_ylabel(r"T$_{\rm b}$ [K]", fontsize=14)
ax.set_xscale("log")
ax.legend(fontsize=8)
ax.tick_params(axis="both", which="major", labelsize=12)
fig.tight_layout()

fs_binned = load_stellar_spectrum(wl)
fig, ax = plt.subplots()
for n in range(n_ph):
    fpfs = (first["rp"] / RS) ** 2 * (fp[n, :] / fs_binned)
    fpfs_occ = (first["rp"] / RS) ** 2 * (fp_occ[n, :] / fs_binned)
    ax.plot(wl, fpfs, label=f"Phase {phase[n]:.4f} (occ {occ_state[n]})")
    ax.plot(wl, fpfs_occ, ls="--", alpha=0.7)
ax.set_xlabel(r"$\lambda$ [$\mu$m]", fontsize=14)
ax.set_ylabel(r"F$_{\rm p}$/F$_{\star}$", fontsize=14)
ax.set_xscale("log")
ax.legend(fontsize=8)
ax.tick_params(axis="both", which="major", labelsize=12)
fig.tight_layout()

plt.show()
