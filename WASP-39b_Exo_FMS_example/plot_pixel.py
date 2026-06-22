import re
from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
from matplotlib import cm
from scipy.io import FortranFile


plt.rc("font", family="sans-serif")
plt.rc("font", serif="Helvetica Neue")
plt.rc("text", usetex="false")


def read_nml_int(name, default):
    path = Path("CMCRT.nml")
    if not path.exists():
        return default
    match = re.search(rf"\b{name}\s*=\s*([0-9]+)", path.read_text())
    return int(match.group(1)) if match else default


def read_emission_file(path):
    head = np.loadtxt(path, max_rows=1)
    data = np.atleast_2d(np.loadtxt(path, skiprows=1))
    ltot = data[:, 3] if data.shape[1] >= 4 else data[:, 2]
    return {
        "wl": data[:, 0],
        "frac": data[:, 1],
        "ltot": ltot,
        "r_top": float(head[2]),
        "viewphi": float(head[3]) if len(head) > 3 else np.nan,
    }


def brightness_temperature(wl_um, flux):
    c0 = 2.99792458e10
    h = 6.626176e-27
    kb = 1.3080662e-16
    wl_cm = wl_um * 1e-4
    return (h * c0) / (kb * wl_cm) / np.log(
        1.0 + ((2.0 * h * c0**2) / (flux / np.pi * wl_cm**5))
    )


def bandpass_weights(wl):
    band_files = sorted(Path(".").glob("*response*.txt"))
    if not band_files:
        return np.ones_like(wl) / len(wl)
    data = np.loadtxt(band_files[0])
    response = np.interp(wl, data[:, 0], data[:, 1], left=0.0, right=0.0)
    norm = np.trapz(response, wl)
    if norm <= 0.0:
        return np.ones_like(wl) / len(wl)
    return response / norm


em_files = sorted(Path(".").glob("Em_*.txt"))
if not em_files:
    raise FileNotFoundError("No Em_*.txt files found")

first = read_emission_file(em_files[0])
wl = first["wl"]
n_wl = len(wl)
xpix = read_nml_int("xpix", 200)
ypix = read_nml_int("ypix", 200)
weights = bandpass_weights(wl)

area = ((first["r_top"] * 1.01) / xpix) * ((first["r_top"] * 1.01) / ypix)

for n, fname in enumerate(em_files):
    em = read_emission_file(fname)
    image_name = f"f_im_{n+1:03d}.txt"
    print(image_name)
    f = FortranFile(image_name, "r")
    bt_band = np.zeros((xpix, ypix))

    for l in range(n_wl):
        f_pix = f.read_record(np.float32).reshape((xpix, ypix), order="F")
        flux_pix = f_pix * em["ltot"][l] / area
        bt_band += brightness_temperature(wl[l], flux_pix) * weights[l]

    fig, ax = plt.subplots()
    cmap = cm.get_cmap("RdYlBu_r").copy()
    cmap.set_under(color="black")
    pmap = ax.imshow(bt_band, vmin=np.nanpercentile(bt_band, 1), vmax=np.nanpercentile(bt_band, 99), cmap=cmap)
    cbar = fig.colorbar(pmap, extend="both")
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(r"T$_{\rm b}$ [K]", fontsize=14)
    ax.set_ylabel("y pixel", fontsize=14)
    ax.set_xlabel("x pixel", fontsize=14)
    ax.set_title(f"Viewphi: {em['viewphi']:.1f}", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    fig.tight_layout()

plt.show()
