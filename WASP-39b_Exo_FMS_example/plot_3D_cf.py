import re
from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
from scipy.io import FortranFile


def read_nml_int(name, default):
    path = Path("CMCRT.nml")
    if not path.exists():
        return default
    match = re.search(rf"\b{name}\s*=\s*([0-9]+)", path.read_text())
    return int(match.group(1)) if match else default


def profile_name():
    path = Path("CMCRT.nml")
    if not path.exists():
        return "FMS.prf"
    match = re.search(r"\bexp_name\s*=\s*'([^']+)'", path.read_text())
    return f"{match.group(1)}.prf" if match else "FMS.prf"


em_file = Path("Em_001.txt")
trans_file = Path("Transmission.txt")
if em_file.exists():
    wl = np.loadtxt(em_file, skiprows=1)[:, 0]
elif trans_file.exists():
    wl = np.loadtxt(trans_file, skiprows=1)[:, 0]
else:
    raise FileNotFoundError("Need Em_001.txt or Transmission.txt to read wavelengths")

n_wl = len(wl)
nlat = read_nml_int("n_theta", 97) - 1
nlon = read_nml_int("n_phi", 193) - 1
nlay = read_nml_int("n_lay", 54)
ni = nlat * nlon * nlay

cf_path = Path("cf_001.txt")
if not cf_path.exists():
    raise FileNotFoundError("No cf_001.txt found")

f = FortranFile(cf_path, "r")
cf = np.zeros((n_wl, nlay, nlon, nlat))
for l in range(n_wl):
    cf_1d = f.read_reals(np.float32)
    if len(cf_1d) != ni:
        raise ValueError(f"cf_001.txt record {l} has {len(cf_1d)} values; expected {ni}")
    denom = np.sum(cf_1d)
    if denom != 0.0:
        cf[l, :, :, :] = np.reshape(cf_1d, [nlay, nlon, nlat], order="F") / denom

prf = np.loadtxt(profile_name(), skiprows=27)
p_1d = prf[:, 1]
p_grid = np.reshape(p_1d, [nlay, nlon, nlat], order="F")

ilat = int(nlat / 2)
ilon = 0
lwl = [0.5, 1.0, 5.0, 10.0, 20.0]
iwl = np.clip(np.searchsorted(wl, lwl), 0, n_wl - 1)

fig, ax = plt.subplots()
for idx in iwl:
    ax.plot(cf[idx, :, ilon, ilat], p_grid[:, ilon, ilat], label=f"{wl[idx]:.3g} um")
ax.set_xlabel("Fractional Contribution", fontsize=14)
ax.set_ylabel("P [bar]", fontsize=14)
ax.set_xscale("log")
ax.set_yscale("log")
ax.invert_yaxis()
ax.legend()
ax.tick_params(axis="both", which="major", labelsize=12)
fig.tight_layout()

fig, ax = plt.subplots()
cmap = sns.color_palette("mako", as_cmap=True)
cf_swap = np.swapaxes(cf[:, :, ilon, ilat], 0, 1)
positive = cf_swap[cf_swap > 0]
floor = -10 if positive.size == 0 else max(-10, np.log10(np.min(positive)))
ceil = -1 if positive.size == 0 else np.log10(np.max(positive))
levels = np.linspace(floor, ceil, 20)
con = ax.contourf(wl, p_grid[:, ilon, ilat], np.log10(np.maximum(cf_swap, 10**floor)), levels=levels, extend="both", cmap=cmap)
fig.colorbar(con, ax=ax)
ax.set_xlabel(r"$\lambda$ [$\mu$m]", fontsize=14)
ax.set_ylabel("P [bar]", fontsize=14)
ax.set_xscale("log")
ax.set_yscale("log")
ax.invert_yaxis()
ax.tick_params(axis="both", which="major", labelsize=12)
fig.tight_layout()

plt.show()
