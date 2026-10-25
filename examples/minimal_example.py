"""Minimal example (numpy + h5py only): whole-Earth spectrum seen from 500 km altitude, for a uniform cutoff.

Total = CXB reflection + CR-induced emission, for a whole Earth with one vertical cutoff (R_c = 7 GV) at solar
modulation phi = 0.6 GV. For a real (position-dependent) cutoff see rome_solar_minimum.py.
"""

import os

import h5py
import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
R_E, d = 6371.0, 6371.0 + 500.0                        # km: observer 500 km above the surface

with h5py.File(os.path.join(DATA, "albedo_cr_emission.h5")) as h:
    e, b = h["e_edges"][:], h["b_edges"][:]
    phis, rcs = h["phi_GV"][:], h["rc_GV"][:]
    i_cr = h["intensity_total"][0, np.argmin(abs(phis - 0.6)), np.argmin(abs(rcs - 7.0))]   # reference model, [E, b]
with h5py.File(os.path.join(DATA, "albedo_cxb_reflection.h5")) as h:
    i_refl = h["gruber99/I_E_b"][:]                                                          # [E, b]

om = np.diff(2 * np.pi * (1 - np.sqrt(1 - np.clip(b / d, 0, 1) ** 2)))   # solid angle of each b bin [sr]
f_refl, f_cr = i_refl @ om, i_cr @ om                                     # ph cm-2 s-1 keV-1
m = (e[:-1] >= 25) & (e[1:] <= 300)
for name, f in (("reflection", f_refl), ("CR", f_cr), ("total", f_refl + f_cr)):
    print(f"25-300 keV {name:10s}: {(f[m] * np.diff(e)[m]).sum():.3f} ph cm-2 s-1")
