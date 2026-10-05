"""Fold the response functions with input spectra of your choice (numpy + h5py only).

1. CXB reflection: the reference CXB (Gruber et al. 1999) reproduces the table of albedo_cxb_reflection.h5; a CXB
   with another shape (Ajello et al. 2008) gives a new I(E, b). Their fit describes the BAT band (15-200 keV);
   extrapolated to MeV energies it falls below Gruber, hence the lower reflection there.
2. Cosmic-ray emission: the reference spectra (force field on the local interstellar spectra of Vos & Potgieter
   2015 for p, Bisschoff et al. 2019 for He, e-, e+) at phi = 0.40 GV and R_c = 6 GV reproduce
   albedo_cr_emission.h5, with all six components: p, He, Z > 2 (on the He response), e-, e+ and re-entrant e+-.
Replace the spectra below with your own to obtain the emission for other inputs.

Usage: python examples/fold_responses.py
"""

import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import albedo_earth as ae  # noqa: E402

D = ae.DATA


def band_flux(i, e, b, d_km, lo, hi):
    """Whole-Earth flux [ph cm-2 s-1] in [lo, hi] keV seen from distance d: sum over b and E of I Omega_b dE."""
    return ae.band(e, i @ ae.solid_angle_per_b(b, d_km), lo, hi)


# --- 1. CXB ----------------------------------------------------------------------------------------------------
def gruber99(e):
    """CXB photon intensity of Gruber et al. (1999) [ph cm-2 s-1 sr-1 keV-1], E in keV."""
    e = np.asarray(e, float)
    lo = 7.877 * e ** -1.29 * np.exp(-e / 41.13)
    hi = 0.0259 * (e / 60) ** -5.5 + 0.504 * (e / 60) ** -1.58 + 0.0288 * (e / 60) ** -1.05
    return np.where(e <= 60, lo, hi / e)


def ajello08(e):
    """CXB of Ajello et al. (2008), Swift/BAT: C / ((E/E_b)^G1 + (E/E_b)^G2)."""
    e = np.asarray(e, float)
    return 0.102 / ((e / 29.99) ** 1.32 + (e / 29.99) ** 2.88)


# --- 2. cosmic rays (reference spectra; E0 = total kinetic energy) -------------------------------------------------
M_P, M_E, M_N_HE = 0.938272, 0.000511, 0.932           # GeV
HEAVY_FACTOR = 0.068 * 14 / 4 + 0.026 * 24 / 4 + 0.0054 * 56 / 4   # Z > 2 per He at equal E/n, x A/4


def _beta2(t, m):
    return 1 - (m / (t + m)) ** 2


def lis_p(t):                    # Vos & Potgieter (2015), (m2 s sr GeV)-1, t [GeV]
    return 2700 * t ** 1.12 / _beta2(t, M_P) * ((t + 0.67) / 1.67) ** -3.93


def lis_he(t):                   # Bisschoff et al. (2019) eq. 15, per GeV/n, t [GeV/n]
    return 163.4 / _beta2(t, M_N_HE) * t ** 1.1 * ((t ** 0.97 + 0.58 ** 0.97) / (1 + 0.58 ** 0.97)) ** -4.0


def lis_e(t):                    # Bisschoff et al. (2019) eq. 12
    return 255.0 / _beta2(t, M_E) * t ** -1 * ((t + 0.63) / 1.63) ** -2.43 + 6.4 * t ** 2 * ((t + 15.0) / 16.0) ** -26.0


def lis_pos(t):                  # Bisschoff et al. (2019) eq. 13
    return (25.0 / _beta2(t, M_E) * t ** 0.1 * ((t ** 1.1 + 0.2 ** 1.1) / (1 + 0.2 ** 1.1)) ** -3.31
            + 23.0 * t ** 0.5 * ((t + 2.2) / 3.2) ** -9.5)


def force_field(lis, t, big_phi, m):
    return lis(t + big_phi) * t * (t + 2 * m) / ((t + big_phi) * (t + big_phi + 2 * m))


def cr_intensity(species, phi, rc):
    """I(E0) [cm-2 s-1 sr-1 keV-1] at 1 AU above the vertical cutoff rc [GV], E0 in keV (total kinetic energy)."""
    def f(e0_kev):
        t = np.asarray(e0_kev, float) * 1e-6                          # GeV
        if species == "proton":
            j, r = force_field(lis_p, t, phi, M_P), np.sqrt(t * (t + 2 * M_P))
        elif species == "alpha":                                      # Phi = (Z/A) phi; per GeV/n -> per GeV
            e_n = t / 4.0
            j, r = force_field(lis_he, e_n, phi / 2.0, M_N_HE) / 4.0, 2.0 * np.sqrt(e_n * (e_n + 2 * M_P))
        else:
            j = force_field(lis_e if species == "e-" else lis_pos, t, phi, M_E)
            r = np.sqrt(t * (t + 2 * M_E))
        return np.where(r >= rc, j * 1e-4 * 1e-6, 0.0)               # m-2 -> cm-2, GeV-1 -> keV-1
    return f


def cr_emission(phi, rc, resp):
    """I(E, b) per component [ph cm-2 s-1 sr-1 keV-1] for a whole Earth with cutoff rc."""
    out = {}
    nuclei = {"p": ("proton", 1.0), "he": ("alpha", 1.0), "heavy": ("alpha", HEAVY_FACTOR)}
    lep = 0.0
    for comp, (sp, k) in nuclei.items():
        f = cr_intensity(sp, phi, rc)
        out[comp] = k * resp[sp].intensity(f)
        lep = lep + k * resp[sp].lepton_rate(f)                         # escaping e+- [s-1 per lepton bin]
    for sp in ("e-", "e+"):
        out[sp] = resp[sp].intensity(cr_intensity(sp, phi, rc))
    i_em, i_ep = ae.reentrant_intensities(lep, resp["proton"].lep_edges, resp["e-"].r_src_km)
    out["re"] = resp["e-"].intensity(i_em) + resp["e+"].intensity(i_ep)  # re-entrant: no cutoff
    return out


def main():
    d = ae.R_E + 600.0
    bands = [(10, 100), (100, 1000), (1000, 10000), (25, 300)]
    print("Fluxes seen from 600 km [ph cm-2 s-1]; " + ", ".join(f"{lo:g}-{hi:g}" for lo, hi in bands) + " keV")

    rc = ae.Response("cxb")
    e, b = rc.e_edges, rc.b_edges
    i_g = rc.intensity(gruber99)
    with h5py.File(os.path.join(D, "albedo_cxb_reflection.h5")) as h:
        tab = h["gruber99/I_E_b"][:]
    print("CXB, Gruber 1999, folded / table:      " + "  ".join(
        f"{band_flux(i_g, e, b, d, lo, hi) / band_flux(tab, e, b, d, lo, hi):8.5f}" for lo, hi in bands))
    i_a = rc.intensity(ajello08)
    print("CXB, Ajello 2008 / Gruber 1999:        " + "  ".join(
        f"{band_flux(i_a, e, b, d, lo, hi) / band_flux(i_g, e, b, d, lo, hi):8.5f}" for lo, hi in bands))

    resp = {sp: ae.Response(sp) for sp in ("proton", "alpha", "e-", "e+")}
    phi, rcut = 0.40, 6.0
    comp = cr_emission(phi, rcut, resp)
    total = sum(comp.values())
    with h5py.File(os.path.join(D, "albedo_cr_emission.h5")) as h:
        j, k = list(h["phi_GV"][:]).index(phi), list(h["rc_GV"][:]).index(rcut)
        tab = h["intensity_total"][0, j, k].astype(float)
    print(f"CR, phi = {phi} GV, R_c = {rcut} GV, folded / table: " + "  ".join(
        f"{band_flux(total, e, b, d, lo, hi) / band_flux(tab, e, b, d, lo, hi):8.5f}" for lo, hi in bands))
    print("CR components, 25-300 keV: " + ", ".join(
        f"{c} {band_flux(i, e, b, d, 25, 300):.4f}" for c, i in comp.items()))


if __name__ == "__main__":
    main()
