"""Step-by-step derivation of the total Earth spectrum (CXB reflection + CR-induced emission) seen 600 km above
Rome (41.9 N, 12.5 E) at solar minimum (phi = 0.40 GV). docs/worked_example.md explains each step.

Usage: python examples/rome_solar_minimum.py [--lat 41.9 --lon 12.5 --alt 600 --phi 0.40 --epoch 2025]
Writes the spectrum to rome_600km_solar_min.txt (and a plot if matplotlib is installed).
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import albedo_earth as ae  # noqa: E402

BANDS = [(1, 10), (10, 100), (100, 1000), (1000, 10000), (25, 300)]


def table(e, rows, ref=None):
    """Band-integrated values [ph cm-2 s-1] of each (label, F) row, or their ratio to ``ref``."""
    out = ["        " + f"{'':14s}" + "".join(f"{f'{lo:g}-{hi:g} keV':>15s}" for lo, hi in BANDS)]
    for label, f in rows:
        val = [ae.band(e, f, lo, hi) / (1.0 if ref is None else ae.band(e, ref, lo, hi)) for lo, hi in BANDS]
        out.append("        " + f"{label:14s}" + "".join(f"{x:15.4g}" if ref is None else f"{x:15.3f}" for x in val))
    return "\n".join(out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--lat", type=float, default=41.9)
    p.add_argument("--lon", type=float, default=12.5)
    p.add_argument("--alt", type=float, default=600.0)
    p.add_argument("--phi", type=float, default=0.40, help="modulation potential [GV], VP15 scale")
    p.add_argument("--epoch", type=float, default=2025.0, help="epoch of the IGRF cutoff grid")
    p.add_argument("--out", default="rome_600km_solar_min.txt")
    a = p.parse_args(argv)

    # Step 1: inputs --------------------------------------------------------------------------------------------
    em = ae.Emission("gruber99")                        # reference hadronic model (em.reference)
    cut = ae.Cutoff(a.epoch)
    e = em.e_edges
    d = ae.R_E + a.alt
    print(f"Step 1  observer {a.alt:g} km above ({a.lat:g} N, {a.lon:g} E): d = {d:g} km; phi = {a.phi:g} GV; "
          f"IGRF {a.epoch:g}; CXB {em.cxb_model}; hadronic model {em.model}")

    # Step 2: geometry ------------------------------------------------------------------------------------------
    th_ground = np.degrees(np.arcsin(ae.R_E / d))
    th_top = np.degrees(np.arcsin(ae.R_TOP / d))
    om = ae.solid_angle_per_b(em.b_edges, d)
    print(f"Step 2  the solid Earth fills nadir angles < {th_ground:.2f} deg ({om[:20].sum():.3f} sr); the "
          f"atmosphere up to 400 km, < {th_top:.2f} deg ({om.sum():.3f} sr)")

    # Step 3: CXB reflection (independent of position and solar cycle) -------------------------------------------
    f_refl = em.refl @ om
    f_refl_err = np.sqrt((em.refl_err ** 2) @ om ** 2)
    print("Step 3  CXB reflection, F(E) = sum_b I(E, b) Omega_b, band integrals [ph cm-2 s-1]")
    print(table(e, [("reflection", f_refl)]))

    # Step 4: local cutoff --------------------------------------------------------------------------------------
    rc0 = float(cut.at(a.lat, a.lon))
    v = ae.lines_of_sight(a.lat, a.lon, a.alt, em.b_edges)
    i_refl, i_cr = ae.sky(em, v, cut, a.phi)
    disk = v["b"] < ae.R_E
    print(f"Step 4  vertical cutoff below the observer R_c = {rc0:.2f} GV; over the visible Earth "
          f"{v['rc'][disk].min():.2f}-{v['rc'][disk].max():.2f} GV (disk), "
          f"{v['rc'].min():.2f}-{v['rc'].max():.2f} GV (with the limb)")

    # Step 5a: CR emission, quick estimate with the nadir cutoff everywhere ------------------------------------
    f_cr_nadir = em.cr_uniform(a.phi, rc0) @ om
    print(f"Step 5a CR emission with the nadir cutoff everywhere (R_c = {rc0:.2f} GV) [ph cm-2 s-1]")
    print(table(e, [("CR (5a)", f_cr_nadir)]))

    # Step 5b: CR emission, cutoff of each line of sight (reference) -------------------------------------------
    f_cr = ae.flux(i_cr, v)
    print("Step 5b CR emission with the cutoff of each line of sight [ph cm-2 s-1]")
    print(table(e, [("CR (5b)", f_cr)]))
    print(table(e, [("ratio 5a/5b", f_cr_nadir)], ref=f_cr))

    # Step 6: total ---------------------------------------------------------------------------------------------
    f_tot = f_refl + f_cr
    f_plane = ae.flux(i_refl + i_cr, v, cosine=True)
    print("Step 6  total = reflection + CR [ph cm-2 s-1]; last row: through a plane facing nadir")
    print(table(e, [("reflection", f_refl), ("CR", f_cr), ("total", f_tot), ("total, plane", f_plane)]))
    print(table(e, [("CR fraction", f_cr)], ref=f_tot))

    # Step 7: model uncertainties -------------------------------------------------------------------------------
    var = {}
    hadronic = [m for m in em.models if m != em.reference]          # alternative hadronic models of the file
    for name, kw in ([(m, dict(model=m)) for m in hadronic]
                     + [("lines", dict(line_scaled=True)), ("churazov07", dict(cxb_model="churazov07"))]):
        emv = ae.Emission(**kw)
        r, c = ae.sky(emv, v, cut, a.phi)
        var[name] = ae.flux(r + c, v)
    print("Step 7  total of the variants relative to the reference")
    print(table(e, list(var.items()), ref=f_tot))
    had_lo = np.minimum.reduce([f_tot] + [var[m] for m in hadronic])
    had_hi = np.maximum.reduce([f_tot] + [var[m] for m in hadronic])
    mb = (e[:-1] >= 25) & (e[1:] <= 300)
    stat = np.sqrt(((f_refl_err * np.diff(e))[mb] ** 2).sum())    # independent bins, added in quadrature
    print(f"        statistical error of the reflection, 25-300 keV: {stat:.2g} ph cm-2 s-1 "
          f"({stat / ae.band(e, f_refl, 25, 300):.2%}); CR: < 1 % in broad bands (event replicas, not tabulated)")

    # Output ----------------------------------------------------------------------------------------------------
    m = (e[:-1] >= 1 - 1e-9) & (e[1:] <= 1e4 * (1 + 1e-9))
    with open(a.out, "w") as f:
        f.write(f"# Earth emission seen {a.alt:g} km above ({a.lat:g} N, {a.lon:g} E); phi = {a.phi:g} GV (VP15); "
                f"IGRF {a.epoch:g}; CXB {em.cxb_model}; CR {em.model}\n")
        f.write("# omnidirectional flux, sum of I dOmega over the visible Earth [ph cm-2 s-1 keV-1]\n")
        f.write("# E_lo E_hi [keV]  F_refl  F_CR  F_total  F_total_hadronic_min  F_total_hadronic_max  "
                "F_total_line_scaled\n")
        for k in np.nonzero(m)[0]:
            f.write(f"{e[k]:10.4g} {e[k + 1]:10.4g} {f_refl[k]:11.4e} {f_cr[k]:11.4e} {f_tot[k]:11.4e} "
                    f"{had_lo[k]:11.4e} {had_hi[k]:11.4e} {var['lines'][k]:11.4e}\n")
    print(f"Written {a.out}")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    ec = em.e_centres[m]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.fill_between(ec, (ec ** 2 * had_lo[m]), (ec ** 2 * had_hi[m]), color="0.8", lw=0, label="Hadronic models")
    ax.step(ec, ec ** 2 * f_tot[m], where="mid", color="k", label="Total")
    ax.step(ec, ec ** 2 * f_refl[m], where="mid", color="tab:blue", label="CXB reflection")
    ax.step(ec, ec ** 2 * f_cr[m], where="mid", color="tab:red", label="Cosmic-ray induced")
    ax.set(xscale="log", yscale="log", xlabel="Energy [keV]", ylabel=r"$E^2 F(E)$ [keV cm$^{-2}$ s$^{-1}$]",
           title=f"{a.alt:g} km above ({a.lat:g} N, {a.lon:g} E), φ = {a.phi:g} GV", ylim=(1e-2, None))
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.splitext(a.out)[0] + ".png", dpi=130)


if __name__ == "__main__":
    main()
