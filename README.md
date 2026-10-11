# earth-albedo-xgamma: the X- and γ-ray emission of the Earth, 1 keV–10 MeV

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23158832.svg)](https://doi.org/10.5281/zenodo.23158832)

Tables of the photon emission of the Earth's atmosphere seen from space, from a Geant4 Monte Carlo simulation
(albedo_sim): the **reflection of the cosmic X-ray background (CXB)** and the **emission induced by galactic cosmic
rays (CR)**, from 1 keV to 10 MeV (photons are recorded up to 10 GeV). The tables give the emergent intensity
against the photon energy and the impact parameter of the ray. With them you can compute the spectrum, and its
angular distribution, seen by any observer above 400 km, at any geographic position and phase of the solar cycle.

The release contains the tables (HDF5), the response functions of the standard production (to fold other CXB or
cosmic-ray spectra), their full description, a small Python module that maps the tables onto the sky of an
observer, and three examples. The method, the validation and the approximations are described in the
accompanying paper (reference to be added).

## Quick start

Requirements: Python ≥ 3.9 with numpy and h5py (matplotlib optional), see `requirements.txt`.

```sh
python examples/minimal_example.py        # whole-Earth spectrum from 500 km, uniform cutoff (numpy + h5py only)
python examples/rome_solar_minimum.py     # total spectrum 600 km above Rome at solar minimum, step by step
python examples/fold_responses.py         # fold the response functions with other CXB or cosmic-ray spectra
```

```python
import albedo_earth as ae

em = ae.Emission()                                   # Gruber CXB, reference CR spectra, FTFP_BERT_HP
cut = ae.Cutoff(2025.0)                              # IGRF-14 vertical cutoffs
v = ae.lines_of_sight(41.9, 12.5, 600.0, em.b_edges) # observer: lat, lon [deg], altitude [km]
i_refl, i_cr = ae.sky(em, v, cut, phi=0.40)         # intensities per (energy, line of sight)
f_total = ae.flux(i_refl + i_cr, v)                  # ph cm-2 s-1 keV-1, on the bins of em.e_edges
```

## Contents

| Path | Content |
|---|---|
| `data/albedo_cxb_reflection.h5` | CXB reflection: I(E, b), albedo A(E) and whole-Earth spectra, for five CXB normalisations |
| `data/albedo_cr_emission.h5` | CR-induced emission (reference spectra): I(E, b) on a grid of solar modulation φ and vertical cutoff R_c, per hadronic model and component |
| `data/albedo_cr_emission_s07.h5` | the same with the CR spectra of Sazonov et al. (2007), for comparison with that work |
| `data/cutoff_igrf_<epoch>.h5` | IGRF-14 vertical cutoff grids (2° × 2°, 35 km altitude); use `cutoff_igrf_2025.h5` for current missions |
| `data/albedo_reference_orbits.h5`, `data/reference_orbits/*.txt` | ready-made orbit averages for mission design: LEO-EQ (600 km, i = 5°) and LEO-SSO (550 km, i = 97.6°), solar minimum and maximum, with a systematic band |
| `data/responses/response_<species>.h5` | response functions of the standard production: CXB photons (`cxb`), protons, 4He (`alpha`), e⁻, e⁺ |
| `docs/format.md` | **every group, dataset and attribute of the HDF5 files, with shapes and units** |
| `docs/worked_example.md` | step-by-step derivation of the total spectrum at a given position and solar-cycle phase |
| `albedo_earth.py` | helper module (numpy + h5py): tables, cutoffs, lines of sight, fluxes |
| `examples/` | `minimal_example.py`, `rome_solar_minimum.py`, `fold_responses.py` |
| `LICENSE`, `CITATION.cff` | licence (CC BY 4.0) and how to cite |
| `MANIFEST.json` | list of the files with sizes, SHA-256 checksums, code version and provenance of the products |

In the HDF5 files every dataset carries `units` and `description` attributes, e.g.
`h5py.File("data/albedo_cr_emission.h5")["intensity_total"].attrs["units"]`.

## The central quantity: I(E, b)

All tables give the emergent **photon intensity** I(E, b), in **ph cm⁻² s⁻¹ sr⁻¹ keV⁻¹**, of a spherically symmetric
Earth, against the photon energy E (keV) and the **impact parameter b** of the ray, i.e., the distance of the ray
line from the centre of the Earth (km). Above the atmosphere photons travel along straight lines and the intensity
is conserved along each line; for a spherically symmetric Earth it depends only on b. The same table therefore
serves **any observer above 400 km** (the top of the model atmosphere):

- a line of sight at nadir angle θ from an observer at distance d from the centre of the Earth has **b = d sin θ**;
- the solid angle of the rays with impact parameter below b is **Ω(< b) = 2π [1 − (1 − b²/d²)^½]**;
- the flux from the whole Earth is **F(E) = Σ_b I(E, b) Ω_b(d)** (ph cm⁻² s⁻¹ keV⁻¹), with Ω_b the solid angle
  of the bin. This is the omnidirectional flux, i.e., the flux through a sphere of unit cross-section; for a flat
  detector facing nadir, weight each ray with cos θ.

Axes, common to all files:
- `e_edges` (141 values): energy bin edges, 1 keV–10 GeV, 20 bins per decade; spectra are bin averages per keV;
- `b_edges` (176 values): 20 disk bins (b < R_E = 6371 km, equal steps in μ = (1 − b²/R_E²)^½), then the limb in
  tangent altitude t = b − R_E: 0–150 km in 1 km steps, 150–400 km in 50 km steps. The limb is finely binned
  because the emission of the Earth peaks there above a few hundred keV.

### CXB reflection

The reflection does not depend on the position of the observer. The solar cycle changes the upper atmosphere, but
its effect on the reflection is below 0.5 % overall and 1–3 % on the 1–5 keV limb below 120 km tangent altitude.
The reflection is proportional to the CXB normalisation: `gruber99` (Gruber et al. 1999) is the reference;
`churazov07` (+10 %) and `ajello08` (+6 %) are the higher normalisations measured by Earth occultation with
INTEGRAL and Swift/BAT. Unscattered CXB photons crossing the limb are not included.

### Cosmic-ray induced emission

The CR emission depends on the local **geomagnetic cutoff rigidity** and on the **solar modulation**. The tables
give the emission of a whole Earth with a single vertical cutoff R_c (`rc_GV`, 0–17 GV) and modulation potential φ
(`phi_GV`, 0.3–1.5 GV). Since the showers develop within a few tens of kilometres of the impact point, the
emission at a place with cutoff R_c is that of the uniform Earth with the same R_c. For a real field, each line of
sight takes I(E, b) at the cutoff where it reaches 35 km (the altitude of the emission), or at its tangent point
if it passes above; the cutoffs come from the IGRF grids (R_c = `r_eq_GV` / `L_bmin`²). Above 17 GV (IGRF
reaches 17.7 GV over South-East Asia) the tables are extrapolated linearly from the 16–17 GV interval (v1.0.2;
earlier versions clamped at 17 GV). `albedo_earth.sky` does this; `docs/worked_example.md` explains it step by step.

- **Solar modulation.** φ is on the scale of the Vos & Potgieter (2015) local interstellar spectrum, i.e., of the
  monthly series of Usoskin et al. (2017) and Väisänen et al. (2023). Solar minimum is about 0.40 GV (annual
  means of 2009 and 2020), solar maximum about 0.80 GV (mean of the maxima of cycles 22–24; 0.66–1.02 GV).
  Interpolate linearly between grid values. The file `albedo_cr_emission_s07.h5` uses the older scale
  (φ_VP15 = 0.645 φ_old + 0.236 GV).
- **Components** (`intensity_components`): protons, He, Z > 2 nuclei, galactic e⁻ and e⁺, and re-entrant e±.
- **Hadronic models** (`intensity_total`): FTFP_BERT_HP (reference), QGSP_BIC_HP and QGSP_INCLXX_HP (the latter two
  below 10 GeV/nucleon, FTFP_BERT_HP above). Their spread is ≤ 4 % at 25–300 keV at low latitude and up to 20 % at
  the poles. The hadronic systematic uncertainty of the reference orbits is larger: it also includes these models
  above 10 GeV/nucleon and the PreCompound de-excitation in the Bertini cascade of FTFP_BERT_HP (folded estimates,
  paper Appendix C): −13 % at low latitude and −20 to −24 % at the poles at 25–300 keV, −26 to −39 % at 1–3 MeV.
  The default Bertini de-excitation emits no discrete nuclear lines, hence the weak 4.44 MeV line (below).
- **Nuclear lines.** Relative to the 511 keV line, the Geant4 yields of the 1.6, 2.3 and 6.1 MeV lines agree with
  the data within 20 %, while the 4.44 MeV line is about ten times too weak. For line work use
  `intensity_total_line_scaled`, which rescales the lines to the SMM and HEAO 3 measurements (511 keV ÷ 1.95,
  1.6, 2.3 and 6.1 MeV ÷ 1.6–2.0, 4.44 MeV × 6.2). v1.0.2 corrects the line fractions of v1.0–v1.0.1, which
  subtracted only half of the continuum under the lines.
  **Caveat:** this assumes the absolute SMM/HEAO 3 scale. The absolute MeV measurements of the atmosphere, all made
  between 1971 and 1989 (SMM, HEAO 3, 1972-076B, balloons), differ from each other by a factor of about five, and the
  nominal model lies within their range (0.5–2.8 times the data at 1–2 MeV); the line ratios and the MeV continuum
  shape (within 40 % of SMM relative to 511 keV) are the firmer results.

### Response functions

The emission tables are folded with fixed input spectra. To use other ones (a different CXB, other cosmic-ray
spectra or modulation models, a solar-energetic-particle event), fold the response functions in
`data/responses/`. For each primary species, `yield[c, i, j, k]` is the number of photons escaping in energy bin
j and impact-parameter bin k per primary of energy bin i entering the top of the atmosphere isotropically. For an
isotropic primary intensity I(E0) (cm⁻² s⁻¹ sr⁻¹ keV⁻¹; E0 is the total kinetic energy for nuclei):

- rate[j, k] = Σ_i π 4π R_src² (∫_i I dE0) Σ_c yield[c, i, j, k] photons s⁻¹, with R_src = R_E + 400 km;
- I(E_j, b_k) = rate[j, k] / (ΔE_j 8π² b Δb_k).

`albedo_earth.Response` does this (`Response("cxb").intensity(my_cxb)`). Nuclei with Z > 2 are folded on the He
response (a nucleus of mass number A at energy per nucleon e_n acts as A/4 He nuclei at e_n), and the cutoff is
applied to the input spectrum. The re-entrant e± are the e± escaping from the nuclear showers (`lepton_yield`),
re-injected isotropically and folded on the e⁻/e⁺ responses. `examples/fold_responses.py` reproduces the released
tables (within 0.02 % in broad bands) and shows a different CXB. Only the standard production is released
(FTFP_BERT_HP, o-Ps fraction 0.035); the other hadronic models are available as emission tables.

### Uncertainties

- **Statistical:** `I_E_b_err` for the reflection. For the CR emission, the photons of one shower are correlated
  across bins, so per-bin errors are not given; band errors from event replicas are in `band_brightness` (typically
  0.2–0.9 % in broad bands) and, for any folded quantity, from `yield_replicas` in the response files.
- **Systematic:** CXB normalisation (about ±10 %), hadronic model (see above), CR spectra (modified ISO 15390 vs the
  reference: +2 to +14 %), geomagnetic field inside the atmosphere (−21 % on the 1–10 MeV limb, ≤ 4 % elsewhere),
  cutoff model, line yields. The reference orbits carry these terms per energy bin (`shifts/`, `sys_up`, `sys_down`)
  and can be used as a guide for other orbits.

## Main approximations

Spherical, static US Standard Atmosphere 1976; sharp vertical cutoffs, the same in all directions (no penumbra, no
east–west effect; the directional cutoffs change the emission by 3 % or less); internal IGRF field only; no
geomagnetic field inside the atmosphere; Z > 2 nuclei by superposition of the He response; re-entrant e± re-injected
isotropically; Geant4 line yields (see the line-scaled variant). **Not included:** trapped particles (SAA), auroral
and precipitating electrons, solar X-rays reflected on the dayside, solar energetic particles, terrestrial gamma-ray
flashes, and the background produced in the instrument itself.

## Provenance

The tables were produced with albedo_sim (Geant4 11.4.1); the code version and the number of simulated primaries
of each product are in `MANIFEST.json` and in the file attributes (`chunks`, `n_primaries`, `settings`).

## Version history

- **v1.0.2** (2026-10-11):
  - *Line-scaled variant corrected.* The line fractions of v1.0–v1.0.1 subtracted only half of the continuum
    under each line; the corrected `line_scale_factor` acts only on the bin containing each line centre, and the
    4.44 MeV line is now scaled up (Geant4 under-produces it; see Nuclear lines). `intensity_total` and the other
    emission tables are unchanged.
  - `band_brightness` now integrates exactly 25–300 keV (partial bins); values are 2.2 % higher than in v1.0.1.
  - Reference orbits: rays split at the b edges of the tables (exact sums; v1.0.1 was 5–9 % low per bin at
    1–4 keV, about 1 % in the bands); cutoffs above 17 GV extrapolated linearly instead of clamped; hadronic
    systematic term widened (see Hadronic models), so `sys_down` is larger above 100 keV.
  - `albedo_earth.lines_of_sight`: the azimuth `psi` now runs from east towards north as documented (v1.0.1 was
    mirrored east–west; integrals over the full azimuth were not affected); cutoffs above 17 GV extrapolated.
  - Documentation: per-bin accuracy of the response functions qualified (0.5 % above 25 keV).
  - IGRF cutoff grid added for epoch 1971.5 (balloon comparison of the paper).
- **v1.0.1** (2026-10-05): Zenodo DOIs in README and `CITATION.cff`; data identical to v1.0.
- **v1.0** (2026-10-05): first release.

## Citation and licence

The content of this release (data, documentation and code) is licensed under the Creative Commons Attribution 4.0
International licence (CC BY 4.0, see `LICENSE`): you may share and adapt it for any purpose, provided that you give
appropriate credit. If you use these data, please cite the accompanying paper and this release (`CITATION.cff`).

The release is archived on Zenodo: [doi:10.5281/zenodo.23158832](https://doi.org/10.5281/zenodo.23158832)
(all versions; version 1.0 is [doi:10.5281/zenodo.23158833](https://doi.org/10.5281/zenodo.23158833), version
1.0.1 [doi:10.5281/zenodo.23159121](https://doi.org/10.5281/zenodo.23159121)).
Paper reference: to be added.
