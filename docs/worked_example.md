# Worked example: the Earth seen 600 km above Rome at solar minimum

This page derives, step by step, the total spectrum of the Earth (CXB reflection + cosmic-ray induced emission)
seen by an observer 600 km above Rome (41.9° N, 12.5° E) at solar minimum. The script
`examples/rome_solar_minimum.py` performs the same steps and prints the numbers quoted here; other positions,
altitudes and modulation levels are set from the command line (`--lat --lon --alt --phi --epoch`).

All fluxes are **omnidirectional**: the intensity integrated over the solid angle of the visible Earth,
F(E) = ∫ I dΩ, in ph cm⁻² s⁻¹ keV⁻¹ (the band integrals below are in ph cm⁻² s⁻¹). Step 6 also gives the flux
through a plane facing nadir.

## Step 1: choose the inputs

| Input | Value | Where it enters |
|---|---|---|
| Position | 41.9° N, 12.5° E (geocentric latitude, east longitude) | the cutoffs of the lines of sight |
| Altitude | 600 km, i.e., d = R_E + 600 = 6971 km from the centre of the Earth | the geometry |
| Solar cycle | solar minimum, φ = 0.40 GV | the CR grid (`phi_GV`) |
| Geomagnetic field | IGRF-14, epoch 2025.0 | `cutoff_igrf_2025.h5` |
| CXB | Gruber et al. (1999) | `albedo_cxb_reflection.h5`, group `gruber99` |
| Hadronic model | FTFP_BERT_HP | `albedo_cr_emission.h5`, `intensity_total[0]` |

**Solar minimum.** The modulation potential φ is on the scale of the Vos & Potgieter (2015) interstellar
spectrum, as the monthly series of Usoskin et al. (2017) and Väisänen et al. (2023). The annual means of the
last two minima (2009 and 2020) are 0.39–0.40 GV; the maxima of cycles 22–24 average 0.80 GV. For a given date,
take φ from that series. φ = 0.40 GV is a grid value; between grid values the tables are interpolated linearly.

```python
import albedo_earth as ae
em = ae.Emission("gruber99", "ftfp")      # reflection I(E, b) and CR grid I(phi, R_c, E, b)
cut = ae.Cutoff(2025.0)                   # IGRF vertical cutoffs
```

## Step 2: geometry

A line of sight at nadir angle θ has impact parameter b = d sin θ. Seen from 600 km, the solid Earth
(b < R_E) fills θ < 66.05° (3.733 sr), and the atmosphere up to 400 km (b < R_E + 400 km) fills θ < 76.24°
(4.789 sr). The solid angle of each impact-parameter bin of the tables is

  Ω_b = 2π [ (1 − b_lo²/d²)^½ − (1 − b_hi²/d²)^½ ].

```python
om = ae.solid_angle_per_b(em.b_edges, ae.R_E + 600.0)   # sr, one value per b bin (175)
```

## Step 3: CXB reflection

The reflection depends only on b (not on the position of the observer, and only marginally on the solar cycle),
so the flux is a single sum over the b bins:

  F_refl(E) = Σ_b I_refl(E, b) Ω_b.

```python
f_refl = em.refl @ om                     # ph cm-2 s-1 keV-1
```

| Band [keV] | 1–10 | 10–100 | 100–1000 | 1000–10 000 | 25–300 |
|---|---|---|---|---|---|
| F_refl [ph cm⁻² s⁻¹] | 0.1085 | 1.667 | 0.1244 | 0.000228 | 1.379 |

The statistical error of the reflection (`I_E_b_err`) is negligible: 0.04 % at 25–300 keV.

## Step 4: the cutoff rigidity of each line of sight

The CR emission depends on the vertical cutoff rigidity R_c where the cosmic rays hit the atmosphere. The
cutoff grid gives R_c = `r_eq_GV` / `L_bmin`² at 35 km altitude, interpolated bilinearly in latitude and
longitude:

```python
rc_nadir = cut.at(41.9, 12.5)             # 5.83 GV below the observer
```

The observer sees a large part of the globe: from 600 km the horizon is 2700 km away (great-circle distance
24°). Each line of sight is therefore assigned the cutoff of the point where it reaches 35 km, or of its tangent
point if it passes above 35 km (limb rays). Above Rome, the cutoffs over the visible disk range from 1.2 GV
(towards the north) to 13.0 GV (towards the south), and 0.6–14.6 GV including the limb.

`ae.lines_of_sight` divides each b bin of the tables into rings of equal width in b (4 by default) and each
ring into azimuth sectors (72), computes the exact solid angle of each ray and the point where its cutoff is
sampled:

```python
v = ae.lines_of_sight(41.9, 12.5, 600.0, em.b_edges)
```

## Step 5: cosmic-ray induced emission

**5a. Quick estimate: one cutoff for the whole visible Earth.** Interpolate the uniform-Earth table at the
nadir cutoff (linearly in φ and R_c) and sum over b as in Step 3:

  F_CR(E) ≈ Σ_b I_CR(E, b; φ, R_c,nadir) Ω_b.

```python
f_cr_nadir = em.cr_uniform(0.40, rc_nadir) @ om
```

**5b. Reference: the cutoff of each line of sight.** For each ray k, interpolate the table at the ray cutoff
R_c,k and at its b bin, and sum with the ray solid angles:

  F_CR(E) = Σ_k I_CR(E, b_k; φ, R_c,k) ΔΩ_k.

```python
i_refl, i_cr = ae.sky(em, v, cut, phi=0.40)   # [E, ray]; also the ray cutoffs in v["rc"]
f_cr = ae.flux(i_cr, v)
```

| Band [keV] | 1–10 | 10–100 | 100–1000 | 1000–10 000 | 25–300 |
|---|---|---|---|---|---|
| F_CR, 5a (nadir cutoff) | 0.00066 | 0.932 | 1.473 | 0.556 | 1.738 |
| F_CR, 5b (per line of sight) | 0.00066 | 0.965 | 1.534 | 0.593 | 1.801 |
| ratio 5a / 5b | 0.997 | 0.966 | 0.960 | 0.938 | 0.965 |

The nadir cutoff underestimates the emission by 3–6 % here, because the northern part of the view has lower
cutoffs and the limb, which dominates at MeV energies, reaches high magnetic latitudes. The difference depends
on the position (it changes sign where the cutoff gradient reverses); use 5b for quantitative work.

## Step 6: the total spectrum

  F(E) = F_refl(E) + F_CR(E).

```python
f_total = f_refl + f_cr                                   # omnidirectional
f_plane = ae.flux(i_refl + i_cr, v, cosine=True)          # through a plane facing nadir
```

| Band [keV] | 1–10 | 10–100 | 100–1000 | 1000–10 000 | 25–300 |
|---|---|---|---|---|---|
| Reflection | 0.1085 | 1.667 | 0.1244 | 0.000228 | 1.379 |
| CR-induced | 0.00066 | 0.965 | 1.534 | 0.593 | 1.801 |
| **Total** | **0.1092** | **2.632** | **1.658** | **0.593** | **3.179** |
| Total through a plane facing nadir | 0.0689 | 1.784 | 1.158 | 0.376 | 2.189 |
| CR fraction of the total | 0.006 | 0.37 | 0.93 | 1.00 | 0.57 |

The reflection dominates below about 50 keV, the CR emission above. The script writes the spectrum per energy
bin to `rome_600km_solar_min.txt` (columns: E_lo, E_hi [keV]; F_refl, F_CR, F_total, the hadronic minimum and
maximum of F_total, and F_total with the line-scaled variant, in ph cm⁻² s⁻¹ keV⁻¹) and, if matplotlib is
installed, plots E² F(E).

The angular distribution, e.g., for an instrument with a limited field of view, is in the same arrays: sum
`i_refl + i_cr` times `v["dom"]` over the rays inside the field of view, or use `ae.nadir_profile` for the
azimuth-averaged intensity against the nadir angle.

## Step 7: uncertainties

Repeat Steps 3–6 with the variants in the tables (`ae.Emission(model=...)`, `line_scaled=True`,
`cxb_model=...`). Total relative to the reference:

| Variant | 1–10 | 10–100 | 100–1000 | 1000–10 000 | 25–300 |
|---|---|---|---|---|---|
| QGSP_BIC_HP | 1.000 | 0.986 | 0.955 | 0.945 | 0.975 |
| QGSP_INCLXX_HP | 0.999 | 0.977 | 0.916 | 0.829 | 0.958 |
| Line-scaled nuclear lines | 1.000 | 1.000 | 0.972 | 1.007 | 1.000 |
| CXB of Churazov et al. (2007), +10 % | 1.099 | 1.063 | 1.008 | 1.000 | 1.043 |

- **Hadronic model:** the spread of the three models (2–4 % at 25–300 keV, up to 17 % at 1–10 MeV here). The
  hadronic systematic of the reference orbits is wider (models above 10 GeV/nucleon and the Bertini
  de-excitation, see the README): about −15 % at 25–300 keV at this cutoff (R_c ≈ 6 GV).
- **Nuclear lines:** the line-scaled variant lowers the 511 keV line by a factor of 2 (3 % of the
  100–1000 keV band); in the 1–10 MeV band the weaker 1.6–6.1 MeV lines and the stronger 4.44 MeV line nearly
  compensate (+0.7 %).
- **CXB normalisation:** the reflection scales with it (+10 % for `churazov07`).
- **Other systematic terms**, not tabulated per position: CR spectra (+2 to +14 %, increasing with R_c), field
  inside the atmosphere (up to −21 % on the MeV limb), cutoff model. The reference orbits list them per energy
  bin (`<orbit>/<phase>/shifts/`); use them as a guide to the size of each term.
- **Statistical:** below 1 % in broad bands for both components.

## Variations

- **Another date:** take φ for that date from the Väisänen et al. (2023) series and the cutoff grid of the
  nearest epoch; e.g., for mid-2013 (annual mean φ = 0.62 GV) `--phi 0.62 --epoch 2013.5`. At 25–300 keV the
  CR emission at solar minimum (0.40 GV) exceeds that at maximum (0.80 GV) by 6 % at R_c = 15 GV, 15 % at 6 GV
  and 48 % at R_c = 0; the reflection does not change.
- **Another altitude:** any altitude above 400 km (`--alt`). Below 400 km the observer is inside the model
  atmosphere and the tables do not apply.
- **An orbit average:** repeat Steps 4–6 at positions along the orbit and average the fluxes; the reference
  orbits (`albedo_reference_orbits.h5`) are precomputed examples (864 positions each).
