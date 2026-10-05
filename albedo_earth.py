"""Earth X/gamma-ray emission seen from an observer, from the released tables (numpy + h5py only).

Total emission = CXB reflection + cosmic-ray (CR) induced emission. Both are tabulated as emergent intensities
I(E, b) [ph cm-2 s-1 sr-1 keV-1] against the photon energy E and the impact parameter b of the ray; the CR emission
is given for a whole Earth with a uniform vertical cutoff R_c, on a grid of (phi, R_c). For an observer above the
atmosphere, each line of sight has b = d sin(theta) (d: distance from the centre of the Earth, theta: nadir angle);
the CR intensity of the ray is taken at the cutoff where the ray reaches the emission layer (35 km), or at its
tangent point for rays passing above it. See docs/worked_example.md for the derivation, step by step.

    import albedo_earth as ae
    em = ae.Emission()                                  # Gruber CXB, modern CR spectra, FTFP_BERT_HP
    cut = ae.Cutoff(2025.0)                             # IGRF-14 vertical cutoffs
    v = ae.lines_of_sight(41.9, 12.5, 600.0, em.b_edges)
    i_refl, i_cr = ae.sky(em, v, cut, phi=0.40)        # [E, ray] intensities
    f_total = ae.flux(i_refl + i_cr, v)                 # ph cm-2 s-1 keV-1

The response functions (data/responses/) give I(E, b) for other input spectra: see :class:`Response` and
examples/fold_responses.py.
"""

from __future__ import annotations

import json
import os

import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
R_E = 6371.0                 # km, Earth radius of the model
R_TOP = R_E + 400.0          # km, top of the model atmosphere
H_EM = 35.0                  # km, CR emission layer, where the cutoff of a ray is sampled


class Emission:
    """Uniform-Earth tables: reflection I(E, b) and the CR grid I(phi, R_c, E, b) of one hadronic model.

    cxb_model: gruber99 (reference), churazov07, ajello08, turler10, frontera07.
    cr_file: albedo_cr_emission.h5 (reference) or albedo_cr_emission_s07.h5.
    model: ftfp (reference), qgsp, inclxx; line_scaled=True uses the line-scaled variant (FTFP_BERT_HP only).
    """

    def __init__(self, cxb_model: str = "gruber99", model: str = "ftfp", line_scaled: bool = False,
                 cr_file: str = "albedo_cr_emission.h5", data: str = DATA):
        import h5py

        with h5py.File(os.path.join(data, cr_file)) as h:
            self.e_edges, self.b_edges = h["e_edges"][:], h["b_edges"][:]
            self.phis, self.rcs = h["phi_GV"][:], h["rc_GV"][:]
            models = json.loads(h.attrs["models"])
            if line_scaled:
                if model != "ftfp":
                    raise ValueError("the line-scaled variant exists for FTFP_BERT_HP only")
                self.cr = h["intensity_total_line_scaled"][:].astype(float)
            else:
                self.cr = h["intensity_total"][models.index(model)].astype(float)    # [phi, rc, E, b]
        with h5py.File(os.path.join(data, "albedo_cxb_reflection.h5")) as h:
            if not (np.array_equal(h["e_edges"][:], self.e_edges) and np.array_equal(h["b_edges"][:], self.b_edges)):
                raise ValueError("CXB and CR tables have different binning")
            self.refl = h[cxb_model]["I_E_b"][:]                                      # [E, b]
            self.refl_err = h[cxb_model]["I_E_b_err"][:]
        self.cxb_model, self.model = cxb_model, model

    @property
    def e_centres(self):
        return np.sqrt(self.e_edges[1:] * self.e_edges[:-1])

    def cr_at_phi(self, phi: float):
        """CR intensity [rc, E, b] at modulation phi [GV], linear interpolation between grid values."""
        if not self.phis[0] - 1e-9 <= phi <= self.phis[-1] + 1e-9:
            raise ValueError(f"phi = {phi} GV outside the grid {self.phis}")
        j = int(np.clip(np.searchsorted(self.phis, phi) - 1, 0, len(self.phis) - 2))
        t = (phi - self.phis[j]) / (self.phis[j + 1] - self.phis[j])
        return (1 - t) * self.cr[j] + t * self.cr[j + 1]

    def cr_at(self, phi: float, rc):
        """CR intensity at modulation phi and cutoff(s) rc [GV] (linear in both): [E, b] for a scalar rc, or
        [E, n] for an array of n cutoffs with a matching array of b-bin indices (see :func:`sky`)."""
        g = self.cr_at_phi(phi)
        rc = np.clip(rc, self.rcs[0], self.rcs[-1])
        j = np.clip(np.searchsorted(self.rcs, rc, side="right") - 1, 0, len(self.rcs) - 2)
        t = (rc - self.rcs[j]) / (self.rcs[j + 1] - self.rcs[j])
        return g, j, t

    def cr_uniform(self, phi: float, rc: float):
        """CR intensity [E, b] for a whole Earth with one cutoff rc [GV]."""
        g, j, t = self.cr_at(phi, float(rc))
        return (1 - t) * g[j] + t * g[j + 1]


class Cutoff:
    """Vertical cutoff rigidity [GV] from the IGRF grid of an epoch: R_c = r_eq / L^2 (bilinear interpolation)."""

    def __init__(self, epoch: float = 2025.0, l_def: str = "bmin", data: str = DATA):
        import h5py

        with h5py.File(os.path.join(data, f"cutoff_igrf_{epoch:g}.h5")) as h:
            self.lat, self.lon = h["lat_deg"][:], h["lon_deg"][:]
            self.rc = float(h.attrs["r_eq_GV"]) / h["L_bmin" if l_def == "bmin" else "L_apex"][:] ** 2
        self.epoch = epoch

    def at(self, lat_deg, lon_deg):
        """R_c [GV] at geocentric latitude and east longitude [deg] (arrays broadcast)."""
        lat, lon = np.asarray(lat_deg, float), np.asarray(lon_deg, float) % 360
        fi = np.clip((lat - self.lat[0]) / (self.lat[1] - self.lat[0]), 0, len(self.lat) - 1 - 1e-9)
        i = np.floor(fi).astype(int)
        ti = fi - i
        i1 = np.minimum(i + 1, len(self.lat) - 1)
        fj = (lon - self.lon[0]) / (self.lon[1] - self.lon[0])
        j = np.floor(fj).astype(int) % len(self.lon)
        tj = fj - np.floor(fj)
        j1 = (j + 1) % len(self.lon)
        g = self.rc
        return (1 - ti) * ((1 - tj) * g[i, j] + tj * g[i, j1]) + ti * ((1 - tj) * g[i1, j] + tj * g[i1, j1])

    def __call__(self, p):
        """R_c [GV] at unit position vectors p (..., 3), Earth-fixed frame (z to the north pole, x to 0 deg E)."""
        p = np.asarray(p, float)
        return self.at(np.degrees(np.arcsin(np.clip(p[..., 2], -1, 1))), np.degrees(np.arctan2(p[..., 1], p[..., 0])))


KM_TO_CM = 1e5


class Step:
    """Piecewise-constant intensity [cm-2 s-1 sr-1 keV-1] on bin edges [keV], with exact bin integrals."""

    def __init__(self, edges, level):
        self.edges, self.level = np.asarray(edges, float), np.asarray(level, float)
        self._cum = np.concatenate([[0.0], np.cumsum(self.level * np.diff(self.edges))])

    def __call__(self, e_kev):
        e = np.asarray(e_kev, float)
        i = np.searchsorted(self.edges, e, side="right") - 1
        ok = (i >= 0) & (i < len(self.level))
        return np.where(ok, self.level[np.clip(i, 0, len(self.level) - 1)], 0.0)

    def bin_integrals(self, edges):
        x = np.clip(np.asarray(edges, float), self.edges[0], self.edges[-1])
        i = np.clip(np.searchsorted(self.edges, x, side="right") - 1, 0, len(self.level) - 1)
        return np.maximum(np.diff(self._cum[i] + self.level[i] * (x - self.edges[i])), 0.0)


class Response:
    """Response functions of one primary species: 'cxb', 'proton', 'alpha' (4He), 'e-', 'e+'.

    yield[c, i, j, k]: escaping photons in (E_j, b_k) per primary of energy bin i entering the source sphere of
    radius R_src isotropically (c = 0 aimed at the Earth disk, 1 grazing). For an isotropic primary intensity
    I(E0) [cm-2 s-1 sr-1 keV-1; E0 = total kinetic energy for nuclei], the rate per (E, b) bin is
    sum_i pi 4 pi R_src^2 (int_i I dE0) yield[i], and I(E, b) = rate / (dE 8 pi^2 b db)."""

    def __init__(self, species: str, data: str = DATA):
        import h5py

        with h5py.File(os.path.join(data, "responses", f"response_{species}.h5")) as h:
            self.e0_edges, self.e_edges, self.b_edges = h["e0_edges"][:], h["e_edges"][:], h["b_edges"][:]
            self.yield_ = h["yield"][:].astype(float)                      # [class, E0, E, b]
            self.yield_var = h["yield_var"][:].astype(float) if "yield_var" in h else None
            self.lepton_yield = h["lepton_yield"][:].astype(float) if "lepton_yield" in h else None
            self.lep_edges = h["lep_edges_keV"][:] if "lep_edges_keV" in h else None
            self.r_src_km = float(h.attrs["R_src_km"])
        self.species = species

    def weights(self, intensity):
        """Primary rate [s-1] per primary-energy bin: pi 4 pi R_src^2 int_i I dE0. ``intensity``: a callable of
        E0 [keV] (integrated with 40 log steps per bin) or an object with ``bin_integrals(edges)``."""
        e = self.e0_edges
        if hasattr(intensity, "bin_integrals"):
            ints = intensity.bin_integrals(e)
        else:
            x = np.logspace(np.log10(e[:-1]), np.log10(e[1:]), 41, axis=1)
            y = intensity(x)
            ints = (np.diff(x, axis=1) * 0.5 * (y[:, 1:] + y[:, :-1])).sum(1)
        return np.pi * 4 * np.pi * (self.r_src_km * KM_TO_CM) ** 2 * ints

    def _norm(self):
        b = self.b_edges * KM_TO_CM
        return np.diff(self.e_edges)[:, None] * (8 * np.pi ** 2 * 0.5 * np.diff(b ** 2))[None, :]

    def rate(self, intensity, variant: str = "all"):
        """Escaping photon rate [s-1] per (E, b) bin; variant 'all' or 'disk' (primaries aimed at the disk)."""
        y = self.yield_[0] if variant == "disk" else self.yield_.sum(0)
        return np.einsum("i,ijk->jk", self.weights(intensity), y)

    def intensity(self, intensity, variant: str = "all"):
        """Emergent intensity I(E, b) [ph cm-2 s-1 sr-1 keV-1]."""
        return self.rate(intensity, variant) / self._norm()

    def intensity_err(self, intensity):
        """Statistical error of :meth:`intensity` (CXB only: Poisson, photons of different primaries independent)."""
        if self.yield_var is None:
            raise ValueError("no Poisson variance for showers: use the event replicas (yield_replicas)")
        w = self.weights(intensity)
        return np.sqrt(np.einsum("i,ijk->jk", w ** 2, self.yield_var.sum(0))) / self._norm()

    def lepton_rate(self, intensity):
        """Escaping e- (0) and e+ (1) rate [s-1] per bin of ``lep_edges`` (CR files)."""
        return np.einsum("i,sil->sl", self.weights(intensity), self.lepton_yield)


def reentrant_intensities(lepton_rate, lep_edges, r_src_km: float):
    """Re-entrant e- and e+ as isotropic intensities on the source sphere (Step objects): the escaping lepton rate
    re-injected over the whole Earth, I = rate / (pi 4 pi R_src^2 dE)."""
    geom = np.pi * 4 * np.pi * (r_src_km * KM_TO_CM) ** 2
    return [Step(lep_edges, np.asarray(r) / (geom * np.diff(lep_edges))) for r in lepton_rate]


def unit(lat_deg, lon_deg):
    """Unit position vector at geocentric latitude and east longitude [deg]."""
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def solid_angle_per_b(b_edges, d_km: float):
    """Solid angle [sr] of each impact-parameter bin seen from distance d [km]: Omega(<b) = 2 pi (1 - cos theta),
    sin theta = b / d."""
    if d_km <= b_edges[-1]:
        raise ValueError("the observer must be above the top of the atmosphere (R_E + 400 km)")
    return np.diff(2 * np.pi * (1 - np.sqrt(1 - (np.asarray(b_edges) / d_km) ** 2)))


def lines_of_sight(lat_deg: float, lon_deg: float, alt_km: float, b_edges, n_sub: int = 4, n_psi: int = 72):
    """Lines of sight to the Earth from an observer at altitude alt_km above (lat, lon).

    Each impact-parameter bin of the tables is split into n_sub rings of equal width in b, each ring into n_psi
    azimuth sectors; ring solid angles are exact, so the sum over the rays of a b bin is its solid angle. Returns a
    dict of flat arrays over the rays: kb (b-bin index), b [km], theta (nadir angle, rad), psi (azimuth from east
    towards north, rad), dom (solid angle, sr), loc (unit vector where the cutoff is sampled: the 35 km crossing,
    or the tangent point of rays passing above 35 km), hit (True if the ray reaches 35 km)."""
    d = R_E + alt_km
    b_edges = np.asarray(b_edges, float)
    if d <= b_edges[-1]:
        raise ValueError("the observer must be above the top of the atmosphere (R_E + 400 km)")
    nb = len(b_edges) - 1
    be = np.concatenate([np.linspace(b_edges[k], b_edges[k + 1], n_sub + 1)[:-1] for k in range(nb)] + [b_edges[-1:]])
    kb_ring = np.repeat(np.arange(nb), n_sub)
    cos_e = np.sqrt(1 - (be / d) ** 2)
    dom_ring = 2 * np.pi * (cos_e[:-1] - cos_e[1:])
    b_ring = 0.5 * (be[1:] + be[:-1])
    th = np.arcsin(b_ring / d)
    psi = (np.arange(n_psi) + 0.5) * 2 * np.pi / n_psi
    p_obs = unit(lat_deg, lon_deg)
    o, down = d * p_obs, -p_obs
    east = np.cross([0, 0, 1.0], down)
    east = east / np.linalg.norm(east) if np.linalg.norm(east) > 1e-12 else np.array([1.0, 0, 0])
    north = np.cross(down, east)
    T, P = np.meshgrid(th, psi, indexing="ij")
    u = (np.cos(T)[..., None] * down
         + np.sin(T)[..., None] * (np.cos(P)[..., None] * east + np.sin(P)[..., None] * north))
    ou = u @ o
    rh = R_E + H_EM
    disc = ou ** 2 - (d ** 2 - rh ** 2)
    hit = disc > 0
    t = np.where(hit, -ou - np.sqrt(np.where(hit, disc, 0.0)), -ou)        # 35 km crossing, else tangent point
    q = o + t[..., None] * u
    loc = q / np.linalg.norm(q, axis=-1, keepdims=True)
    rays = dict(kb=np.repeat(kb_ring, n_psi), b=np.repeat(b_ring, n_psi), theta=T.ravel(), psi=P.ravel(),
                dom=np.repeat(dom_ring / n_psi, n_psi), loc=loc.reshape(-1, 3), hit=hit.ravel())
    rays.update(d_km=d, lat_deg=lat_deg, lon_deg=lon_deg, alt_km=alt_km)
    return rays


def sky(em: Emission, v: dict, cutoff, phi: float):
    """Reflection and CR intensities [ph cm-2 s-1 sr-1 keV-1] per (E, ray); also stores the ray cutoffs in v['rc']."""
    kb = v["kb"]
    rc = cutoff(v["loc"])
    g, j, t = em.cr_at(phi, rc)
    i_cr = (1 - t)[None] * g[j, :, kb].T + t[None] * g[j + 1, :, kb].T
    v["rc"] = rc
    return em.refl[:, kb], i_cr


def flux(i, v, cosine: bool = False):
    """Flux [ph cm-2 s-1 keV-1]: sum over the rays of I dOmega (omnidirectional), or of I cos(theta) dOmega
    (through a plane facing nadir) with cosine=True."""
    w = v["dom"] * (np.cos(v["theta"]) if cosine else 1.0)
    return i @ w


def band(e_edges, y, lo: float, hi: float):
    """Integrate a per-keV quantity (energy on axis 0) over the bins inside [lo, hi] keV."""
    e = np.asarray(e_edges)
    m = (e[:-1] >= lo * (1 - 1e-9)) & (e[1:] <= hi * (1 + 1e-9))
    return np.tensordot(np.diff(e)[m], np.asarray(y)[m], axes=(0, 0))


def nadir_profile(i, v):
    """Azimuth-averaged intensity per impact-parameter ring: nadir angle [deg], mean I (E, ring)."""
    th = np.unique(v["theta"])
    k = np.searchsorted(th, v["theta"])
    w = np.zeros((i.shape[0], len(th)))
    om = np.bincount(k, weights=v["dom"], minlength=len(th))
    for r in range(i.shape[0]):
        w[r] = np.bincount(k, weights=i[r] * v["dom"], minlength=len(th))
    return np.degrees(th), w / om
