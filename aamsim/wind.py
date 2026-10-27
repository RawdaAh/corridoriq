"""Spatiotemporal wind fields and the wind-to-link projection of the TRR paper.

A :class:`WindField` stores the horizontal wind vector (u east/along-corridor,
v north/lateral, both mph) at every corridor node for every time bin.

Direction convention: ``direction_deg`` is the direction the wind blows
*toward*, counter-clockwise from the corridor axis (+x). Hence 0 deg is a pure
tailwind for a flight from the entry face to the exit face and 180 deg is a pure
headwind; this matches the experiment definitions in the TRR paper.

Link projection (paper Eqs. 14-25): nodal components are averaged on the link,
the relative angle gamma between link heading and wind is folded into
[0, 180], and the along-track component psi = w cos(gamma) changes the ground
speed.  ``model='paper'`` reproduces the paper's piecewise rule (tail if
gamma <= 45, head if gamma >= 135, otherwise no effect); ``model='continuous'``
uses v = v_uav + psi for every link.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .corridor import Corridor

MIN_GROUND_SPEED_MPH = 1.0  # links slower than this are treated as impassable


@dataclass
class WindField:
    u: np.ndarray          # (B, N) along-corridor component, mph
    v: np.ndarray          # (B, N) lateral component, mph
    bin_sec: float = 180.0  # duration of one wind time bin
    name: str = "wind"

    @property
    def n_bins(self) -> int:
        return self.u.shape[0]

    def bin_at(self, t_sec):
        b = np.floor_divide(np.asarray(t_sec, dtype=float), self.bin_sec).astype(int)
        return np.clip(b, 0, self.n_bins - 1)

    def speed(self) -> np.ndarray:
        return np.hypot(self.u, self.v)

    def scaled(self, factor: float, name: str | None = None) -> "WindField":
        return WindField(self.u * factor, self.v * factor, self.bin_sec, name or self.name)

    def renamed(self, name: str) -> "WindField":
        return WindField(self.u, self.v, self.bin_sec, name)

    # --------------------------------------------------------------- link factors
    def link_ground_speed(self, corr: Corridor, v_uav_mph: float, model: str = "paper"):
        """Wind-adjusted ground speed (mph) for every (bin, arc). Shape (B, A)."""
        d = corr.pos[corr.arc_to] - corr.pos[corr.arc_from]
        horiz = np.hypot(d[:, 0], d[:, 1])
        theta = np.degrees(np.arctan2(d[:, 1], d[:, 0])) % 360.0          # Eq. 14
        ub = 0.5 * (self.u[:, corr.arc_from] + self.u[:, corr.arc_to])    # Eq. 16
        vb = 0.5 * (self.v[:, corr.arc_from] + self.v[:, corr.arc_to])
        w = np.hypot(ub, vb)                                              # Eq. 17a
        alpha = np.degrees(np.arctan2(vb, ub)) % 360.0                    # Eq. 17b
        gamma = (alpha - theta[None, :]) % 360.0                          # Eq. 18
        gamma = np.where(gamma > 180.0, 360.0 - gamma, gamma)
        psi = w * np.cos(np.radians(gamma))                               # Eq. 19
        if model == "paper":                                              # Eq. 20
            vg = np.where(gamma <= 45.0, v_uav_mph + np.abs(psi),
                          np.where(gamma >= 135.0, v_uav_mph - np.abs(psi), v_uav_mph))
        elif model == "continuous":
            vg = v_uav_mph + psi
        elif model == "crab":
            # hold the track: point into the cross component, then add the along component
            cross = w * np.sin(np.radians(gamma))
            vg = np.sqrt(np.maximum(v_uav_mph ** 2 - cross ** 2, 0.0)) + psi
        else:
            raise ValueError("model must be 'paper', 'continuous' or 'crab'")
        # Purely vertical links have no horizontal track: wind does not act on them.
        vg = np.where(horiz[None, :] < 1e-9, v_uav_mph, vg)
        return vg

    def link_factors(self, corr: Corridor, v_uav_mph: float, model: str = "paper"):
        """Return (lambda, beta) of Eqs. 24-25 for every (bin, arc).

        lambda = v_uav / v_ij multiplies travel time, beta = lambda**2 multiplies
        energy. Impassable links (ground speed below 1 mph) get +inf.
        """
        vg = self.link_ground_speed(corr, v_uav_mph, model)
        with np.errstate(divide="ignore"):
            lam = np.where(vg >= MIN_GROUND_SPEED_MPH, v_uav_mph / vg, np.inf)
        return lam, lam ** 2


# ------------------------------------------------------------------- constructors

def _uv(speed, direction_deg):
    r = np.radians(direction_deg)
    return speed * np.cos(r), speed * np.sin(r)


def zero_wind(corr: Corridor, n_bins: int = 1, bin_sec: float = 180.0) -> WindField:
    z = np.zeros((n_bins, corr.n_nodes))
    return WindField(z, z.copy(), bin_sec, "zero")


def from_speed_function(corr: Corridor, speed_fn, direction_deg, n_bins: int = 1,
                        bin_sec: float = 180.0, name: str = "wind") -> WindField:
    """Build a field from ``speed_fn(x, y, z, t_sec) -> mph``.

    ``direction_deg`` may be a scalar or a function with the same signature.
    """
    x, y, z = corr.pos.T
    u = np.zeros((n_bins, corr.n_nodes))
    v = np.zeros_like(u)
    for b in range(n_bins):
        t = b * bin_sec
        s = np.broadcast_to(np.asarray(speed_fn(x, y, z, t), dtype=float), x.shape)
        d = direction_deg(x, y, z, t) if callable(direction_deg) else direction_deg
        u[b], v[b] = _uv(s, d)
    return WindField(u, v, bin_sec, name)


def trr_shear_profile(z):
    """Altitude-varying speeds of TRR Sec. 6.2 (mph): 10 / 15 / 20 / 30."""
    z = np.asarray(z, dtype=float)
    return np.select([z < 200, z < 300, z < 400], [10.0, 15.0, 20.0], 30.0)


def trr_weak_profile(z):
    """Weaker profile of TRR Sec. 6.3, blocks B1-B2 (mph): 10 / 12 / 13 / 15."""
    z = np.asarray(z, dtype=float)
    return np.select([z < 200, z < 300, z < 400], [10.0, 12.0, 13.0], 15.0)


def trr_altitude_wind(corr, direction_deg):
    """TRR experiment 6.2: shear profile, uniform along the corridor."""
    return from_speed_function(corr, lambda x, y, z, t: trr_shear_profile(z),
                               direction_deg, name=f"shear_{direction_deg:g}")


def trr_spatial_wind(corr, direction_deg, n_blocks: int = 5, weak_blocks=(0, 1)):
    """TRR experiment 6.3: weak profile in B1-B2, strong profile in B3-B5."""
    length = corr.pos[:, 0].max()

    def speed(x, y, z, t):
        block = np.minimum((x / length * n_blocks).astype(int), n_blocks - 1)
        # a node on a block boundary belongs to the downstream block, except the last
        weak = np.isin(block, weak_blocks)
        return np.where(weak, trr_weak_profile(z), trr_shear_profile(z))
    return from_speed_function(corr, speed, direction_deg,
                               name=f"spatial_{direction_deg:g}")


def trr_temporal_wind(corr, direction_deg, trend: str, n_bins: int = 20,
                      bin_sec: float = 180.0, window_bins: int = 10,
                      lo: float = 2.0, hi: float = 22.0):
    """TRR experiment 6.4: speed changes linearly over the 30-min dispatch window
    (10 x 3-min intervals) and then stays at its final value until the horizon."""
    def speed(x, y, z, t):
        frac = min(t / (bin_sec * (window_bins - 1)), 1.0)
        s = lo + (hi - lo) * frac if trend == "calm_to_strong" else hi - (hi - lo) * frac
        return np.full_like(x, s)
    return from_speed_function(corr, speed, direction_deg, n_bins, bin_sec,
                               name=f"{trend}_{direction_deg:g}")


# --------------------------------------------- synthetic "truth" and planning beliefs

def synthetic_urban_wind(corr: Corridor, seed: int, ref_speed_mph: float = 15.0,
                         direction_deg: float = 0.0, n_bins: int = 20,
                         bin_sec: float = 180.0, heterogeneity: float = 0.35,
                         veer_deg: float = 25.0, drift: float = 0.3) -> WindField:
    """A reproducible stand-in for a high-fidelity (e.g. CFD or reanalysis) wind field.

    * vertical shear: power law  s(z) = s_ref (z / 400 ft)^0.25 (ground uses 50 ft)
    * horizontal heterogeneity: smooth random bumps along x and y (building wakes)
    * directional veer: smooth random rotation along the corridor
    * temporal drift: linear speed trend plus a slow oscillation over the horizon

    This is synthetic and only meant to give the four planning beliefs something
    structured to be wrong about. Replace it with real data via
    :func:`wind_from_arrays` for a paper result.
    """
    rng = np.random.default_rng(seed)
    x, y, z = corr.pos.T
    L = max(x.max(), 1.0)
    W = max(y.max(), 1.0)
    zz = np.maximum(z, 50.0)
    shear = (zz / 400.0) ** 0.25

    # smooth spatial multiplier from a few Gaussian bumps
    n_bump = 6
    cx = rng.uniform(0, L, n_bump)
    cy = rng.uniform(0, W, n_bump)
    amp = rng.normal(0, heterogeneity, n_bump)
    sx = rng.uniform(0.08, 0.2, n_bump) * L
    low_alt_weight = np.clip(1.2 - z / 400.0, 0.2, 1.2)  # wakes are stronger near the ground
    bump = np.zeros_like(x)
    for k in range(n_bump):
        bump += amp[k] * np.exp(-0.5 * ((x - cx[k]) / sx[k]) ** 2
                                - 0.5 * ((y - cy[k]) / (0.6 * W + 1)) ** 2)
    spatial = np.clip(1.0 + bump * low_alt_weight, 0.3, 1.8)

    phase = rng.uniform(0, 2 * np.pi, 2)
    veer = veer_deg * np.sin(2 * np.pi * x / L + phase[0]) * low_alt_weight
    trend = rng.uniform(-drift, drift)
    u = np.zeros((n_bins, corr.n_nodes))
    v = np.zeros_like(u)
    for b in range(n_bins):
        f = b / max(n_bins - 1, 1)
        temporal = 1.0 + trend * (f - 0.5) + 0.1 * np.sin(2 * np.pi * f + phase[1])
        s = ref_speed_mph * shear * spatial * temporal
        u[b], v[b] = _uv(s, direction_deg + veer)
    return WindField(u, v, bin_sec, "ground_truth")


def profile_belief(truth: WindField, corr: Corridor) -> WindField:
    """Height-only belief: the vector mean of the truth over each altitude layer,
    applied uniformly over the layer (analogous to GeoWind2Plan's 'Profile wind').
    Uses the first time bin only, i.e. a single sounding at planning time."""
    k = corr.ijk[:, 2]
    u = np.zeros_like(truth.u)
    v = np.zeros_like(truth.v)
    for layer in np.unique(k):
        m = k == layer
        u[:, m] = truth.u[0, m].mean()
        v[:, m] = truth.v[0, m].mean()
    return WindField(u, v, truth.bin_sec, "profile")


def forecast_belief(truth: WindField, corr: Corridor, seed: int,
                    speed_bias: float = 0.15, dir_error_deg: float = 15.0,
                    noise: float = 0.15, lag_bins: int = 2) -> WindField:
    """A forecast-like belief: the truth with a speed bias, a direction error,
    smooth spatial noise and a timing lag (the forecast is 'late' by ``lag_bins``)."""
    rng = np.random.default_rng(seed + 7919)
    bias = 1.0 + rng.uniform(-speed_bias, speed_bias)
    rot = np.radians(rng.normal(0, dir_error_deg))
    x = corr.pos[:, 0]
    L = max(x.max(), 1.0)
    ph = rng.uniform(0, 2 * np.pi)
    mult = bias * (1.0 + noise * np.sin(2 * np.pi * x / (0.7 * L) + ph))
    idx = np.clip(np.arange(truth.n_bins) - lag_bins, 0, truth.n_bins - 1)
    u0, v0 = truth.u[idx], truth.v[idx]
    u = (np.cos(rot) * u0 - np.sin(rot) * v0) * mult
    v = (np.sin(rot) * u0 + np.cos(rot) * v0) * mult
    return WindField(u, v, truth.bin_sec, "forecast")


def planning_beliefs(truth: WindField, corr: Corridor, seed: int) -> dict:
    """The four planning winds, mirroring GeoWind2Plan's four methods."""
    return {
        "ground_truth": truth.renamed("ground_truth"),
        "forecast": forecast_belief(truth, corr, seed),
        "profile": profile_belief(truth, corr),
        "wind_agnostic": zero_wind(corr, truth.n_bins, truth.bin_sec).renamed("wind_agnostic"),
    }


def wind_from_arrays(u, v, bin_sec: float, name: str = "data") -> WindField:
    """Wrap externally supplied node wind (e.g. interpolated HRRR/CFD) of shape (B, N)."""
    u = np.atleast_2d(np.asarray(u, dtype=float))
    v = np.atleast_2d(np.asarray(v, dtype=float))
    if u.shape != v.shape:
        raise ValueError("u and v must have the same shape (bins, nodes)")
    return WindField(u, v, bin_sec, name)
