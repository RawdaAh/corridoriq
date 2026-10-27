"""Fast building-aware urban wind surrogate (NOT CFD).

For every altitude layer z the buildings taller than z are solid cells on a
horizontal grid. The flow on that layer is

    1. potential flow around the solid cells: Laplace's equation for the velocity
       potential with the free-stream value on the domain boundary and zero normal
       gradient at building walls (red-black SOR). Flow accelerates through gaps
       and around corners and turns around blocks.
    2. a wake deficit: speed is reduced behind buildings, decaying with the
       upstream distance to the nearest solid cell (length scale ``wake_len``).
    3. a power-law shear profile of the background wind, s(z) = s_ref (z/120 m)^alpha.

This keeps the qualitative features GeoWind2Plan shows (canyon speed-up, wakes,
turning) at millisecond cost. ``demo/geo_engine.js`` implements the same
algorithm in the browser; tests/test_geo.py checks the two agree.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

MPH_PER_MS = 2.236936


@dataclass
class WindParams:
    speed_ms: float = 4.0          # background speed at the reference height
    direction_deg: float = 0.0     # direction the wind blows TOWARD, CCW from the corridor axis
    ref_height: float = 120.0
    alpha: float = 0.25            # shear exponent (urban terrain)
    iterations: int = 300
    omega: float = 1.85
    wake_len: float = 90.0         # m
    wake_depth: float = 0.65       # max fractional speed deficit just behind a wall
    wake_reach: int = 24           # cells sampled upstream
    max_speedup: float = 1.8


def layer_speed(p: WindParams, z: float) -> float:
    return p.speed_ms * (max(z, 10.0) / p.ref_height) ** p.alpha


def solve_layer(solid: np.ndarray, res: float, x0: float, y0: float, ux: float, uy: float,
                p: WindParams):
    """Velocity (u, v) on one layer. ``solid`` is (nx, ny) bool."""
    nx, ny = solid.shape
    gx = x0 + (np.arange(nx) + 0.5) * res
    gy = y0 + (np.arange(ny) + 0.5) * res
    X, Y = np.meshgrid(gx, gy, indexing="ij")
    phi = ux * X + uy * Y
    fluid = ~solid
    inner = np.zeros_like(fluid)
    inner[1:-1, 1:-1] = fluid[1:-1, 1:-1]
    # neighbour availability (solid neighbour -> mirror, i.e. use own value)
    fE = np.zeros_like(fluid); fE[:-1] = fluid[1:]
    fW = np.zeros_like(fluid); fW[1:] = fluid[:-1]
    fN = np.zeros_like(fluid); fN[:, :-1] = fluid[:, 1:]
    fS = np.zeros_like(fluid); fS[:, 1:] = fluid[:, :-1]
    I, J = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    colour = (I + J) % 2
    masks = [inner & (colour == 0), inner & (colour == 1)]
    for _ in range(p.iterations):
        for m in masks:
            E = np.empty_like(phi); E[:-1] = phi[1:]; E[-1] = phi[-1]
            W = np.empty_like(phi); W[1:] = phi[:-1]; W[0] = phi[0]
            N = np.empty_like(phi); N[:, :-1] = phi[:, 1:]; N[:, -1] = phi[:, -1]
            S = np.empty_like(phi); S[:, 1:] = phi[:, :-1]; S[:, 0] = phi[:, 0]
            E = np.where(fE, E, phi); W = np.where(fW, W, phi)
            N = np.where(fN, N, phi); S = np.where(fS, S, phi)
            avg = 0.25 * (E + W + N + S)
            phi = np.where(m, phi + p.omega * (avg - phi), phi)
    # gradient: central where both neighbours are fluid, one-sided otherwise
    pE = np.empty_like(phi); pE[:-1] = phi[1:]; pE[-1] = phi[-1]
    pW = np.empty_like(phi); pW[1:] = phi[:-1]; pW[0] = phi[0]
    pN = np.empty_like(phi); pN[:, :-1] = phi[:, 1:]; pN[:, -1] = phi[:, -1]
    pS = np.empty_like(phi); pS[:, 1:] = phi[:, :-1]; pS[:, 0] = phi[:, 0]
    bE = fE.copy(); bE[-1] = False
    bW = fW.copy(); bW[0] = False
    bN = fN.copy(); bN[:, -1] = False
    bS = fS.copy(); bS[:, 0] = False
    u = np.where(bE & bW, (pE - pW) / (2 * res),
                 np.where(bE, (pE - phi) / res, np.where(bW, (phi - pW) / res, 0.0)))
    v = np.where(bN & bS, (pN - pS) / (2 * res),
                 np.where(bN, (pN - phi) / res, np.where(bS, (phi - pS) / res, 0.0)))
    # domain-edge cells: free stream
    edge = np.zeros_like(fluid); edge[0] = edge[-1] = True; edge[:, 0] = edge[:, -1] = True
    u = np.where(edge & fluid, ux, u); v = np.where(edge & fluid, uy, v)
    # wake deficit from the nearest upstream solid cell
    U = math.hypot(ux, uy)
    if U > 1e-9 and solid.any():
        dx, dy = ux / U, uy / U
        shelter = np.zeros_like(phi)
        for k in range(1, p.wake_reach + 1):
            si = np.floor(I - dx * k + 0.5).astype(int)   # round half up (matches JS)
            sj = np.floor(J - dy * k + 0.5).astype(int)
            ok = (si >= 0) & (si < nx) & (sj >= 0) & (sj < ny)
            hit = np.zeros_like(fluid)
            hit[ok] = solid[si[ok], sj[ok]]
            shelter = np.maximum(shelter, hit * math.exp(-k * res / p.wake_len))
        # one 3x3 box-blur pass to soften the wake edges
        pad = np.pad(shelter, 1, mode="edge")
        shelter = sum(pad[1 + a:nx + 1 + a, 1 + b:ny + 1 + b] for a in (-1, 0, 1) for b in (-1, 0, 1)) / 9.0
        f = 1.0 - p.wake_depth * shelter
        u, v = u * f, v * f
    sp = np.hypot(u, v)
    cap = p.max_speedup * U
    scale = np.where(sp > cap, cap / np.maximum(sp, 1e-12), 1.0)
    u, v = u * scale * fluid, v * scale * fluid
    return u, v


@dataclass
class UrbanWind:
    """Layered wind on the corridor grid. u, v: (n_layers, nx, ny) in m/s."""
    z: np.ndarray
    u: np.ndarray
    v: np.ndarray
    x0: float
    y0: float
    res: float

    def sample(self, x, y, z):
        """Bilinear in x, y and linear in z. Returns (u, v) m/s."""
        x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
        nl, nx, ny = self.u.shape
        fx = np.clip((x - self.x0) / self.res - 0.5, 0, nx - 1.000001)
        fy = np.clip((y - self.y0) / self.res - 0.5, 0, ny - 1.000001)
        i, j = fx.astype(int), fy.astype(int)
        ax, ay = fx - i, fy - j
        k = np.clip(np.searchsorted(self.z, z) - 1, 0, nl - 2) if nl > 1 else np.zeros_like(i)
        az = np.clip((z - self.z[k]) / (self.z[k + 1] - self.z[k]), 0, 1) if nl > 1 else 0.0

        def bil(F, kk):
            return ((1 - ax) * (1 - ay) * F[kk, i, j] + ax * (1 - ay) * F[kk, i + 1, j]
                    + (1 - ax) * ay * F[kk, i, j + 1] + ax * ay * F[kk, i + 1, j + 1])
        if nl == 1:
            return bil(self.u, k), bil(self.v, k)
        return ((1 - az) * bil(self.u, k) + az * bil(self.u, k + 1),
                (1 - az) * bil(self.v, k) + az * bil(self.v, k + 1))


def solve_urban_wind(grid: dict, z_levels, p: WindParams, layer_winds=None) -> UrbanWind:
    """grid: dict with H, x0, y0, res (from GeoCorridor.grid).

    ``layer_winds`` optionally gives the background (speed m/s, toward-deg relative to the
    corridor) per level, e.g. from :func:`aamsim.weather.layer_winds`; otherwise the
    power-law profile of ``p`` is used."""
    H = grid["H"]
    zs = np.array([max(float(z), 10.0) for z in z_levels])
    us, vs = [], []
    for k, z in enumerate(zs):
        if layer_winds is not None:
            s, d = layer_winds[k]
        else:
            s, d = layer_speed(p, z), p.direction_deg
        th = math.radians(d)
        ux, uy = s * math.cos(th), s * math.sin(th)
        u, v = solve_layer(H >= z, grid["res"], grid["x0"], grid["y0"], ux, uy, p)
        us.append(u); vs.append(v)
    return UrbanWind(zs, np.array(us), np.array(vs), grid["x0"], grid["y0"], grid["res"])


def wind_at_nodes(uw: UrbanWind, corr):
    """WindField (mph, one time bin) at the corridor nodes."""
    from .wind import WindField
    u, v = uw.sample(corr.pos[:, 0], corr.pos[:, 1], corr.pos[:, 2])
    return WindField(np.atleast_2d(u * MPH_PER_MS), np.atleast_2d(v * MPH_PER_MS), 3600.0, "urban")


def geo_beliefs(uw: UrbanWind, corr, p: WindParams, dir_error_deg=15.0, speed_bias=-0.15):
    """The four planning winds at corridor nodes (mph).

    ground_truth : the building-aware field
    forecast     : coarse weather-model wind: free-stream shear profile, no buildings,
                   with a direction error and a speed bias
    profile      : true field averaged over each altitude layer (height-only)
    wind_agnostic: zero
    """
    from .wind import WindField
    truth = wind_at_nodes(uw, corr)
    z = corr.pos[:, 2]
    th = math.radians(p.direction_deg + dir_error_deg)
    s = np.array([layer_speed(p, zz) for zz in z]) * (1 + speed_bias) * MPH_PER_MS
    fc = WindField(np.atleast_2d(s * math.cos(th)), np.atleast_2d(s * math.sin(th)), 3600.0, "forecast")
    pu, pv = truth.u.copy(), truth.v.copy()
    for lev in np.unique(corr.ijk[:, 2]):
        m = corr.ijk[:, 2] == lev
        pu[:, m] = truth.u[:, m].mean(); pv[:, m] = truth.v[:, m].mean()
    prof = WindField(pu, pv, 3600.0, "profile")
    zero = WindField(np.zeros_like(truth.u), np.zeros_like(truth.v), 3600.0, "wind_agnostic")
    return {"ground_truth": truth.renamed("ground_truth"), "forecast": fc, "profile": prof,
            "wind_agnostic": zero}
