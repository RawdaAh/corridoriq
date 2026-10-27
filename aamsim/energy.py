"""Link travel-time and energy models (TRR paper Eqs. 10-13 and 26-33)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .corridor import Corridor, mph_to_fps

J_PER_KWH = 3.6e6


@dataclass
class VehicleModel:
    """Small multirotor parameters used across the three papers."""

    cruise_mph: float = 35.0          # v_UAV
    energy_j_per_ft: float = 20.42    # Rodrigues et al.: 67 J/m at 0.5 kg payload
    phase_dependent: bool = False     # TRR Sec. 5.4 climb / take-off terms
    k_climb: float = 0.2              # alpha_climb = 1 + k * sin(ascent angle)
    alpha_desc: float = 1.0
    takeoff_fraction: float = 0.10    # E_TO as a fraction of the no-wind reference energy
    takeoff_reference_ft: float = 30000.0
    battery_kwh: float = np.inf       # energy feasibility (state of charge) limit
    length_unit_ft: float = 1.0       # corridor length unit in feet (3.28084 for metre geometry)

    @property
    def cruise_fps(self) -> float:
        return mph_to_fps(self.cruise_mph)

    def takeoff_energy_j(self) -> float:
        if not self.phase_dependent:
            return 0.0
        return self.takeoff_fraction * self.energy_j_per_ft * self.takeoff_reference_ft


def geo_vehicle(**kw) -> VehicleModel:
    """Vehicle for corridors built in metres (real-city geometry)."""
    return VehicleModel(length_unit_ft=3.28084, **kw)


def baseline_link_costs(corr: Corridor, veh: VehicleModel):
    """No-wind travel time (s) and energy (J) of every arc (Eqs. 11-12)."""
    L = corr.arc_len * veh.length_unit_ft
    t = L / veh.cruise_fps
    e = veh.energy_j_per_ft * L
    if veh.phase_dependent:
        dz = corr.pos[corr.arc_to, 2] - corr.pos[corr.arc_from, 2]
        sin_phi = np.where(corr.arc_len > 0, dz / corr.arc_len, 0.0)        # Eq. 31
        mult = np.where(dz > 0, 1.0 + veh.k_climb * sin_phi,                 # Eqs. 29, 32
                        np.where(dz < 0, veh.alpha_desc, 1.0))               # Eqs. 28, 30
        e = e * mult
    return t, e


def wind_link_costs(corr: Corridor, veh: VehicleModel, wind, model: str = "paper"):
    """Time-dependent link costs (Eqs. 26-27). Returns arrays of shape (B, A)."""
    t0, e0 = baseline_link_costs(corr, veh)
    lam, beta = wind.link_factors(corr, veh.cruise_mph, model)
    return t0[None, :] * lam, e0[None, :] * beta
