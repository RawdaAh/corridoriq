"""AAM Corridor Lab: one simulator for wind-aware routing (TRR), shared-intent
reservation (TRB/SIRE) and rolling-horizon corridor scheduling."""
from .corridor import Corridor, CorridorSpec, build_corridor, trr_wind_corridor, trb_sire_corridor
from .energy import VehicleModel

__version__ = "0.1.0"
