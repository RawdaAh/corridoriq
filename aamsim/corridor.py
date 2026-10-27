"""Modular 3D corridor lattice shared by all three studies.

The corridor is built from identical rectangular prism modules. Within each
module every pair of the 8 corners is connected (12 prism edges, 12 face
diagonals, 4 body diagonals). Links shared by neighbouring modules are merged.
Ground-to-ground links are always removed (no taxiing). Optionally, ground
access can be restricted to purely vertical take-off / landing links.

Coordinates are in feet: x longitudinal, y lateral, z altitude above ground.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

import numpy as np

FT_PER_MILE = 5280.0
SEC_PER_HOUR = 3600.0


def mph_to_fps(mph: float) -> float:
    return mph * FT_PER_MILE / SEC_PER_HOUR


@dataclass
class CorridorSpec:
    """Geometry of a modular corridor.

    n_long, n_lat, n_alt : number of modules along x, y, z.
    module_length, module_width : module size along x and y (ft).
    z_levels : explicit altitude of each node layer (ft). If None, layers are
        evenly spaced by ``module_height`` starting at ground (0 ft).
    vertical_ground_only : if True, a link touching the ground must be purely
        vertical (take-off / landing straight up). If False, diagonal
        ground-air links are kept (the setting used in the TRR and TRB papers).
    """

    n_long: int = 30
    n_lat: int = 2
    n_alt: int = 4
    module_length: float = 1000.0
    module_width: float = 200.0
    module_height: float = 100.0
    z_levels: tuple | None = None
    vertical_ground_only: bool = False

    def altitudes(self) -> np.ndarray:
        if self.z_levels is not None:
            z = np.asarray(self.z_levels, dtype=float)
            if len(z) != self.n_alt + 1:
                raise ValueError("z_levels must have n_alt + 1 entries")
            if z[0] != 0 or np.any(np.diff(z) <= 0):
                raise ValueError("z_levels must start at 0 and increase")
            return z
        return np.arange(self.n_alt + 1) * self.module_height


@dataclass
class Corridor:
    """Node/arc arrays for a corridor. Built with :func:`build_corridor`."""

    spec: CorridorSpec
    pos: np.ndarray            # (N, 3) node coordinates in ft
    ijk: np.ndarray            # (N, 3) integer module indices
    arc_from: np.ndarray       # (A,) directed arcs
    arc_to: np.ndarray
    arc_len: np.ndarray        # (A,) 3D length in ft
    arc_edge: np.ndarray       # (A,) undirected edge id of each arc
    edge_nodes: np.ndarray     # (E, 2) node pair of each undirected edge (u < v)
    origins: list = field(default_factory=list)       # ground nodes on entry face
    destinations: list = field(default_factory=list)  # ground nodes on exit face

    @property
    def n_nodes(self) -> int:
        return len(self.pos)

    @property
    def n_arcs(self) -> int:
        return len(self.arc_from)

    @property
    def n_edges(self) -> int:
        return len(self.edge_nodes)

    def node(self, i: int, j: int, k: int) -> int:
        s = self.spec
        return (i * (s.n_lat + 1) + j) * (s.n_alt + 1) + k

    @property
    def ground(self) -> np.ndarray:
        return self.ijk[:, 2] == 0

    def out_arcs(self) -> list:
        """Adjacency list: out_arcs()[u] -> array of arc ids leaving u."""
        if not hasattr(self, "_out"):
            order = np.argsort(self.arc_from, kind="stable")
            counts = np.bincount(self.arc_from, minlength=self.n_nodes)
            splits = np.cumsum(counts)[:-1]
            self._out = np.split(order, splits)
        return self._out

    def arc_index(self) -> dict:
        if not hasattr(self, "_aidx"):
            self._aidx = {(int(u), int(v)): a for a, (u, v) in
                          enumerate(zip(self.arc_from, self.arc_to))}
        return self._aidx

    def heading_deg(self) -> np.ndarray:
        """Horizontal heading of each arc, degrees counter-clockwise from +x."""
        d = self.pos[self.arc_to] - self.pos[self.arc_from]
        return np.degrees(np.arctan2(d[:, 1], d[:, 0])) % 360.0

    def summary(self) -> dict:
        return dict(nodes=self.n_nodes, arcs=self.n_arcs, edges=self.n_edges,
                    origins=len(self.origins), destinations=len(self.destinations),
                    length_ft=float(self.pos[:, 0].max()),
                    width_ft=float(self.pos[:, 1].max()),
                    ceiling_ft=float(self.pos[:, 2].max()))


def build_corridor(spec: CorridorSpec | None = None, **kw) -> Corridor:
    """Generate the modular lattice described in the three papers."""
    spec = spec or CorridorSpec(**kw)
    z = spec.altitudes()
    NL, NW, NH = spec.n_long, spec.n_lat, spec.n_alt
    ii, jj, kk = np.meshgrid(np.arange(NL + 1), np.arange(NW + 1), np.arange(NH + 1),
                             indexing="ij")
    ijk = np.stack([ii.ravel(), jj.ravel(), kk.ravel()], axis=1)
    pos = np.stack([ijk[:, 0] * spec.module_length,
                    ijk[:, 1] * spec.module_width,
                    z[ijk[:, 2]]], axis=1).astype(float)

    def nid(i, j, k):
        return (i * (NW + 1) + j) * (NH + 1) + k

    edges = set()
    for i in range(NL):
        for j in range(NW):
            for k in range(NH):
                corners = [nid(i + a, j + b, k + c)
                           for a in (0, 1) for b in (0, 1) for c in (0, 1)]
                for u, v in combinations(corners, 2):
                    if pos[u, 2] == 0 and pos[v, 2] == 0:
                        continue  # no ground-to-ground movement
                    if spec.vertical_ground_only and (pos[u, 2] == 0 or pos[v, 2] == 0):
                        if not np.allclose(pos[u, :2], pos[v, :2]):
                            continue
                    edges.add((min(u, v), max(u, v)))
    edge_nodes = np.array(sorted(edges), dtype=np.int64)
    a_from = np.concatenate([edge_nodes[:, 0], edge_nodes[:, 1]])
    a_to = np.concatenate([edge_nodes[:, 1], edge_nodes[:, 0]])
    a_edge = np.concatenate([np.arange(len(edge_nodes))] * 2)
    a_len = np.linalg.norm(pos[a_to] - pos[a_from], axis=1)

    origins = [nid(0, j, 0) for j in range(NW + 1)]
    destinations = [nid(NL, j, 0) for j in range(NW + 1)]
    return Corridor(spec=spec, pos=pos, ijk=ijk, arc_from=a_from, arc_to=a_to,
                    arc_len=a_len, arc_edge=a_edge, edge_nodes=edge_nodes,
                    origins=origins, destinations=destinations)


# Named geometries used in the papers -------------------------------------------------

def trr_wind_corridor(n_lat: int = 2) -> Corridor:
    """Paper 1 (TRR, wind-aware routing): 5 x 6,000 ft blocks, layers 0/200/300/400 ft,
    100 ft lateral lane spacing, diagonal take-off allowed."""
    return build_corridor(CorridorSpec(n_long=5, n_lat=n_lat, n_alt=3,
                                       module_length=6000.0, module_width=100.0,
                                       z_levels=(0.0, 200.0, 300.0, 400.0)))


def trb_sire_corridor(n_long: int = 30, n_lat: int = 2, n_alt: int = 4) -> Corridor:
    """Paper 2 (TRB, SIRE): 30 x 2 x 4 modules of 1,000 x 200 x 100 ft."""
    return build_corridor(CorridorSpec(n_long=n_long, n_lat=n_lat, n_alt=n_alt,
                                       module_length=1000.0, module_width=200.0,
                                       module_height=100.0))
