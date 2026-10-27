"""Real-city geometry: OpenStreetMap buildings and a corridor laid over them.

Workflow
    1. ``fetch_osm(bbox, path)``        download an .osm extract (OSM API, small areas)
    2. ``load_buildings(path, origin)`` footprints in local metres + heights
    3. ``build_geo_corridor(...)``      the modular corridor of the papers, aligned with a
                                        start->goal axis, with every node and link that
                                        intersects a building (plus clearance) removed

Local frame: metres, x along the corridor axis (start -> goal), y to its left,
z up. ``to_local`` / ``to_corridor`` convert between lat/lon, an east-north frame
and the corridor frame.

Building heights: ``height`` tag (m or ft), else ``building:levels`` x 3.5 m,
else ``min_height``-aware default by building type. Building parts are used when
present because tall towers in OSM are usually modelled as parts.
"""
from __future__ import annotations

import math
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .corridor import Corridor

M_PER_DEG_LAT = 111132.0
LEVEL_M = 3.5
DEFAULT_HEIGHT = {"house": 7.0, "residential": 9.0, "garage": 4.0, "garages": 4.0, "shed": 3.0,
                  "roof": 5.0, "parking": 12.0, "commercial": 12.0, "retail": 7.0,
                  "industrial": 9.0, "warehouse": 9.0, "office": 20.0, "apartments": 15.0,
                  "church": 15.0, "school": 10.0, "hotel": 30.0, "yes": 10.0}


def fetch_osm(bbox, path):
    """Download raw OSM XML for bbox = (min_lon, min_lat, max_lon, max_lat)."""
    url = ("https://api.openstreetmap.org/api/0.6/map?bbox="
           + ",".join(f"{v:.6f}" for v in bbox))
    req = urllib.request.Request(url, headers={"User-Agent": "aam-corridor-lab/0.1"})
    with urllib.request.urlopen(req, timeout=180) as r:
        raw = r.read()
    if str(path).endswith(".gz"):
        import gzip
        raw = gzip.compress(raw)
    Path(path).write_bytes(raw)
    return path


def _parse_height(tags):
    def num(s):
        s = s.strip().lower().replace(",", ".")
        try:
            if s.endswith("ft") or s.endswith("'"):
                return float(s.rstrip("ft'").strip()) * 0.3048
            return float(s.split()[0].rstrip("m"))
        except (ValueError, IndexError):
            return None
    h = num(tags["height"]) if "height" in tags else None
    if h is None and "building:levels" in tags:
        lv = num(tags["building:levels"])
        h = lv * LEVEL_M if lv else None
    if h is None:
        h = DEFAULT_HEIGHT.get(tags.get("building", tags.get("building:part", "yes")), 10.0)
    base = num(tags["min_height"]) if "min_height" in tags else None
    if base is None and "building:min_level" in tags:
        ml = num(tags["building:min_level"])
        base = ml * LEVEL_M if ml else 0.0
    return float(h), float(base or 0.0), ("height" in tags or "building:levels" in tags)


@dataclass
class Building:
    xy: np.ndarray          # (k, 2) footprint in east-north metres (closed ring not repeated)
    height: float
    base: float = 0.0
    measured: bool = False  # height came from a tag, not a default
    osm_id: str = ""


@dataclass
class City:
    name: str
    origin: tuple            # (lat, lon) of the local frame origin
    buildings: list = field(default_factory=list)

    def to_local(self, lat, lon):
        lat0, lon0 = self.origin
        return np.array([(lon - lon0) * M_PER_DEG_LAT * math.cos(math.radians(lat0)),
                         (lat - lat0) * M_PER_DEG_LAT])


def _xml_root(path):
    import gzip
    path = str(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return ET.parse(f).getroot()
    return ET.parse(path).getroot()


def load_buildings(path, origin=None, name="city") -> City:
    root = _xml_root(path)
    nodes = {n.get("id"): (float(n.get("lat")), float(n.get("lon"))) for n in root.iter("node")}
    if origin is None:
        b = root.find("bounds")
        origin = ((float(b.get("minlat")) + float(b.get("maxlat"))) / 2,
                  (float(b.get("minlon")) + float(b.get("maxlon"))) / 2)
    city = City(name, origin)
    ways = {}
    for w in root.iter("way"):
        ways[w.get("id")] = w
    parts, whole = [], []
    used_as_outline = set()
    for rel in root.iter("relation"):  # building outlines that have parts: skip the outline
        tags = {t.get("k"): t.get("v") for t in rel.iter("tag")}
        if tags.get("type") == "building":
            for m in rel.iter("member"):
                if m.get("role") == "outline":
                    used_as_outline.add(m.get("ref"))
    for wid, w in ways.items():
        tags = {t.get("k"): t.get("v") for t in w.iter("tag")}
        refs = [nd.get("ref") for nd in w.iter("nd")]
        if len(refs) < 4 or refs[0] != refs[-1]:
            continue
        if "building:part" in tags:
            parts.append((wid, tags, refs))
        elif "building" in tags and wid not in used_as_outline:
            whole.append((wid, tags, refs))
    # outlines fully covered by parts are dropped (towers modelled with parts)
    part_polys = []
    for wid, tags, refs in parts + whole:
        try:
            ll = np.array([nodes[r] for r in refs[:-1]])
        except KeyError:
            continue
        xy = np.stack([(ll[:, 1] - origin[1]) * M_PER_DEG_LAT * math.cos(math.radians(origin[0])),
                       (ll[:, 0] - origin[0]) * M_PER_DEG_LAT], axis=1)
        h, base, measured = _parse_height(tags)
        if h - base < 2.0:
            continue
        b = Building(xy, h, base, measured, wid)
        (part_polys if (wid, tags, refs) in parts else city.buildings).append(b)
    if part_polys:
        # drop an outline if a part's centroid lies inside it and the outline has no measured height
        keep = []
        for b in city.buildings:
            inside = any(point_in_poly(p.xy.mean(axis=0), b.xy) for p in part_polys)
            if not (inside and not b.measured):
                keep.append(b)
        city.buildings = keep + part_polys
    return city


def point_in_poly(p, poly):
    x, y = p
    inside = False
    n = len(poly)
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi:
            inside = not inside
        j = i
    return inside


# ---------------------------------------------------------------------------- corridor

@dataclass
class GeoFrame:
    """Rigid transform between east-north metres and the corridor frame."""
    start_en: np.ndarray
    bearing: float          # radians, angle of the corridor axis measured from east (CCW)

    def to_corridor(self, en):
        d = np.asarray(en, float) - self.start_en
        c, s = math.cos(self.bearing), math.sin(self.bearing)
        return np.stack([d[..., 0] * c + d[..., 1] * s, -d[..., 0] * s + d[..., 1] * c], axis=-1)

    def to_en(self, xy):
        xy = np.asarray(xy, float)
        c, s = math.cos(self.bearing), math.sin(self.bearing)
        return np.stack([xy[..., 0] * c - xy[..., 1] * s, xy[..., 0] * s + xy[..., 1] * c],
                        axis=-1) + self.start_en


def rasterize(city: City, frame: GeoFrame, x0, y0, nx, ny, res):
    """Building-top height (m) on a corridor-frame grid; cell (i, j) centre at
    (x0 + (i + .5) res, y0 + (j + .5) res)."""
    H = np.zeros((nx, ny))
    gx = x0 + (np.arange(nx) + 0.5) * res
    gy = y0 + (np.arange(ny) + 0.5) * res
    for b in city.buildings:
        poly = frame.to_corridor(b.xy)
        lo, hi = poly.min(0), poly.max(0)
        i0, i1 = np.searchsorted(gx, lo[0]), np.searchsorted(gx, hi[0])
        j0, j1 = np.searchsorted(gy, lo[1]), np.searchsorted(gy, hi[1])
        if i1 <= i0 or j1 <= j0:
            # tiny footprint: mark nearest cell
            i = int(np.clip((poly[:, 0].mean() - x0) // res, 0, nx - 1))
            j = int(np.clip((poly[:, 1].mean() - y0) // res, 0, ny - 1))
            H[i, j] = max(H[i, j], b.height)
            continue
        X, Y = np.meshgrid(gx[i0:i1], gy[j0:j1], indexing="ij")
        inside = _pip_grid(X.ravel(), Y.ravel(), poly).reshape(X.shape)
        sub = H[i0:i1, j0:j1]
        sub[inside] = np.maximum(sub[inside], b.height)
    return np.round(H, 1)


def load_roads(path, origin):
    """Street centre lines (highway=*) as lists of east-north points, for map context."""
    root = _xml_root(path)
    nodes = {n.get("id"): (float(n.get("lat")), float(n.get("lon"))) for n in root.iter("node")}
    keep = {"motorway", "trunk", "primary", "secondary", "tertiary", "residential", "unclassified",
            "motorway_link", "trunk_link", "primary_link", "secondary_link", "service"}
    roads = []
    c = math.cos(math.radians(origin[0]))
    for w in root.iter("way"):
        tags = {t.get("k"): t.get("v") for t in w.iter("tag")}
        if tags.get("highway") not in keep:
            continue
        pts = [nodes[r.get("ref")] for r in w.iter("nd") if r.get("ref") in nodes]
        if len(pts) < 2:
            continue
        ll = np.array(pts)
        en = np.stack([(ll[:, 1] - origin[1]) * M_PER_DEG_LAT * c, (ll[:, 0] - origin[0]) * M_PER_DEG_LAT], 1)
        major = tags["highway"] in ("motorway", "trunk", "primary", "motorway_link", "trunk_link")
        roads.append((en, major))
    return roads


def _pip_grid(px, py, poly):
    inside = np.zeros(px.shape, bool)
    n = len(poly)
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        cond = (yi > py) != (yj > py)
        xint = (xj - xi) * (py - yi) / (yj - yi + 1e-12) + xi
        inside ^= cond & (px < xint)
        j = i
    return inside


@dataclass
class GeoCorridor:
    city: City
    frame: GeoFrame
    corridor: Corridor
    grid: dict               # x0, y0, nx, ny, res, H (building heights on the corridor grid)
    removed_nodes: int = 0
    removed_edges: int = 0


def build_geo_corridor(city: City, start_latlon, goal_latlon, module_length=150.0,
                       module_width=50.0, n_lat=6, z_levels=(0, 30, 60, 90, 120),
                       clearance=8.0, res=10.0, margin=200.0, portal_depth=2) -> GeoCorridor:
    """Lay the modular corridor along start->goal and drop everything inside buildings.

    Defaults: 150 m x 50 m modules, 7 lanes (350 m wide), layers every 30 m up to the
    120 m (400 ft) Part 107 ceiling. The start and goal are the centre lanes' ground
    nodes; ground portals inside a building footprint are removed.
    """
    from .corridor import CorridorSpec, build_corridor
    a = city.to_local(*start_latlon)
    b = city.to_local(*goal_latlon)
    L = float(np.hypot(*(b - a)))
    frame = GeoFrame(a, math.atan2(b[1] - a[1], b[0] - a[0]))
    n_long = max(1, int(round(L / module_length)))
    spec = CorridorSpec(n_long=n_long, n_lat=n_lat, n_alt=len(z_levels) - 1,
                        module_length=L / n_long, module_width=module_width,
                        z_levels=tuple(float(z) for z in z_levels), vertical_ground_only=True)
    full = build_corridor(spec)
    # corridor frame: centre the lanes on the axis
    half = n_lat * module_width / 2
    pos = full.pos.copy()
    pos[:, 1] -= half
    x0, y0 = -margin, -half - margin
    nx_, ny_ = int((L + 2 * margin) / res) + 1, int((2 * half + 2 * margin) / res) + 1
    H = rasterize(city, frame, x0, y0, nx_, ny_, res)

    def height_at(x, y):
        i = np.clip(((x - x0) / res).astype(int), 0, nx_ - 1)
        j = np.clip(((y - y0) / res).astype(int), 0, ny_ - 1)
        return H[i, j]

    # obstacle test with lateral clearance: max building height within +-clearance
    pad = int(math.ceil(clearance / res))
    Hpad = H.copy()
    for di in range(-pad, pad + 1):
        for dj in range(-pad, pad + 1):
            Hpad = np.maximum(Hpad, np.roll(np.roll(H, di, 0), dj, 1))

    def blocked(x, y, z):
        i = np.clip(((x - x0) / res).astype(int), 0, nx_ - 1)
        j = np.clip(((y - y0) / res).astype(int), 0, ny_ - 1)
        return z <= Hpad[i, j] + clearance

    node_ok = ~blocked(pos[:, 0], pos[:, 1], np.maximum(pos[:, 2], 0.5))
    ground = full.ijk[:, 2] == 0
    gi = np.clip(((pos[:, 0] - x0) / res).astype(int), 0, nx_ - 1)
    gj = np.clip(((pos[:, 1] - y0) / res).astype(int), 0, ny_ - 1)
    node_ok = np.where(ground, Hpad[gi, gj] <= 0.0, node_ok)  # ground pads: no building nearby
    edge_ok = node_ok[full.edge_nodes[:, 0]] & node_ok[full.edge_nodes[:, 1]]
    # sample along each edge
    P0, P1 = pos[full.edge_nodes[:, 0]], pos[full.edge_nodes[:, 1]]
    for f in np.linspace(0.1, 0.9, 9):
        P = P0 + (P1 - P0) * f
        vertical = np.all(np.abs(P1[:, :2] - P0[:, :2]) < 1e-6, axis=1)
        hit = blocked(P[:, 0], P[:, 1], np.maximum(P[:, 2], 0.5)) & ~vertical
        # vertical take-off / landing columns must be free of buildings over their full height
        hit |= vertical & (Hpad[np.clip(((P[:, 0] - x0) / res).astype(int), 0, nx_ - 1),
                                np.clip(((P[:, 1] - y0) / res).astype(int), 0, ny_ - 1)] > 0)
        edge_ok &= ~hit
    keep_edges = full.edge_nodes[edge_ok]
    # reindex nodes that still touch an edge
    used = np.zeros(full.n_nodes, bool)
    used[keep_edges.ravel()] = True
    new_id = -np.ones(full.n_nodes, int)
    new_id[used] = np.arange(used.sum())
    en = keep_edges
    e2 = new_id[en]
    a_from = np.concatenate([e2[:, 0], e2[:, 1]])
    a_to = np.concatenate([e2[:, 1], e2[:, 0]])
    p2 = pos[used]
    corr = Corridor(spec=spec, pos=p2, ijk=full.ijk[used], arc_from=a_from, arc_to=a_to,
                    arc_len=np.linalg.norm(p2[a_to] - p2[a_from], axis=1),
                    arc_edge=np.concatenate([np.arange(len(e2))] * 2), edge_nodes=e2,
                    origins=[int(new_id[n]) for n in range(full.n_nodes)
                             if used[n] and ground[n] and full.ijk[n, 0] < portal_depth],
                    destinations=[int(new_id[n]) for n in range(full.n_nodes)
                                  if used[n] and ground[n] and full.ijk[n, 0] > n_long - portal_depth])
    return GeoCorridor(city, frame, corr,
                       dict(x0=x0, y0=y0, nx=nx_, ny=ny_, res=res, H=H, half=half, L=L),
                       removed_nodes=int(full.n_nodes - used.sum()),
                       removed_edges=int(full.n_edges - len(e2)))


def nearest_portal(gc: GeoCorridor, which: str, lane: int | None = None):
    """Ground node on the start ('origin') or goal face closest to the centre lane."""
    c = gc.corridor
    ids = c.origins if which == "origin" else c.destinations
    if not ids:
        raise ValueError(f"no open ground portal on the {which} face")
    end = 0.0 if which == "origin" else c.pos[:, 0].max()
    return min(ids, key=lambda n: (abs(c.pos[n, 0] - end), abs(c.pos[n, 1])) if lane is None
               else abs(c.ijk[n, 1] - lane))


DALLAS = dict(name="Downtown Dallas", bbox=(-96.8148, 32.7723, -96.7851, 32.7950),
              start=(32.7768, -96.7896), goal=(32.7905, -96.8103),
              start_name="Dallas Farmers Market", goal_name="American Airlines Center")
