"""Everything the app draws or lets people pick, from one OpenStreetMap extract.

    scene = load_scene("data/osm/dallas.osm.gz")
    scene["buildings"]  footprints (east-north metres), height, kind (office, residential, ...)
    scene["roads"]      centre lines with a carriageway width by road class
    scene["areas"]      parks, grass, water, plazas, parking, rail yards
    scene["trees"]      individual trees (natural=tree)
    scene["places"]     named places people can pick as origin or destination
"""
from __future__ import annotations

import math

import numpy as np

from .geo import M_PER_DEG_LAT, _xml_root, load_buildings

ROAD_WIDTH = {"motorway": 22, "trunk": 18, "primary": 15, "secondary": 13, "tertiary": 11,
              "residential": 9, "unclassified": 8, "service": 5, "motorway_link": 9, "trunk_link": 8,
              "primary_link": 8, "secondary_link": 7, "living_street": 6, "pedestrian": 6, "footway": 2.5}
AREA = [("water", lambda t: t.get("natural") == "water" or t.get("waterway") == "riverbank"
         or t.get("water") is not None),
        ("park", lambda t: t.get("leisure") in ("park", "garden", "dog_park", "playground")),
        ("grass", lambda t: t.get("landuse") in ("grass", "meadow", "recreation_ground", "village_green")
         or t.get("leisure") == "pitch" or t.get("natural") in ("scrub", "wood")),
        ("plaza", lambda t: t.get("highway") == "pedestrian" and t.get("area") == "yes"
         or t.get("place") == "square"),
        ("parking", lambda t: t.get("amenity") == "parking"),
        ("rail", lambda t: t.get("landuse") == "railway")]
PLACE_KEYS = ("building", "amenity", "tourism", "leisure", "public_transport", "shop", "office")
SKIP_PLACE = {"bench", "waste_basket", "bicycle_parking", "parking_entrance", "parking_space",
              "drinking_water", "toilets", "fountain", "shelter", "loading_dock", "picnic_table",
              "artwork", "vending_machine", "atm", "charging_station", "post_box"}


def load_scene(path, name="city"):
    city = load_buildings(path, name=name)
    root = _xml_root(path)
    lat0, lon0 = city.origin
    kx, ky = M_PER_DEG_LAT * math.cos(math.radians(lat0)), M_PER_DEG_LAT
    nodes = {}
    trees, places = [], []
    for n in root.iter("node"):
        ll = (float(n.get("lat")), float(n.get("lon")))
        nodes[n.get("id")] = ll
        tags = {t.get("k"): t.get("v") for t in n.iter("tag")}
        if tags.get("natural") == "tree":
            trees.append(((ll[1] - lon0) * kx, (ll[0] - lat0) * ky))
        elif "name" in tags and any(k in tags for k in PLACE_KEYS) and tags.get("amenity") not in SKIP_PLACE \
                and not _is_stop(tags):
            places.append(dict(name=tags["name"], kind=_kind(tags),
                               en=[(ll[1] - lon0) * kx, (ll[0] - lat0) * ky]))
    en = lambda ll: ((ll[1] - lon0) * kx, (ll[0] - lat0) * ky)  # noqa: E731
    roads, areas = [], []
    for w in root.iter("way"):
        tags = {t.get("k"): t.get("v") for t in w.iter("tag")}
        refs = [r.get("ref") for r in w.iter("nd") if r.get("ref") in nodes]
        if len(refs) < 2:
            continue
        pts = np.array([en(nodes[r]) for r in refs])
        closed = refs[0] == refs[-1] and len(refs) >= 4
        hw = tags.get("highway")
        if hw in ROAD_WIDTH and not (closed and tags.get("area") == "yes"):
            roads.append(dict(p=pts, w=ROAD_WIDTH[hw], major=hw in ("motorway", "trunk", "primary",
                                                                    "motorway_link", "trunk_link")))
        if closed:
            for kind, test in AREA:
                if test(tags):
                    areas.append(dict(p=pts[:-1], kind=kind))
                    break
            if "name" in tags and any(k in tags for k in PLACE_KEYS) and tags.get("amenity") not in SKIP_PLACE:
                places.append(dict(name=tags["name"], kind=_kind(tags), en=list(pts[:-1].mean(0))))
    # de-duplicate places by name (keep the first, usually the building)
    seen, uniq = set(), []
    for p in places:
        k = p["name"].strip().lower()
        if k in seen or len(k) < 3:
            continue
        seen.add(k)
        uniq.append(p)
    return dict(city=city, roads=roads, areas=areas, trees=trees, places=uniq)


def _is_stop(tags):
    return (tags.get("highway") == "bus_stop" or tags.get("public_transport") in ("platform", "stop_position")) \
        and tags.get("railway") not in ("station", "halt")


def _kind(tags):
    for key, label in (("tourism", "attraction"), ("leisure", "park"), ("public_transport", "transit"),
                       ("railway", "transit"), ("amenity", None), ("office", "office"), ("shop", "shop")):
        if key in tags:
            v = tags[key]
            if key == "amenity":
                return {"hospital": "hospital", "parking": "parking", "school": "school",
                        "university": "school", "restaurant": "food", "cafe": "food",
                        "fast_food": "food", "bar": "food"}.get(v, "amenity")
            return label
    b = tags.get("building", "")
    return {"hotel": "hotel", "office": "office", "commercial": "office", "apartments": "residential",
            "residential": "residential", "church": "landmark", "civic": "civic",
            "government": "civic", "stadium": "landmark", "parking": "parking"}.get(b, "building")
