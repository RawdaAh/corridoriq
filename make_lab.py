"""Build demo/lab.html - the AAM Corridor Lab (all three papers on a real OSM city).

    python make_lab.py                 # Dallas, weather snapshot from NWS/Open-Meteo at build time
    python make_lab.py --osm data/osm/mycity.osm.gz --start LAT LON --goal LAT LON --name "My city"

Needs demo/rh_cases.json (python export_rh_cases.py) for the rolling-horizon tab.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from aamsim.geo import DALLAS, load_buildings, load_roads
from aamsim.urbanwind import WindParams

# Published values used for validation columns ------------------------------------------
REF = dict(
    trr=dict(  # Table 1 / Table 3: [time min, energy kWh, paths, ascents, descents]
        t1={"none": [9.74, 0.170, 1, 1, 1], "0": [5.50, 0.054, 1, 3, 3], "45": [6.30, 0.070, 3, 3, 3],
            "90": [9.74, 0.170, 1, 1, 1], "135": [9.88, 0.174, 1, 1, 1], "180": [16.29, 0.477, 1, 1, 1],
            "225": [10.04, 0.175, 1, 1, 1], "270": [9.74, 0.170, 1, 1, 1], "315": [6.38, 0.073, 1, 3, 3]},
        t3={"none": [9.74, 0.170, 1, 1, 1], "0": [5.95, 0.063, 2, 3, 3], "45": [6.71, 0.079, 4, 3, 3],
            "90": [9.74, 0.170, 1, 1, 1], "135": [9.88, 0.174, 1, 1, 1], "180": [17.51, 0.559, 1, 1, 1],
            "225": [10.08, 0.176, 2, 2, 2], "270": [9.74, 0.170, 1, 1, 1], "315": [6.78, 0.082, 2, 3, 3]},
        f8=[[0, 100], [10, 100], [20, 100], [30, 100], [40, 100], [50, 100], [60, 100], [70, 100], [80, 100], [90, 40]]),
    trb=dict(sr={"200": 0.9957, "400": 0.9684, "600": 0.8968, "800": 0.8066, "1000": 0.6997,
                 "2000": 0.3686, "4000": 0.1768},
             dir={"1": [36.86, 39.02], "0.75": [42.76, 62.28], "0.5": [44.20, 72.17], "0.25": [42.77, 62.37]}))


def weather_snapshot(lat, lon):
    from aamsim.weather import fetch_profile
    try:
        p = fetch_profile(lat, lon)
        p["source"] += " · snapshot at build time"
        return p
    except Exception as e:  # noqa: BLE001
        print("weather fetch failed:", e)
        return dict(source="default (no network at build time)", time="", lat=lat, lon=lon,
                    heights_m=[10, 80, 120, 180], speed_ms=[3.0, 4.5, 5.0, 5.5], dir_from_deg=[180, 185, 190, 190])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--osm", default="data/osm/dallas.osm.gz")
    ap.add_argument("--name", default=DALLAS["name"])
    ap.add_argument("--start", type=float, nargs=2, default=DALLAS["start"])
    ap.add_argument("--goal", type=float, nargs=2, default=DALLAS["goal"])
    ap.add_argument("--out", default="demo/lab.html")
    a = ap.parse_args()
    city = load_buildings(a.osm, name=a.name)
    roads = load_roads(a.osm, city.origin)
    r1 = lambda v: [round(float(x), 1) for x in np.asarray(v).ravel()]  # noqa: E731
    data = dict(
        city=dict(name=a.name, origin=list(city.origin),
                  buildings=[dict(p=r1(b.xy), h=round(b.height, 1), m=int(b.measured)) for b in city.buildings],
                  roads=[dict(p=r1(en), M=int(m)) for en, m in roads],
                  start_en=r1(city.to_local(*a.start)), goal_en=r1(city.to_local(*a.goal))),
        windp={k: v for k, v in WindParams().__dict__.items()},
        weather=weather_snapshot(*city.origin), ref=REF,
        rh=json.loads(Path("demo/rh_cases.json").read_text()))
    eng = Path("demo/lab_engine.js").read_text()
    html = (Path("demo/lab_template.html").read_text()
            .replace("/*__ENGINE__*/", eng)
            .replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")).replace("</", "<\\/"))
            .replace('/*__ENGINE_SRC__*/""', json.dumps(eng).replace("</", "<\\/")))
    Path(a.out).write_text(html)
    print(f"{a.out}: {len(html) / 1e6:.2f} MB · {len(data['city']['buildings'])} buildings · weather: "
          f"{data['weather']['source']} {data['weather']['speed_ms']}")


if __name__ == "__main__":
    main()
