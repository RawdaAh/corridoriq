"""Build CorridorIQ (demo/corridoriq.html): the product app on top of the research engine.

    python make_app.py                                   # Dallas, live forecast at build time
    python make_app.py --osm data/osm/mycity.osm.gz --name "My city"

The page also reads a new area from any .osm file (Change area) and, when opened
from disk, fetches the live forecast itself.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from aamsim.osm_scene import load_scene
from aamsim.urbanwind import WindParams
from make_lab import REF


def rnd(a, q=0.5):
    return [round(float(x) / q) * q for x in a]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--osm", default="data/osm/dallas.osm.gz")
    ap.add_argument("--name", default="Downtown Dallas")
    ap.add_argument("--out", default="demo/corridoriq.html")
    a = ap.parse_args()
    sc = load_scene(a.osm, a.name)
    city = sc["city"]
    try:
        from aamsim.weather import fetch_hours
        wx = fetch_hours(*city.origin, hours=24)
        wx["source"] = wx["source"].split(" (")[0]
        wx["fetched"] = "captured at build time"
    except Exception as e:  # noqa: BLE001
        print("forecast unavailable:", e)
        wx = dict(source="default profile (no forecast at build time)", heights_m=[10, 80, 120, 180],
                  hours=[dict(time="", speed_ms=[3, 4.5, 5, 5.5], dir_from_deg=[180, 185, 190, 190])])
    data = dict(
        city=dict(name=a.name, origin=list(city.origin),
                  buildings=[dict(p=rnd(b.xy.ravel(), 0.1), h=round(b.height, 1), k=b.kind)
                             for b in city.buildings],
                  roads=[dict(p=rnd(r["p"].ravel()), w=r["w"], M=int(r["major"])) for r in sc["roads"]],
                  areas=[dict(p=rnd(x["p"].ravel()), k=x["kind"]) for x in sc["areas"]],
                  trees=[rnd(t, 0.5) for t in sc["trees"]],
                  places=[dict(n=p["name"], k=p["kind"], e=rnd(p["en"])) for p in sc["places"]]),
        weather=wx, windp=dict(WindParams().__dict__), ref=REF,
        rh=json.loads(Path("demo/rh_cases.json").read_text()))
    eng = Path("demo/lab_engine.js").read_text()
    # three.js is inlined so the single file also runs where outside scripts are blocked (SharePoint, Outlook, offline)
    vendor = {"https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js": "demo/vendor/three.min.js",
              "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js": "demo/vendor/OrbitControls.js"}
    html = (Path("demo/app_template.html").read_text()
            .replace("/*__ENGINE__*/", eng)
            .replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")).replace("</", "<\\/"))
            .replace('/*__ENGINE_SRC__*/""', json.dumps(eng).replace("</", "<\\/")))
    for url, f in vendor.items():
        tag = f'<script src="{url}"></script>'
        assert tag in html, url
        html = html.replace(tag, "<script>" + Path(f).read_text().replace("</script", "<\\/script") + "</script>")
    Path(a.out).write_text(html)
    if a.out == "demo/corridoriq.html":          # GitHub Pages serves docs/index.html
        Path("docs").mkdir(exist_ok=True)
        Path("docs/index.html").write_text(html)
    print(f"{a.out}: {len(html) / 1e6:.1f} MB · {len(data['city']['buildings'])} buildings · "
          f"{len(data['city']['places'])} places · {len(data['city']['trees'])} trees · {wx['source']}")


if __name__ == "__main__":
    main()
