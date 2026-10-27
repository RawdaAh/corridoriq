"""CorridorIQ build: placeholders filled, embedded data complete."""
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_app_builds(tmp_path):
    out = tmp_path / "app.html"
    subprocess.run([sys.executable, "make_app.py", "--out", str(out)], cwd=ROOT, check=True)
    html = out.read_text()
    assert "/*__ENGINE__*/" not in html and "/*__DATA__*/" not in html and "/*__ENGINE_SRC__*/" not in html
    data = json.loads(re.search(r"const DATA=(\{.*?\}); const ENGINE_SRC", html, re.S).group(1))
    city = data["city"]
    assert len(city["buildings"]) > 1000 and len(city["places"]) > 100
    names = {p["n"] for p in city["places"]}
    assert {"The Shed", "American Airlines Center"} <= names
    assert data["weather"]["hours"] and len(data["weather"]["hours"][0]["speed_ms"]) == 4
    assert len(data["rh"]["cases"]) >= 1 and "trr" in data["ref"]
