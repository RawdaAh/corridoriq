"""Live or archived weather as a wind profile for the corridor.

Sources
  open-meteo : hourly wind speed and direction at 10, 80, 120 and 180 m above ground
               (no key; https://open-meteo.com). Best for the 30-120 m corridor layers.
  nws        : U.S. National Weather Service grid forecast, surface (10 m) wind only;
               extended upward with a power-law shear profile.

A profile is a dict:
  {"source", "time", "lat", "lon", "heights_m": [...], "speed_ms": [...], "dir_from_deg": [...]}
with meteorological directions (where the wind blows FROM, clockwise from north).

``layer_winds(profile, z_levels, bearing_deg)`` converts it to per-layer
(speed m/s, direction the wind blows TOWARD measured CCW from the corridor axis),
the convention used everywhere else in this package.
"""
from __future__ import annotations

import json
import math
import urllib.parse
import urllib.request

HEIGHTS = [10, 80, 120, 180]


def _get(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "aam-corridor-lab/0.1 (research)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def fetch_open_meteo(lat, lon, hour_index=None):
    q = dict(latitude=lat, longitude=lon, wind_speed_unit="ms", forecast_days=1, timezone="UTC",
             current="wind_speed_10m,wind_direction_10m",
             hourly=",".join([f"wind_speed_{h}m" for h in HEIGHTS] + [f"wind_direction_{h}m" for h in HEIGHTS]))
    d = _get("https://api.open-meteo.com/v1/forecast?" + urllib.parse.urlencode(q))
    if d.get("error"):
        raise RuntimeError(d.get("reason", "open-meteo error"))
    hr = d["hourly"]
    if hour_index is None:  # hour closest to 'current'
        now = d.get("current", {}).get("time", hr["time"][0])[:13]
        hour_index = next((i for i, t in enumerate(hr["time"]) if t[:13] == now), 0)
    return dict(source="Open-Meteo", time=hr["time"][hour_index] + "Z", lat=lat, lon=lon,
                heights_m=HEIGHTS,
                speed_ms=[hr[f"wind_speed_{h}m"][hour_index] for h in HEIGHTS],
                dir_from_deg=[hr[f"wind_direction_{h}m"][hour_index] for h in HEIGHTS])


def fetch_nws(lat, lon, alpha=0.25):
    pt = _get(f"https://api.weather.gov/points/{lat:.4f},{lon:.4f}")
    fc = _get(pt["properties"]["forecastHourly"])["properties"]["periods"][0]
    mph = float(str(fc["windSpeed"]).split()[0])
    compass = {"N": 0, "NNE": 22.5, "NE": 45, "ENE": 67.5, "E": 90, "ESE": 112.5, "SE": 135,
               "SSE": 157.5, "S": 180, "SSW": 202.5, "SW": 225, "WSW": 247.5, "W": 270,
               "WNW": 292.5, "NW": 315, "NNW": 337.5}
    d_from = compass.get(fc["windDirection"], 0.0)
    s10 = mph * 0.44704
    return dict(source="NWS (surface + power law)", time=fc["startTime"], lat=lat, lon=lon,
                heights_m=HEIGHTS, speed_ms=[s10 * (h / 10) ** alpha for h in HEIGHTS],
                dir_from_deg=[d_from] * len(HEIGHTS))


def fetch_profile(lat, lon, source="auto"):
    errors = []
    for s in (["open-meteo", "nws"] if source == "auto" else [source]):
        try:
            return fetch_open_meteo(lat, lon) if s == "open-meteo" else fetch_nws(lat, lon)
        except Exception as e:  # noqa: BLE001 - report and try the next source
            errors.append(f"{s}: {e}")
    raise RuntimeError("; ".join(errors))


def profile_at(profile, z):
    """Speed (m/s) and direction-from (deg) at height z: linear in log-height between
    levels, power law below the lowest level, vector-interpolated direction."""
    h = profile["heights_m"]
    s = profile["speed_ms"]
    d = profile["dir_from_deg"]
    z = max(float(z), 2.0)
    if z <= h[0]:
        return s[0] * (z / h[0]) ** 0.25, d[0]
    if z >= h[-1]:
        return s[-1], d[-1]
    k = max(i for i in range(len(h) - 1) if h[i] <= z)
    f = (math.log(z) - math.log(h[k])) / (math.log(h[k + 1]) - math.log(h[k]))
    sp = s[k] + f * (s[k + 1] - s[k])
    u = (1 - f) * math.sin(math.radians(d[k])) + f * math.sin(math.radians(d[k + 1]))
    v = (1 - f) * math.cos(math.radians(d[k])) + f * math.cos(math.radians(d[k + 1]))
    return sp, math.degrees(math.atan2(u, v)) % 360


def toward_relative(dir_from_deg, bearing_deg):
    """Meteorological 'from' (clockwise from north) -> 'toward', CCW from the corridor axis.
    ``bearing_deg`` is the corridor axis measured CCW from east."""
    toward_compass = (dir_from_deg + 180.0) % 360
    toward_math = 90.0 - toward_compass
    return (toward_math - bearing_deg) % 360


def layer_winds(profile, z_levels, bearing_deg):
    out = []
    for z in z_levels:
        sp, df = profile_at(profile, max(z, 10.0))
        out.append((sp, toward_relative(df, bearing_deg)))
    return out


def fetch_hours(lat, lon, hours=24):
    """Hourly wind profiles for the next ``hours`` hours (Open-Meteo, NWS fallback).
    Returns dict(source, hours=[{time, speed_ms[4], dir_from_deg[4]}], heights_m)."""
    try:
        q = dict(latitude=lat, longitude=lon, wind_speed_unit="ms", forecast_days=2, timezone="UTC",
                 hourly=",".join([f"wind_speed_{h}m" for h in HEIGHTS] + [f"wind_direction_{h}m" for h in HEIGHTS]))
        d = _get("https://api.open-meteo.com/v1/forecast?" + urllib.parse.urlencode(q))
        if d.get("error"):
            raise RuntimeError(d.get("reason"))
        hr = d["hourly"]
        import datetime as _dt
        now = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H")
        i0 = next((i for i, t in enumerate(hr["time"]) if t[:13] >= now), 0)
        out = [dict(time=hr["time"][i] + "Z", speed_ms=[hr[f"wind_speed_{h}m"][i] for h in HEIGHTS],
                    dir_from_deg=[hr[f"wind_direction_{h}m"][i] for h in HEIGHTS])
               for i in range(i0, min(i0 + hours, len(hr["time"])))]
        return dict(source="Open-Meteo forecast", heights_m=HEIGHTS, hours=out)
    except Exception as e1:  # noqa: BLE001
        pt = _get(f"https://api.weather.gov/points/{lat:.4f},{lon:.4f}")
        per = _get(pt["properties"]["forecastHourly"])["properties"]["periods"][:hours]
        compass = {"N": 0, "NNE": 22.5, "NE": 45, "ENE": 67.5, "E": 90, "ESE": 112.5, "SE": 135,
                   "SSE": 157.5, "S": 180, "SSW": 202.5, "SW": 225, "WSW": 247.5, "W": 270,
                   "WNW": 292.5, "NW": 315, "NNW": 337.5}
        out = []
        for p in per:
            s10 = float(str(p["windSpeed"]).split()[0]) * 0.44704
            d = compass.get(p["windDirection"], 0.0)
            out.append(dict(time=p["startTime"], speed_ms=[s10 * (h / 10) ** 0.25 for h in HEIGHTS],
                            dir_from_deg=[d] * len(HEIGHTS)))
        return dict(source=f"U.S. National Weather Service hourly forecast (Open-Meteo unavailable: {e1})",
                    heights_m=HEIGHTS, hours=out)
