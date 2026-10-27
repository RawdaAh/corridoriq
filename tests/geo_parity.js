// node tests/geo_parity.js  -> prints JSON with parity stats (used by tests/test_geo.py)
const fs = require("fs");
const G = require("../demo/geo_engine.js");
const D = JSON.parse(fs.readFileSync(__dirname + "/data/geo_data.json"));
const R = JSON.parse(fs.readFileSync(__dirname + "/data/geo_reference.json"));
const buf = Buffer.from(D.grid.H, "base64");
const H = new Float64Array(D.grid.nx * D.grid.ny);
for (let i = 0; i < H.length; i++) H[i] = buf.readUInt16LE(2 * i) / 10;
const grid = Object.assign({}, D.grid, { H });
const p = Object.assign({}, D.wind, R.params);
const corr = G.makeCorridor(D.corridor);
const t0 = Date.now();
const W = G.solveWind(grid, D.corridor.z_levels, p);
const tw = Date.now() - t0;
const nw = G.nodeWind(W, corr);
let du = 0;
for (let n = 0; n < nw.u.length; n++) du = Math.max(du, Math.abs(nw.u[n] - R.node_u[n]), Math.abs(nw.v[n] - R.node_v[n]));
const plan = G.planMission(corr, W, p, D.vehicle, R.o, R.d);
const routes = {};
for (const m of G.METHODS) {
  const r = plan.routes[m], ref = R.routes[m];
  routes[m] = { same_nodes: JSON.stringify(r.nodes) === JSON.stringify(ref.nodes),
    dPlannedE_rel: Math.abs(r.planE[r.planE.length - 1] - ref.PE) / ref.PE,
    dE: Math.abs(r.E[r.E.length - 1] - ref.E), dT: Math.abs(r.T[r.T.length - 1] - ref.T) };
}
const trips = R.trips.map((t) => ({ id: t.trip_id, op: t.operator, o: t.origin, d: t.destination, dispatch: t.dispatch }));
const s = G.sire(corr, plan.truthCosts.t, 10, 120, trips, "undirected");
const acc = s.map((r) => (r.ok ? 1 : 0));
const same = acc.every((a, i) => a === R.sire_accepted[i]);
console.log(JSON.stringify({ wind_ms: tw, max_node_wind_diff_mph: du, routes, sire_same: same,
  sire_accepted: acc.reduce((a, b) => a + b, 0) }));
