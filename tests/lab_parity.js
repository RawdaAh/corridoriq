// node tests/lab_parity.js -> JSON with browser-engine reproduction checks (read by tests/test_lab.py)
const fs = require("fs");
const G = require("../demo/lab_engine.js");
const FT = G.FT, out = {};
// ---- Paper 1: TRR corridor (5 x 6,000 ft, 100 ft lanes, layers 0/200/300/400 ft, diagonal take-off)
const trr = G.buildCorridor({ L: 30000 / FT, module_length: 6000 / FT, module_width: 100 / FT, n_lat: 2,
  z_levels: [0, 200, 300, 400].map((z) => z / FT), vertical_ground_only: false }, null);
const o = trr.node(0, 0, 0), d = trr.node(5, 0, 0);
const veh = { cruise_mph: 35, energy_j_per_ft: 20.42, phase_dependent: false, k_climb: 0.2, alpha_desc: 1, takeoff_fraction: 0.1 };
const run = (w, v) => { const c = G.linkCosts(trr, w, v); const r = G.pareto(trr, c, w, v, o, d).routes;
  return { n: r.length, T: r.reduce((s, x) => s + x.time, 0) / r.length / 60, E: r.reduce((s, x) => s + x.energy, 0) / r.length / 3.6e6 }; };
out.trr = { nodes: trr.pos.length, none: run(G.windNone(trr), veh), tail: run(G.windShear(trr, 0), veh), d45: run(G.windShear(trr, 45), veh),
  head: run(G.windShear(trr, 180), veh), spatial_tail: run(G.windSpatial(trr, 0), veh),
  phase_none: run(G.windNone(trr), Object.assign({}, veh, { phase_dependent: true })),
  phase_head: run(G.windShear(trr, 180), Object.assign({}, veh, { phase_dependent: true })) };
// temporal: depart 0 and 27 min, tailwind calm->strong
const wt = G.windTemporal(trr, 0, "calm_to_strong"), ct = G.linkCosts(trr, wt, veh);
out.trr.temporal = [0, 9].map((k) => { const r = G.pareto(trr, ct, wt, veh, o, d, { depart: 180 * k }).routes; return r.reduce((s, x) => s + x.time, 0) / r.length / 60; });
// ---- Paper 2: TRB corridor, trips exported by Python (D = 2000, seed 1)
const trb = G.buildCorridor({ L: 30000 / FT, module_length: 1000 / FT, module_width: 200 / FT, n_lat: 2, z_levels: [0, 100, 200, 300, 400].map((z) => z / FT), vertical_ground_only: false }, null);
const ref = JSON.parse(fs.readFileSync(__dirname + "/data/sire_reference.json"));
const tau = Array.from(trb.len, (L) => L * FT / (35 * 5280 / 3600));
const res = G.sire(trb, tau, ref.trips, { slot: 60, H: 120, rule: "undirected", groundTransit: true });
const k = G.kpis(ref.trips, res);
out.trb = { nodes: trb.pos.length, edges: trb.edges.length, sr: k.sr, accepted: k.accepted, tt: k.tt, gini: k.gini, cascade: k.cascade,
  same_decisions: res.every((r, i) => (r.ok ? 1 : 0) === ref.accepted[i]) };
console.log(JSON.stringify(out));
