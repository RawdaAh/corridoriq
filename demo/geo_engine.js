/* AAM Corridor Lab - browser engine (twin of aamsim.urbanwind / geo_plan / sire).
   Pure functions, no DOM. Works in the browser (window.GeoEngine) and in node (module.exports). */
(function (root) {
  "use strict";
  const MPH_PER_MS = 2.236936;

  function layerSpeed(p, z) { return p.speed_ms * Math.pow(Math.max(z, 10) / p.ref_height, p.alpha); }

  // ---------------------------------------------------------------- urban wind (one layer)
  function solveLayer(solid, nx, ny, res, x0, y0, ux, uy, p) {
    const n = nx * ny, phi = new Float64Array(n);
    const id = (i, j) => i * ny + j;
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++)
      phi[id(i, j)] = ux * (x0 + (i + 0.5) * res) + uy * (y0 + (j + 0.5) * res);
    const fluid = (i, j) => !solid[id(i, j)];
    for (let it = 0; it < p.iterations; it++) {
      for (let colour = 0; colour < 2; colour++) {
        for (let i = 1; i < nx - 1; i++) {
          for (let j = 1 + ((i + 1 + colour) % 2 === 0 ? 0 : 0); j < ny - 1; j++) {
            if ((i + j) % 2 !== colour) continue;
            const c = id(i, j);
            if (solid[c]) continue;
            const pc = phi[c];
            const E = fluid(i + 1, j) ? phi[id(i + 1, j)] : pc;
            const W = fluid(i - 1, j) ? phi[id(i - 1, j)] : pc;
            const N = fluid(i, j + 1) ? phi[id(i, j + 1)] : pc;
            const S = fluid(i, j - 1) ? phi[id(i, j - 1)] : pc;
            const avg = 0.25 * (E + W + N + S);
            phi[c] = pc + p.omega * (avg - pc);
          }
        }
      }
    }
    const u = new Float64Array(n), v = new Float64Array(n);
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
      const c = id(i, j);
      if (solid[c]) continue;
      if (i === 0 || i === nx - 1 || j === 0 || j === ny - 1) { u[c] = ux; v[c] = uy; continue; }
      const bE = fluid(i + 1, j), bW = fluid(i - 1, j), bN = fluid(i, j + 1), bS = fluid(i, j - 1);
      const pc = phi[c];
      u[c] = bE && bW ? (phi[id(i + 1, j)] - phi[id(i - 1, j)]) / (2 * res)
        : bE ? (phi[id(i + 1, j)] - pc) / res : bW ? (pc - phi[id(i - 1, j)]) / res : 0;
      v[c] = bN && bS ? (phi[id(i, j + 1)] - phi[id(i, j - 1)]) / (2 * res)
        : bN ? (phi[id(i, j + 1)] - pc) / res : bS ? (pc - phi[id(i, j - 1)]) / res : 0;
    }
    const U = Math.hypot(ux, uy);
    let anySolid = false; for (let c = 0; c < n; c++) if (solid[c]) { anySolid = true; break; }
    if (U > 1e-9 && anySolid) {
      const dx = ux / U, dy = uy / U;
      const sh = new Float64Array(n);
      for (let k = 1; k <= p.wake_reach; k++) {
        const w = Math.exp(-k * res / p.wake_len);
        for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
          const si = Math.floor(i - dx * k + 0.5), sj = Math.floor(j - dy * k + 0.5);
          if (si < 0 || si >= nx || sj < 0 || sj >= ny) continue;
          if (solid[id(si, sj)] && w > sh[id(i, j)]) sh[id(i, j)] = w;
        }
      }
      const bl = new Float64Array(n);
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
        let s = 0;
        for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) {
          const ii = Math.min(nx - 1, Math.max(0, i + a)), jj = Math.min(ny - 1, Math.max(0, j + b));
          s += sh[id(ii, jj)];
        }
        bl[id(i, j)] = s / 9;
      }
      for (let c = 0; c < n; c++) { const f = 1 - p.wake_depth * bl[c]; u[c] *= f; v[c] *= f; }
    }
    const cap = p.max_speedup * U;
    for (let c = 0; c < n; c++) {
      if (solid[c]) { u[c] = 0; v[c] = 0; continue; }
      const sp = Math.hypot(u[c], v[c]);
      if (sp > cap) { const s = cap / Math.max(sp, 1e-12); u[c] *= s; v[c] *= s; }
    }
    return { u, v };
  }

  function solveWind(grid, zLevels, p) {
    const zs = zLevels.map((z) => Math.max(z, 10));
    const th = p.direction_deg * Math.PI / 180, layers = [];
    for (const z of zs) {
      const s = layerSpeed(p, z), solid = new Uint8Array(grid.nx * grid.ny);
      for (let c = 0; c < solid.length; c++) solid[c] = grid.H[c] >= z ? 1 : 0;
      layers.push(solveLayer(solid, grid.nx, grid.ny, grid.res, grid.x0, grid.y0, s * Math.cos(th), s * Math.sin(th), p));
    }
    return { z: zs, layers, nx: grid.nx, ny: grid.ny, res: grid.res, x0: grid.x0, y0: grid.y0 };
  }

  function sample(W, x, y, z, out) {
    const nx = W.nx, ny = W.ny, nl = W.z.length;
    const fx = Math.min(Math.max((x - W.x0) / W.res - 0.5, 0), nx - 1.000001);
    const fy = Math.min(Math.max((y - W.y0) / W.res - 0.5, 0), ny - 1.000001);
    const i = Math.floor(fx), j = Math.floor(fy), ax = fx - i, ay = fy - j;
    let k = 0;
    while (k < nl - 2 && W.z[k + 1] < z) k++;           // == searchsorted(z) - 1, clipped
    if (nl > 1 && z <= W.z[0]) k = 0;
    const az = nl > 1 ? Math.min(Math.max((z - W.z[k]) / (W.z[k + 1] - W.z[k]), 0), 1) : 0;
    const bil = (F) => (1 - ax) * (1 - ay) * F[i * ny + j] + ax * (1 - ay) * F[(i + 1) * ny + j]
      + (1 - ax) * ay * F[i * ny + j + 1] + ax * ay * F[(i + 1) * ny + j + 1];
    const L0 = W.layers[k], L1 = W.layers[Math.min(k + 1, nl - 1)];
    out = out || [0, 0];
    out[0] = (1 - az) * bil(L0.u) + az * bil(L1.u);
    out[1] = (1 - az) * bil(L0.v) + az * bil(L1.v);
    return out;
  }

  // ---------------------------------------------------------------- corridor & beliefs
  function makeCorridor(C) {
    const E = C.edges.length, A = 2 * E;
    const af = new Int32Array(A), at = new Int32Array(A), edge = new Int32Array(A), len = new Float64Array(A);
    for (let e = 0; e < E; e++) {
      const [u, v] = C.edges[e];
      af[e] = u; at[e] = v; af[e + E] = v; at[e + E] = u; edge[e] = e; edge[e + E] = e;
    }
    for (let a = 0; a < A; a++) {
      const p = C.pos[af[a]], q = C.pos[at[a]];
      len[a] = Math.hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2]);
    }
    const out = Array.from({ length: C.pos.length }, () => []);
    for (let a = 0; a < A; a++) out[af[a]].push(a);    // ascending arc id, like argsort(stable)
    return Object.assign({}, C, { af, at, edge, len, out, A, ground: C.pos.map((p) => p[2] === 0) });
  }

  function nodeWind(W, corr) {      // mph at nodes
    const u = new Float64Array(corr.pos.length), v = new Float64Array(corr.pos.length), o = [0, 0];
    corr.pos.forEach((p, n) => { sample(W, p[0], p[1], p[2], o); u[n] = o[0] * MPH_PER_MS; v[n] = o[1] * MPH_PER_MS; });
    return { u, v };
  }

  function beliefs(W, corr, p, dirErr = 15, bias = -0.15) {
    const truth = nodeWind(W, corr), N = corr.pos.length;
    const fu = new Float64Array(N), fv = new Float64Array(N), th = (p.direction_deg + dirErr) * Math.PI / 180;
    corr.pos.forEach((q, n) => { const s = layerSpeed(p, q[2]) * (1 + bias) * MPH_PER_MS; fu[n] = s * Math.cos(th); fv[n] = s * Math.sin(th); });
    const pu = new Float64Array(N), pv = new Float64Array(N), sums = {};
    corr.ijk.forEach((k, n) => { const l = k[2]; (sums[l] = sums[l] || [0, 0, 0]); sums[l][0] += truth.u[n]; sums[l][1] += truth.v[n]; sums[l][2]++; });
    corr.ijk.forEach((k, n) => { const s = sums[k[2]]; pu[n] = s[0] / s[2]; pv[n] = s[1] / s[2]; });
    return { ground_truth: truth, forecast: { u: fu, v: fv }, profile: { u: pu, v: pv },
      wind_agnostic: { u: new Float64Array(N), v: new Float64Array(N) } };
  }

  // TRR link projection with the 'crab' rule, vehicle in metre geometry
  function linkCosts(corr, w, veh) {
    const A = corr.A, t = new Float64Array(A), e = new Float64Array(A);
    const vU = veh.cruise_mph, fps = vU * 5280 / 3600;
    for (let a = 0; a < A; a++) {
      const p = corr.pos[corr.af[a]], q = corr.pos[corr.at[a]];
      const dxh = q[0] - p[0], dyh = q[1] - p[1], horiz = Math.hypot(dxh, dyh);
      const L = corr.len[a] * veh.length_unit_ft;
      let vg = vU;
      if (horiz >= 1e-9) {
        const theta = ((Math.atan2(dyh, dxh) * 180 / Math.PI) % 360 + 360) % 360;
        const ub = 0.5 * (w.u[corr.af[a]] + w.u[corr.at[a]]), vb = 0.5 * (w.v[corr.af[a]] + w.v[corr.at[a]]);
        const ws = Math.hypot(ub, vb), alpha = ((Math.atan2(vb, ub) * 180 / Math.PI) % 360 + 360) % 360;
        let g = ((alpha - theta) % 360 + 360) % 360; if (g > 180) g = 360 - g;
        const psi = ws * Math.cos(g * Math.PI / 180), cross = ws * Math.sin(g * Math.PI / 180);
        vg = Math.sqrt(Math.max(vU * vU - cross * cross, 0)) + psi;
      }
      const lam = vg >= 1 ? vU / vg : Infinity;
      t[a] = L / fps * lam; e[a] = veh.energy_j_per_ft * L * lam * lam;
    }
    return { t, e };
  }

  // binary heap on [k1, k2, node]
  function Heap() { this.a = []; }
  Heap.prototype.less = function (x, y) { return x[0] < y[0] || (x[0] === y[0] && (x[1] < y[1] || (x[1] === y[1] && x[2] < y[2]))); };
  Heap.prototype.push = function (x) { const a = this.a; a.push(x); let i = a.length - 1; while (i > 0) { const p = (i - 1) >> 1; if (this.less(a[i], a[p])) { [a[i], a[p]] = [a[p], a[i]]; i = p; } else break; } };
  Heap.prototype.pop = function () { const a = this.a, top = a[0], last = a.pop(); if (a.length) { a[0] = last; let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i; if (l < a.length && this.less(a[l], a[m])) m = l; if (r < a.length && this.less(a[r], a[m])) m = r; if (m === i) break; [a[i], a[m]] = [a[m], a[i]]; i = m; } } return top; };

  function minEnergyRoute(corr, e, t, o, d) {
    const N = corr.pos.length, bE = new Float64Array(N).fill(Infinity), bT = new Float64Array(N).fill(Infinity), pred = new Int32Array(N).fill(-1);
    bE[o] = 0; bT[o] = 0;
    const h = new Heap(); h.push([0, 0, o]);
    while (h.a.length) {
      const [ce, ct, u] = h.pop();
      if (ce > bE[u] || (ce === bE[u] && ct > bT[u])) continue;
      if (u === d) break;
      for (const a of corr.out[u]) {
        const v = corr.at[a];
        if (corr.ground[v] && v !== d) continue;
        if (!isFinite(e[a])) continue;
        const ne = ce + e[a], nt = ct + t[a];
        if (ne < bE[v] || (ne === bE[v] && nt < bT[v])) { bE[v] = ne; bT[v] = nt; pred[v] = a; h.push([ne, nt, v]); }
      }
    }
    if (!isFinite(bE[d])) return null;
    const arcs = []; let n = d;
    while (n !== o) { const a = pred[n]; arcs.push(a); n = corr.af[a]; }
    return arcs.reverse();
  }

  function fly(corr, arcs, c) {
    const T = [0], E = [0];
    for (const a of arcs) { T.push(T[T.length - 1] + c.t[a]); E.push(E[E.length - 1] + c.e[a]); }
    return { T, E };
  }

  const METHODS = ["ground_truth", "forecast", "profile", "wind_agnostic"];
  function planMission(corr, W, p, veh, o, d) {
    const B = beliefs(W, corr, p), truthC = linkCosts(corr, B.ground_truth, veh), out = {};
    for (const m of METHODS) {
      const c = linkCosts(corr, B[m], veh), arcs = minEnergyRoute(corr, c.e, c.t, o, d);
      if (!arcs) { out[m] = null; continue; }
      const plan = fly(corr, arcs, c), real = fly(corr, arcs, truthC);
      const nodes = [corr.af[arcs[0]], ...arcs.map((a) => corr.at[a])];
      out[m] = { arcs, nodes, planT: plan.T, planE: plan.E, T: real.T, E: real.E,
        dist: arcs.reduce((s, a) => s + corr.len[a], 0) };
    }
    return { routes: out, beliefs: B, truthCosts: truthC };
  }

  // ---------------------------------------------------------------- SIRE registry
  function sire(corr, tauSec, slotSec, H, trips, rule) {
    const A = corr.A, q = new Int32Array(A);
    for (let a = 0; a < A; a++) q[a] = Math.max(1, Math.ceil(tauSec[a] / slotSec - 1e-9));
    const nk = rule === "directed" ? A : corr.edges.length, key = (a) => rule === "directed" ? a : corr.edge[a];
    let qmax = 1; for (let a = 0; a < A; a++) if (isFinite(tauSec[a]) && q[a] > qmax) qmax = q[a];
    const W = H + qmax + 2, blocked = new Uint8Array(nk * W), N = corr.pos.length;
    const order = trips.map((t, i) => i).sort((i, j) => trips[i].dispatch - trips[j].dispatch || trips[i].op - trips[j].op || trips[i].id - trips[j].id);
    const best = new Float64Array((H + 1) * N), pred = new Int32Array((H + 1) * N), predT = new Int32Array((H + 1) * N);
    const results = new Array(trips.length);
    for (const ti of order) {
      const tr = trips[ti], t0 = tr.dispatch, o = tr.o, d = tr.d;
      best.fill(Infinity); pred.fill(-1);
      best[t0 * N + o] = 0;
      let bd = Infinity, bt = -1, fe = t0;
      for (let t = t0; t < H && t <= fe; t++) {
        const c0 = best[t * N + d]; if (c0 < bd) { bd = c0; bt = t; }
        for (let u = 0; u < N; u++) {
          const cu = best[t * N + u];
          if (cu === Infinity || cu >= bd || u === d) continue;
          for (const a of corr.out[u]) {
            if (!isFinite(tauSec[a])) continue;
            const v = corr.at[a];
            if (corr.ground[v] && v !== d) continue;
            const ta = t + q[a]; if (ta >= H) continue;
            const k = key(a) * W; let free = true;
            for (let s = t; s < ta; s++) if (blocked[k + s]) { free = false; break; }
            if (!free) continue;
            const cand = cu + tauSec[a];
            if (cand < best[ta * N + v]) { best[ta * N + v] = cand; pred[ta * N + v] = a; predT[ta * N + v] = t; if (ta > fe) fe = ta; }
          }
        }
      }
      for (let t = t0; t <= H; t++) if (best[t * N + d] < bd) { bd = best[t * N + d]; bt = t; }
      if (!isFinite(bd)) { results[ti] = { ok: false }; continue; }
      const legs = []; let n = d, t = bt;
      while (!(n === o && t === t0)) { const a = pred[t * N + n], s = predT[t * N + n]; legs.push([a, s, t]); n = corr.af[a]; t = s; }
      legs.reverse();
      for (const [a, s, e] of legs) for (let x = s; x < e; x++) blocked[key(a) * W + x] = 1;
      results[ti] = { ok: true, legs, travel: bd };
    }
    return results;
  }

  const api = { MPH_PER_MS, layerSpeed, solveLayer, solveWind, sample, makeCorridor, nodeWind, beliefs,
    linkCosts, minEnergyRoute, fly, planMission, sire, METHODS };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.GeoEngine = api;
})(typeof window !== "undefined" ? window : globalThis);
