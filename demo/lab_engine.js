/* AAM Corridor Lab - browser engine. Line-for-line twin of the Python package:
     corridor.py / geo.py      -> buildCorridor, rasterize, parseOSM
     wind.py (TRR experiments) -> windNone, windShear, windSpatial, windTemporal
     urbanwind.py              -> solveLayer, solveUrban, sample
     energy.py + pareto.py     -> linkCosts (TRR Eqs. 10-33), pareto (Algorithms 1-2)
     sire.py                   -> demand, sire, kpis (TRB Eqs. 1-19)
   Units: geometry in metres, wind at nodes in mph (as in the TRR equations), vehicle in mph and J/ft.
   No DOM access except parseOSM (DOMParser). Works in the browser, a Web Worker, and node. */
(function (root) {
  "use strict";
  const FT = 3.28084, MPH_PER_MS = 2.236936, TOL = 1e-6;

  // ================================================================= geometry
  function corridorFrame(startEN, goalEN) {
    const dx = goalEN[0] - startEN[0], dy = goalEN[1] - startEN[1];
    return { s: startEN, b: Math.atan2(dy, dx), L: Math.hypot(dx, dy) };
  }
  function toCorr(F, x, y) { const c = Math.cos(F.b), s = Math.sin(F.b), dx = x - F.s[0], dy = y - F.s[1]; return [dx * c + dy * s, -dx * s + dy * c]; }
  function toEN(F, x, y) { const c = Math.cos(F.b), s = Math.sin(F.b); return [x * c - y * s + F.s[0], x * s + y * c + F.s[1]]; }

  function pip(px, py, poly) {      // poly: flat [x0,y0,x1,y1,...]
    let inside = false; const n = poly.length / 2;
    for (let i = 0, j = n - 1; i < n; j = i++) {
      const xi = poly[2 * i], yi = poly[2 * i + 1], xj = poly[2 * j], yj = poly[2 * j + 1];
      if ((yi > py) !== (yj > py) && px < (xj - xi) * (py - yi) / (yj - yi + 1e-12) + xi) inside = !inside;
    }
    return inside;
  }
  const bisectLeft = (arr0, res, n, v) => { // first k with arr0 + (k + .5) res >= v
    let lo = 0, hi = n; while (lo < hi) { const m = (lo + hi) >> 1; if (arr0 + (m + 0.5) * res < v) lo = m + 1; else hi = m; } return lo; };

  // buildings: [{p:[x,y,...] (corridor frame), h}] -> H (Float64Array nx*ny), building-top height
  function rasterize(buildings, x0, y0, nx, ny, res) {
    const H = new Float64Array(nx * ny);
    for (const b of buildings) {
      const p = b.p; let lx = Infinity, ly = Infinity, hx = -Infinity, hy = -Infinity, mx = 0, my = 0;
      for (let k = 0; k < p.length; k += 2) { lx = Math.min(lx, p[k]); hx = Math.max(hx, p[k]); ly = Math.min(ly, p[k + 1]); hy = Math.max(hy, p[k + 1]); mx += p[k]; my += p[k + 1]; }
      const i0 = bisectLeft(x0, res, nx, lx), i1 = bisectLeft(x0, res, nx, hx), j0 = bisectLeft(y0, res, ny, ly), j1 = bisectLeft(y0, res, ny, hy);
      if (i1 <= i0 || j1 <= j0) {
        const m = p.length / 2;
        const i = Math.min(nx - 1, Math.max(0, Math.floor((mx / m - x0) / res))), j = Math.min(ny - 1, Math.max(0, Math.floor((my / m - y0) / res)));
        H[i * ny + j] = Math.max(H[i * ny + j], b.h); continue;
      }
      for (let i = i0; i < i1; i++) for (let j = j0; j < j1; j++)
        if (pip(x0 + (i + 0.5) * res, y0 + (j + 0.5) * res, p)) H[i * ny + j] = Math.max(H[i * ny + j], b.h);
    }
    for (let c = 0; c < H.length; c++) H[c] = Math.round(H[c] * 10) / 10;
    return H;
  }

  /* spec: {L, module_length, module_width, n_lat, z_levels (m), vertical_ground_only, clearance,
            res, margin, portal_depth}; buildings in corridor frame or null (open field). */
  function buildCorridor(spec, buildings) {
    const L = spec.L, NL = Math.max(1, Math.round(L / spec.module_length)), ML = L / NL;
    const NW = spec.n_lat, NH = spec.z_levels.length - 1, W = spec.module_width, half = NW * W / 2;
    const nid = (i, j, k) => (i * (NW + 1) + j) * (NH + 1) + k;
    const N = (NL + 1) * (NW + 1) * (NH + 1), pos = [], ijk = [];
    for (let i = 0; i <= NL; i++) for (let j = 0; j <= NW; j++) for (let k = 0; k <= NH; k++) { pos.push([i * ML, j * W - half, spec.z_levels[k]]); ijk.push([i, j, k]); }
    const set = new Set(), edges = [];
    for (let i = 0; i < NL; i++) for (let j = 0; j < NW; j++) for (let k = 0; k < NH; k++) {
      const cs = []; for (const a of [0, 1]) for (const b of [0, 1]) for (const c of [0, 1]) cs.push(nid(i + a, j + b, k + c));
      for (let x = 0; x < 8; x++) for (let y = x + 1; y < 8; y++) {
        const u = cs[x], v = cs[y];
        if (pos[u][2] === 0 && pos[v][2] === 0) continue;
        if (spec.vertical_ground_only && (pos[u][2] === 0 || pos[v][2] === 0) && !(Math.abs(pos[u][0] - pos[v][0]) < 1e-9 && Math.abs(pos[u][1] - pos[v][1]) < 1e-9)) continue;
        const a = Math.min(u, v), b = Math.max(u, v), key = a * N + b;
        if (!set.has(key)) { set.add(key); edges.push([a, b]); }
      }
    }
    edges.sort((p, q) => p[0] - q[0] || p[1] - q[1]);
    const ground = pos.map((p) => p[2] === 0), depth = spec.portal_depth || 1;
    let keep = edges.map(() => true), grid = null, removed = 0;
    if (buildings && buildings.length) {
      const res = spec.res || 10, margin = spec.margin || 200, clr = spec.clearance ?? 8;
      const x0 = -margin, y0 = -half - margin, nx = Math.floor((L + 2 * margin) / res) + 1, ny = Math.floor((2 * half + 2 * margin) / res) + 1;
      const H = rasterize(buildings, x0, y0, nx, ny, res), pad = Math.ceil(clr / res), Hp = new Float64Array(nx * ny);
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
        let m = H[i * ny + j];
        for (let a = -pad; a <= pad; a++) for (let b = -pad; b <= pad; b++) {
          const ii = ((i - a) % nx + nx) % nx, jj = ((j - b) % ny + ny) % ny;   // np.roll semantics
          if (H[ii * ny + jj] > m) m = H[ii * ny + jj];
        }
        Hp[i * ny + j] = m;
      }
      const cell = (x, y) => [Math.min(nx - 1, Math.max(0, Math.trunc((x - x0) / res))), Math.min(ny - 1, Math.max(0, Math.trunc((y - y0) / res)))];
      const blocked = (x, y, z) => { const [i, j] = cell(x, y); return z <= Hp[i * ny + j] + clr; };
      const nodeOk = pos.map((p, n) => { if (ground[n]) { const [i, j] = cell(p[0], p[1]); return Hp[i * ny + j] <= 0; } return !blocked(p[0], p[1], Math.max(p[2], 0.5)); });
      keep = edges.map(([u, v]) => {
        if (!nodeOk[u] || !nodeOk[v]) return false;
        const P0 = pos[u], P1 = pos[v], vert = Math.abs(P1[0] - P0[0]) < 1e-6 && Math.abs(P1[1] - P0[1]) < 1e-6;
        for (let s = 0; s < 9; s++) {
          const f = 0.1 + 0.1 * s, x = P0[0] + (P1[0] - P0[0]) * f, y = P0[1] + (P1[1] - P0[1]) * f, z = P0[2] + (P1[2] - P0[2]) * f;
          if (!vert && blocked(x, y, Math.max(z, 0.5))) return false;
          if (vert) { const [i, j] = cell(x, y); if (Hp[i * ny + j] > 0) return false; }
        }
        return true;
      });
      grid = { x0, y0, nx, ny, res, H, half, L };
    }
    const kept = edges.filter((e, k) => keep[k]);
    removed = edges.length - kept.length;
    const used = new Uint8Array(N); for (const [u, v] of kept) { used[u] = 1; used[v] = 1; }
    const nidMap = new Int32Array(N).fill(-1); let c = 0; for (let n = 0; n < N; n++) if (used[n]) nidMap[n] = c++;
    const P2 = [], IJK = [], origins = [], destinations = [];
    for (let n = 0; n < N; n++) if (used[n]) {
      P2.push(pos[n]); IJK.push(ijk[n]);
      if (ground[n] && ijk[n][0] < depth) origins.push(nidMap[n]);
      if (ground[n] && ijk[n][0] > NL - depth) destinations.push(nidMap[n]);
    }
    return finish({ pos: P2, ijk: IJK, edges: kept.map(([u, v]) => [nidMap[u], nidMap[v]]), origins, destinations,
      NL, NW, NH, half, L, ML, W, z_levels: spec.z_levels.slice(), grid, removed_links: removed, full_links: edges.length });
  }
  function finish(C) {
    const E = C.edges.length, A = 2 * E;
    const af = new Int32Array(A), at = new Int32Array(A), edge = new Int32Array(A), len = new Float64Array(A);
    for (let e = 0; e < E; e++) { const [u, v] = C.edges[e]; af[e] = u; at[e] = v; af[e + E] = v; at[e + E] = u; edge[e] = e; edge[e + E] = e; }
    for (let a = 0; a < A; a++) { const p = C.pos[af[a]], q = C.pos[at[a]]; len[a] = Math.hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2]); }
    const out = Array.from({ length: C.pos.length }, () => []);
    for (let a = 0; a < A; a++) out[af[a]].push(a);
    const node = (i, j, k) => C.ijk.findIndex((q) => q[0] === i && q[1] === j && q[2] === k);
    return Object.assign(C, { af, at, edge, len, out, A, ground: C.pos.map((p) => p[2] === 0), node });
  }

  // ================================================================= wind at nodes (mph), bins
  const field = (C, nb, bin_sec, fn, name) => {
    const bins = [];
    for (let b = 0; b < nb; b++) { const u = new Float64Array(C.pos.length), v = new Float64Array(C.pos.length);
      C.pos.forEach((p, n) => { const [s, d] = fn(p, b * bin_sec); u[n] = s * Math.cos(d * Math.PI / 180); v[n] = s * Math.sin(d * Math.PI / 180); });
      bins.push({ u, v }); }
    return { bins, bin_sec, name };
  };
  const zft = (z) => z * FT;
  const shear = (zf) => zf < 200 ? 10 : zf < 300 ? 15 : zf < 400 ? 20 : 30;     // TRR 6.2 (mph)
  const weak = (zf) => zf < 200 ? 10 : zf < 300 ? 12 : zf < 400 ? 13 : 15;      // TRR 6.3 B1-B2
  const windNone = (C) => field(C, 1, 180, () => [0, 0], "none");
  const windShear = (C, dir) => field(C, 1, 180, (p) => [shear(zft(p[2])), dir], `shear ${dir}`);
  const windSpatial = (C, dir, nBlocks = 5) => field(C, 1, 180, (p) => {
    const blk = Math.min(Math.floor(p[0] / C.L * nBlocks), nBlocks - 1);
    return [blk === 0 || blk === 1 ? weak(zft(p[2])) : shear(zft(p[2])), dir]; }, `spatial ${dir}`);
  const windTemporal = (C, dir, trend, nb = 20, bin = 180, win = 10, lo = 2, hi = 22) => field(C, nb, bin, (p, t) => {
    const f = Math.min(t / (bin * (win - 1)), 1); return [trend === "calm_to_strong" ? lo + (hi - lo) * f : hi - (hi - lo) * f, dir]; }, `${trend} ${dir}`);
  function windFromUrban(C, W) {   // one bin, mph
    const u = new Float64Array(C.pos.length), v = new Float64Array(C.pos.length), o = [0, 0];
    C.pos.forEach((p, n) => { sample(W, p[0], p[1], p[2], o); u[n] = o[0] * MPH_PER_MS; v[n] = o[1] * MPH_PER_MS; });
    return { bins: [{ u, v }], bin_sec: 3600, name: "urban" };
  }

  // ================================================================= urban wind surrogate
  function solveLayer(solid, nx, ny, res, x0, y0, ux, uy, p) {
    const n = nx * ny, phi = new Float64Array(n), id = (i, j) => i * ny + j, fl = (i, j) => !solid[id(i, j)];
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) phi[id(i, j)] = ux * (x0 + (i + 0.5) * res) + uy * (y0 + (j + 0.5) * res);
    for (let it = 0; it < p.iterations; it++) for (let colour = 0; colour < 2; colour++)
      for (let i = 1; i < nx - 1; i++) for (let j = 1 + ((i + 1 + colour) & 1); j < ny - 1; j += 2) {
        const c = id(i, j); if (solid[c]) continue; const pc = phi[c];
        const E = fl(i + 1, j) ? phi[id(i + 1, j)] : pc, Wv = fl(i - 1, j) ? phi[id(i - 1, j)] : pc;
        const N = fl(i, j + 1) ? phi[id(i, j + 1)] : pc, S = fl(i, j - 1) ? phi[id(i, j - 1)] : pc;
        phi[c] = pc + p.omega * (0.25 * (E + Wv + N + S) - pc);
      }
    const u = new Float64Array(n), v = new Float64Array(n);
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
      const c = id(i, j); if (solid[c]) continue;
      if (i === 0 || i === nx - 1 || j === 0 || j === ny - 1) { u[c] = ux; v[c] = uy; continue; }
      const bE = fl(i + 1, j), bW = fl(i - 1, j), bN = fl(i, j + 1), bS = fl(i, j - 1), pc = phi[c];
      u[c] = bE && bW ? (phi[id(i + 1, j)] - phi[id(i - 1, j)]) / (2 * res) : bE ? (phi[id(i + 1, j)] - pc) / res : bW ? (pc - phi[id(i - 1, j)]) / res : 0;
      v[c] = bN && bS ? (phi[id(i, j + 1)] - phi[id(i, j - 1)]) / (2 * res) : bN ? (phi[id(i, j + 1)] - pc) / res : bS ? (pc - phi[id(i, j - 1)]) / res : 0;
    }
    const U = Math.hypot(ux, uy); let any = false; for (let c = 0; c < n; c++) if (solid[c]) { any = true; break; }
    if (U > 1e-9 && any) {
      const dx = ux / U, dy = uy / U, sh = new Float64Array(n);
      for (let k = 1; k <= p.wake_reach; k++) { const w = Math.exp(-k * res / p.wake_len);
        for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) { const si = Math.floor(i - dx * k + 0.5), sj = Math.floor(j - dy * k + 0.5);
          if (si < 0 || si >= nx || sj < 0 || sj >= ny) continue; if (solid[id(si, sj)] && w > sh[id(i, j)]) sh[id(i, j)] = w; } }
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) { let s = 0;
        for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) s += sh[id(Math.min(nx - 1, Math.max(0, i + a)), Math.min(ny - 1, Math.max(0, j + b)))];
        const f = 1 - p.wake_depth * s / 9; u[id(i, j)] *= f; v[id(i, j)] *= f; }
    }
    const cap = p.max_speedup * U;
    for (let c = 0; c < n; c++) { if (solid[c]) { u[c] = 0; v[c] = 0; continue; } const s = Math.hypot(u[c], v[c]); if (s > cap) { const k = cap / Math.max(s, 1e-12); u[c] *= k; v[c] *= k; } }
    return { u, v };
  }
  const layerSpeed = (p, z) => p.speed_ms * Math.pow(Math.max(z, 10) / p.ref_height, p.alpha);
  // layers: optional [[speed m/s, toward deg rel. corridor], ...] per z level (weather profile)
  function solveUrban(grid, zLevels, p, layers) {
    const zs = zLevels.map((z) => Math.max(z, 10)), out = [];
    zs.forEach((z, k) => {
      const [s, d] = layers ? layers[k] : [layerSpeed(p, z), p.direction_deg], th = d * Math.PI / 180;
      const solid = new Uint8Array(grid.nx * grid.ny); for (let c = 0; c < solid.length; c++) solid[c] = grid.H[c] >= z ? 1 : 0;
      out.push(solveLayer(solid, grid.nx, grid.ny, grid.res, grid.x0, grid.y0, s * Math.cos(th), s * Math.sin(th), p));
    });
    return { z: zs, layers: out, nx: grid.nx, ny: grid.ny, res: grid.res, x0: grid.x0, y0: grid.y0 };
  }
  function sample(W, x, y, z, out) {
    const nx = W.nx, ny = W.ny, nl = W.z.length;
    const fx = Math.min(Math.max((x - W.x0) / W.res - 0.5, 0), nx - 1.000001), fy = Math.min(Math.max((y - W.y0) / W.res - 0.5, 0), ny - 1.000001);
    const i = Math.floor(fx), j = Math.floor(fy), ax = fx - i, ay = fy - j;
    let k = 0; while (k < nl - 2 && W.z[k + 1] < z) k++;
    const az = nl > 1 ? Math.min(Math.max((z - W.z[k]) / (W.z[k + 1] - W.z[k]), 0), 1) : 0;
    const bil = (F) => (1 - ax) * (1 - ay) * F[i * ny + j] + ax * (1 - ay) * F[(i + 1) * ny + j] + (1 - ax) * ay * F[i * ny + j + 1] + ax * ay * F[(i + 1) * ny + j + 1];
    const L0 = W.layers[k], L1 = W.layers[Math.min(k + 1, nl - 1)];
    out = out || [0, 0]; out[0] = (1 - az) * bil(L0.u) + az * bil(L1.u); out[1] = (1 - az) * bil(L0.v) + az * bil(L1.v); return out;
  }

  // ================================================================= TRR link costs (Eqs. 10-33)
  /* veh: {cruise_mph, energy_j_per_ft, phase_dependent, k_climb, alpha_desc, takeoff_fraction}
     model: 'paper' (Eq. 20 piecewise) | 'continuous' (v = v_uav + psi). Returns {t:[B][A] s, e:[B][A] J}. */
  function linkCosts(C, wind, veh, model = "paper") {
    const A = C.A, vU = veh.cruise_mph, fps = vU * 5280 / 3600, B = wind.bins.length, T = [], E = [];
    const t0 = new Float64Array(A), e0 = new Float64Array(A), ec = new Float64Array(A), horiz = new Float64Array(A), theta = new Float64Array(A);
    for (let a = 0; a < A; a++) {
      const p = C.pos[C.af[a]], q = C.pos[C.at[a]], Lft = C.len[a] * FT, dz = q[2] - p[2];
      t0[a] = Lft / fps;                                                   // Eq. 11
      let mult = 1;
      if (veh.phase_dependent) { const sinPhi = C.len[a] > 0 ? dz / C.len[a] : 0;  // Eq. 31
        mult = dz > 0 ? 1 + veh.k_climb * sinPhi : dz < 0 ? veh.alpha_desc : 1; }  // Eqs. 28-30, 32
      ec[a] = veh.energy_j_per_ft * Lft;                                  // Eq. 12
      e0[a] = ec[a] * mult;
      horiz[a] = Math.hypot(q[0] - p[0], q[1] - p[1]);
      theta[a] = ((Math.atan2(q[1] - p[1], q[0] - p[0]) * 180 / Math.PI) % 360 + 360) % 360;   // Eq. 14
    }
    for (let b = 0; b < B; b++) {
      const w = wind.bins[b], tb = new Float64Array(A), eb = new Float64Array(A);
      for (let a = 0; a < A; a++) {
        let vg = vU;
        if (horiz[a] >= 1e-9) {
          const ub = 0.5 * (w.u[C.af[a]] + w.u[C.at[a]]), vb = 0.5 * (w.v[C.af[a]] + w.v[C.at[a]]);   // Eq. 16
          const ws = Math.hypot(ub, vb), al = ((Math.atan2(vb, ub) * 180 / Math.PI) % 360 + 360) % 360;   // Eq. 17
          let g = ((al - theta[a]) % 360 + 360) % 360; if (g > 180) g = 360 - g;                      // Eq. 18
          const psi = ws * Math.cos(g * Math.PI / 180);                                                // Eq. 19
          vg = model === "paper" ? (g <= 45 ? vU + Math.abs(psi) : g >= 135 ? vU - Math.abs(psi) : vU) : vU + psi;   // Eq. 20
        }
        const lam = vg >= 1 ? vU / vg : Infinity;                          // Eq. 25
        tb[a] = t0[a] * lam; eb[a] = e0[a] * lam * lam;                    // Eqs. 26-27 (beta = lambda^2, Eq. 24)
      }
      T.push(tb); E.push(eb);
    }
    return { t: T, e: E, e0, ec, t0 };
  }
  function noWindRouteEnergy(C, e0, o, d) {
    const best = new Float64Array(C.pos.length).fill(Infinity); best[o] = 0; const h = new Heap(); h.push([0, 0, o]);
    while (h.a.length) { const [c, , u] = h.pop(); if (u === d) return c; if (c > best[u]) continue;
      for (const a of C.out[u]) { const v = C.at[a]; if (C.ground[v] && v !== d) continue; const nc = c + e0[a]; if (nc < best[v]) { best[v] = nc; h.push([nc, 0, v]); } } }
    return Infinity;
  }

  // ================================================================= Pareto label-correcting (Algorithms 1-2)
  function pareto(C, costs, wind, veh, o, d, opts = {}) {
    const depart = opts.depart || 0, horizon = opts.horizon || 3600, mask = opts.mask, gt = !!opts.groundTransit, maxL = opts.maxLabels || 400000;
    const B = costs.t.length, binAt = (t) => Math.min(Math.max(Math.floor(t / wind.bin_sec), 0), B - 1);
    const eTO = veh.phase_dependent ? veh.takeoff_fraction * noWindRouteEnergy(C, costs.ec, o, d) : 0;   // Eq. 33, E_TO
    const Lt = [0], Le = [eTO], Ln = [o], Lp = [-1], La = [-1], alive = [true];
    const onPath = (lid, v) => { for (let x = lid; x !== -1; x = Lp[x]) if (Ln[x] === v) return true; return false; };
    const buckets = new Map(); let front = []; const q = [0]; let head = 0;
    while (head < q.length) {
      const lid = q[head++]; if (!alive[lid]) continue;
      const u = Ln[lid], tu = Lt[lid], eu = Le[lid], b = binAt(depart + tu);
      for (const a of C.out[u]) {
        if (mask && !mask[a]) continue;
        const v = C.at[a];
        if (onPath(lid, v)) continue;
        if (!gt && C.ground[v] && v !== d) continue;
        const dt = costs.t[b][a], de = costs.e[b][a]; if (!isFinite(dt)) continue;
        const tv = tu + dt, ev = eu + de; if (tv > horizon) continue;
        let dom = false; for (const f of front) if (Lt[f] <= tv + TOL && Le[f] <= ev + TOL) { dom = true; break; } if (dom) continue;
        if (v === d) { front = front.filter((f) => !(tv <= Lt[f] + TOL && ev <= Le[f] + TOL)); Lt.push(tv); Le.push(ev); Ln.push(v); Lp.push(lid); La.push(a); alive.push(true); front.push(Lt.length - 1); continue; }
        const key = v * 4096 + binAt(depart + tv); let bk = buckets.get(key); if (!bk) { bk = []; buckets.set(key, bk); }
        let d2 = false; for (const x of bk) if (alive[x] && Lt[x] <= tv + TOL && Le[x] <= ev + TOL) { d2 = true; break; } if (d2) continue;
        for (const x of bk) if (alive[x] && tv <= Lt[x] + TOL && ev <= Le[x] + TOL) alive[x] = false;
        const nb = bk.filter((x) => alive[x]);
        Lt.push(tv); Le.push(ev); Ln.push(v); Lp.push(lid); La.push(a); alive.push(true); nb.push(Lt.length - 1); buckets.set(key, nb); q.push(Lt.length - 1);
        if (Lt.length > maxL) return { routes: [], overflow: true };
      }
    }
    front.sort((x, y) => Lt[x] - Lt[y] || Le[x] - Le[y]);
    const routes = front.map((f) => { const nodes = [], arcs = [], T = [], E = []; for (let x = f; x !== -1; x = Lp[x]) { nodes.push(Ln[x]); T.push(Lt[x]); E.push(Le[x]); if (La[x] >= 0) arcs.push(La[x]); }
      nodes.reverse(); arcs.reverse(); T.reverse(); E.reverse(); return { nodes, arcs, T, E, time: Lt[f], energy: Le[f] }; });
    return { routes, labels: Lt.length, eTO };
  }
  function evaluate(C, costs, wind, arcs, depart, eTO = 0) {   // fly a fixed arc list in (possibly another) wind
    const B = costs.t.length, binAt = (t) => Math.min(Math.max(Math.floor(t / wind.bin_sec), 0), B - 1);
    const T = [0], E = [eTO]; for (const a of arcs) { const b = binAt(depart + T[T.length - 1]); T.push(T[T.length - 1] + costs.t[b][a]); E.push(E[E.length - 1] + costs.e[b][a]); }
    return { T, E };
  }
  function routeProfile(C, r, nBlocks = 5) {   // ascents/descents, % time & energy by altitude band and block
    const zf = r.nodes.map((n) => Math.round(C.pos[n][2] * FT)), x = r.nodes.map((n) => C.pos[n][0]);
    const at = {}, ae = {}, bt = new Array(nBlocks).fill(0), be = new Array(nBlocks).fill(0); let asc = 0, des = 0, tt = 0, te = 0;
    for (let k = 0; k < r.nodes.length - 1; k++) {
      const dt = r.T[k + 1] - r.T[k], de = r.E[k + 1] - r.E[k]; tt += dt; te += de;
      if (zf[k + 1] > zf[k]) asc++; if (zf[k + 1] < zf[k]) des++;
      const lab = zf[k] === zf[k + 1] ? `${zf[k]}` : `${Math.min(zf[k], zf[k + 1])}-${Math.max(zf[k], zf[k + 1])}`;
      at[lab] = (at[lab] || 0) + dt; ae[lab] = (ae[lab] || 0) + de;
      const blk = Math.min(Math.floor(0.5 * (x[k] + x[k + 1]) / C.L * nBlocks), nBlocks - 1); bt[blk] += dt; be[blk] += de;
    }
    const pct = (o, s) => Object.fromEntries(Object.entries(o).map(([k, v]) => [k, 100 * v / s]));
    return { ascents: asc, descents: des, altTime: pct(at, tt), altEnergy: pct(ae, te), blockTime: bt.map((v) => 100 * v / tt), blockEnergy: be.map((v) => 100 * v / te) };
  }

  // ================================================================= routing helpers
  function Heap() { this.a = []; }
  Heap.prototype.less = function (x, y) { return x[0] < y[0] || (x[0] === y[0] && (x[1] < y[1] || (x[1] === y[1] && x[2] < y[2]))); };
  Heap.prototype.push = function (x) { const a = this.a; a.push(x); let i = a.length - 1; while (i > 0) { const p = (i - 1) >> 1; if (this.less(a[i], a[p])) { [a[i], a[p]] = [a[p], a[i]]; i = p; } else break; } };
  Heap.prototype.pop = function () { const a = this.a, top = a[0], last = a.pop(); if (a.length) { a[0] = last; let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i; if (l < a.length && this.less(a[l], a[m])) m = l; if (r < a.length && this.less(a[r], a[m])) m = r; if (m === i) break; [a[i], a[m]] = [a[m], a[i]]; i = m; } } return top; };

  // ================================================================= SIRE (TRB paper)
  function rng(seed) { let s = seed >>> 0; return () => { s = s + 0x6D2B79F5 | 0; let t = Math.imul(s ^ s >>> 15, 1 | s); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
  /* Balanced demand (TRB): per-operator demand split evenly; forward share entry->exit, the rest
     reverse; each direction spread evenly over its OD pairs, remainder to the first pairs. */
  function demand(C, D, opts = {}) {
    const ops = opts.operators || 3, win = opts.window || 60, share = opts.forward_share ?? 1, R = rng(opts.seed || 1), pattern = opts.pattern || "uniform";
    const F = [], Rv = []; for (const o of C.origins) for (const d of C.destinations) { F.push([o, d]); Rv.push([d, o]); }
    const trips = []; let id = 0;
    for (let m = 1; m <= ops; m++) {
      const dm = Math.floor(D / ops) + (m - 1 < D % ops ? 1 : 0), nf = Math.round(share * dm), ods = [];
      for (const [pairs, n] of [[F, nf], [Rv, dm - nf]]) { if (!pairs.length) continue; const base = Math.floor(n / pairs.length), rem = n % pairs.length; pairs.forEach((od, i) => { for (let k = 0; k < base + (i < rem ? 1 : 0); k++) ods.push(od); }); }
      const times = ods.map(() => { if (pattern === "peak") { const g = Math.sqrt(-2 * Math.log(R() + 1e-12)) * Math.cos(2 * Math.PI * R()); return Math.min(win - 1, Math.max(0, Math.round(Math.floor(win / 2) + g * Math.max(1, win / 4)))); } return Math.floor(R() * win); });
      for (let i = times.length - 1; i > 0; i--) { const j = Math.floor(R() * (i + 1)); [times[i], times[j]] = [times[j], times[i]]; }
      ods.forEach((od, i) => trips.push({ id: ++id, op: m, o: od[0], d: od[1], dispatch: times[i] }));
    }
    return trips;
  }
  /* tauSec: per-arc traversal seconds; opts: {slot, H, rule: 'undirected'|'directed', groundTransit, tieBreak, seed, onProgress} */
  function sire(C, tauSec, trips, opts) {
    tauSec = Array.from(tauSec, (x) => isFinite(x) ? Math.round(x * 1e6) / 1e6 : Infinity);
    const slot = opts.slot, H = opts.H, rule = opts.rule || "undirected", gt = opts.groundTransit !== false, A = C.A, q = new Int32Array(A);
    for (let a = 0; a < A; a++) q[a] = Math.max(1, Math.ceil(tauSec[a] / slot - 1e-9));
    const nk = rule === "directed" ? A : C.edges.length, key = (a) => rule === "directed" ? a : C.edge[a];
    let qmax = 1; for (let a = 0; a < A; a++) if (isFinite(tauSec[a]) && q[a] > qmax) qmax = q[a];
    const Wd = H + qmax + 2, blocked = new Uint8Array(nk * Wd), N = C.pos.length;
    const R = rng((opts.seed || 0) + 104729), nOp = Math.max(...trips.map((t) => t.op));
    const rk = trips.map((t) => opts.tieBreak === "random" ? R() : opts.tieBreak === "rotating" ? ((t.op - 1 - t.dispatch) % nOp + nOp) % nOp : t.op);
    const order = trips.map((t, i) => i).sort((i, j) => trips[i].dispatch - trips[j].dispatch || rk[i] - rk[j] || trips[i].id - trips[j].id);
    const best = new Float64Array((H + 1) * N), pred = new Int32Array((H + 1) * N), predT = new Int32Array((H + 1) * N), out = new Array(trips.length);
    order.forEach((ti, pos) => {
      const tr = trips[ti], t0 = tr.dispatch, o = tr.o, d = tr.d;
      best.fill(Infinity); best[t0 * N + o] = 0; let bd = Infinity, bt = -1, fe = t0;
      for (let t = t0; t < H && t <= fe; t++) {
        if (best[t * N + d] < bd) { bd = best[t * N + d]; bt = t; }
        for (let u = 0; u < N; u++) { const cu = best[t * N + u]; if (cu === Infinity || cu >= bd || u === d) continue;
          for (const a of C.out[u]) { if (!isFinite(tauSec[a])) continue; const v = C.at[a]; if (!gt && C.ground[v] && v !== d) continue;
            const ta = t + q[a]; if (ta >= H) continue; const k = key(a) * Wd; let free = true; for (let s = t; s < ta; s++) if (blocked[k + s]) { free = false; break; } if (!free) continue;
            const c = cu + tauSec[a]; if (c < best[ta * N + v]) { best[ta * N + v] = c; pred[ta * N + v] = a; predT[ta * N + v] = t; if (ta > fe) fe = ta; } } }
      }
      for (let t = t0; t <= H; t++) if (best[t * N + d] < bd) { bd = best[t * N + d]; bt = t; }
      if (!isFinite(bd)) { out[ti] = { ok: false, order: pos }; } else {
        const legs = []; let n = d, t = bt; while (!(n === o && t === t0)) { const a = pred[t * N + n], s = predT[t * N + n]; legs.push([a, s, t]); n = C.af[a]; t = s; }
        legs.reverse(); for (const [a, s, e] of legs) for (let x = s; x < e; x++) blocked[key(a) * Wd + x] = 1;
        out[ti] = { ok: true, legs, travel: bd, order: pos };
      }
      if (opts.onProgress && pos % 50 === 0) opts.onProgress(pos / trips.length);
    });
    return out;
  }
  function gini(a) { const n = a.length, m = a.reduce((s, x) => s + x, 0) / n; if (!n || m === 0) return 0; let s = 0; for (const x of a) for (const y of a) s += Math.abs(x - y); return s / (2 * n * n * m); }   // Eq. 18
  function kpis(trips, res, W = 5) {
    const D = trips.length, acc = res.filter((r) => r.ok), ops = {};
    trips.forEach((t, i) => { (ops[t.op] = ops[t.op] || [0, 0]); ops[t.op][0]++; if (res[i].ok) ops[t.op][1]++; });
    const opSR = Object.keys(ops).sort((a, b) => a - b).map((k) => ops[k][1] / ops[k][0]);
    const F = new Array(D); res.forEach((r) => { F[r.order] = r.ok ? 0 : 1; });
    let rej = 0, foll = 0; for (let k = 0; k < D; k++) if (F[k]) { rej++; for (let qq = 1; qq <= W && k + qq < D; qq++) if (F[k + qq]) { foll++; break; } }   // Eq. 19
    return { demand: D, accepted: acc.length, sr: acc.length / D, tt: acc.reduce((s, r) => s + r.travel, 0) / Math.max(1, acc.length),
      gini: gini(opSR), cascade: rej ? foll / rej : 0, opSR };
  }

  // ================================================================= OpenStreetMap import (browser only)
  const DEF_H = { house: 7, residential: 9, garage: 4, garages: 4, shed: 3, roof: 5, parking: 12, commercial: 12, retail: 7, industrial: 9, warehouse: 9, office: 20, apartments: 15, church: 15, school: 10, hotel: 30, yes: 10 };
  function parseHeight(tags) {
    const num = (s) => { if (s == null) return null; s = String(s).trim().toLowerCase().replace(",", "."); const f = parseFloat(s); if (!isFinite(f)) return null; return /ft|'/.test(s) ? f * 0.3048 : f; };
    let h = num(tags.height); if (h == null && tags["building:levels"]) { const l = num(tags["building:levels"]); h = l ? l * 3.5 : null; }
    if (h == null) h = DEF_H[tags.building || tags["building:part"] || "yes"] || 10;
    return h;
  }
  function parseOSM(text) {
    const doc = new DOMParser().parseFromString(text, "application/xml");
    const bnd = doc.querySelector("bounds"), nodes = new Map();
    doc.querySelectorAll("node").forEach((n) => nodes.set(n.getAttribute("id"), [+n.getAttribute("lat"), +n.getAttribute("lon")]));
    let lat0, lon0;
    if (bnd) { lat0 = (+bnd.getAttribute("minlat") + +bnd.getAttribute("maxlat")) / 2; lon0 = (+bnd.getAttribute("minlon") + +bnd.getAttribute("maxlon")) / 2; }
    else { let s = [0, 0]; nodes.forEach((v) => { s[0] += v[0]; s[1] += v[1]; }); lat0 = s[0] / nodes.size; lon0 = s[1] / nodes.size; }
    const kx = 111132 * Math.cos(lat0 * Math.PI / 180), ky = 111132, en = ([la, lo]) => [(lo - lon0) * kx, (la - lat0) * ky];
    const buildings = [], roads = [];
    doc.querySelectorAll("way").forEach((w) => {
      const tags = {}; w.querySelectorAll("tag").forEach((t) => { tags[t.getAttribute("k")] = t.getAttribute("v"); });
      const refs = [...w.querySelectorAll("nd")].map((n) => n.getAttribute("ref")).filter((r) => nodes.has(r));
      if ((tags.building || tags["building:part"]) && refs.length >= 4 && refs[0] === refs[refs.length - 1]) {
        const p = []; refs.slice(0, -1).forEach((r) => { const e = en(nodes.get(r)); p.push(e[0], e[1]); });
        buildings.push({ p, h: parseHeight(tags), m: tags.height || tags["building:levels"] ? 1 : 0 });
      } else if (tags.highway && refs.length >= 2) {
        const p = []; refs.forEach((r) => { const e = en(nodes.get(r)); p.push(e[0], e[1]); });
        roads.push({ p, M: /motorway|trunk|primary/.test(tags.highway) ? 1 : 0 });
      }
    });
    return { origin: [lat0, lon0], buildings, roads };
  }

  const api = { FT, MPH_PER_MS, corridorFrame, toCorr, toEN, rasterize, buildCorridor, finish,
    windNone, windShear, windSpatial, windTemporal, windFromUrban, solveLayer, solveUrban, sample, layerSpeed,
    linkCosts, noWindRouteEnergy, pareto, evaluate, routeProfile, Heap, rng, demand, sire, gini, kpis, parseOSM };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.LabEngine = api;
})(typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : globalThis);
