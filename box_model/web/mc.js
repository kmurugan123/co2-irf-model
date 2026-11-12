/* Monte Carlo experiment engine for the box model (runs inside a Web Worker; also testable in Node).
   Needs model.js (makeParams, runModel). Same method as experiments3.monte_carlo in Python:
   draw the uncertain inputs, turn the drawn ocean/land data into model rates with calibration grids
   (built from concentration-driven runs, never the CO2 record), then run the emission-driven model. */
function mulberry32(seed) { let a = seed >>> 0; return () => { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
function makeRng(seed) {
  const u = mulberry32(seed); let spare = null;
  return { uniform: (a, b) => a + (b - a) * u(),
    normal: (m, s) => { if (spare !== null) { const z = spare; spare = null; return m + s * z; }
      let x, y, r; do { x = 2 * u() - 1; y = 2 * u() - 1; r = x * x + y * y; } while (r >= 1 || r === 0);
      const f = Math.sqrt(-2 * Math.log(r) / r); spare = y * f; return m + s * x * f; } };
}
function interp(x, xs, ys) {                       // linear, clamped; xs increasing
  if (x <= xs[0]) return ys[0]; const n = xs.length; if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0, hi = n - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  const w = (x - xs[lo]) / (xs[hi] - xs[lo]); return ys[lo] + w * (ys[hi] - ys[lo]);
}
function k1For(inventory, k2, CAL) {
  const lk2 = CAL.K2.map(Math.log), lk1 = CAL.K1.map(Math.log);
  const col = CAL.inv3.map(row => interp(Math.log(k2), lk2, row));      // inventory vs k1 at this k2
  return { k1: Math.exp(interp(inventory, col, lk1)), atLimit: inventory >= col[col.length - 1] };
}
function kdFor(inventory, CAL) { return Math.exp(interp(inventory, CAL.inv2, CAL.KD.map(Math.log))); }
function betaFor(land, CAL) { return interp(land, CAL.land, CAL.B); }
function scenarioFactor(sc, t, custom = 0) {
  return sc === "flat" ? 1 : sc === "cut2" ? 0.98 ** t : sc === "netzero2050" ? Math.max(0, 1 - t / 26) : sc === "stop" ? 0 : (1 + custom) ** t;
}
/* cfg: { n, seed, layers, feedback, scenario, endYear, k_a, T,
          luc:{on,sd}, fossil:{on,sd}, inv:{on,mean,sd}, age:{on,min,max}, land:{on,mean,sd}, tcre:{on,min,max} } */
function runMC(cfg, D, CAL, onProgress) {
  const rng = makeRng(cfg.seed || 1), n = cfg.n;
  const end = cfg.scenario && cfg.scenario !== "none" ? (cfg.endYear || 2100) : 2024;
  const N = end - 1850 + 1, hist = D.ice.concat(D.mlo.slice(0, 2024 - D.mloY0 + 1));
  const ppm = new Float32Array(n * N), dT = new Float32Array(n);
  const inp = { f_luc: [], f_fossil: [], inventory: [], deep_age: [], land_target: [], tcre: [] };
  const rates = { k1: [], k2: [], k_d: [], beta: [], atLimit: [] };
  const peakYear = new Float32Array(n), peakVal = new Float32Array(n);
  const share = { air: new Float32Array(n), ocean: new Float32Array(n), land: new Float32Array(n) };   // of 1850-2024 emissions
  for (let j = 0; j < n; j++) {
    const fl = cfg.luc.on ? Math.max(0, rng.normal(1, cfg.luc.sd)) : 1;
    const ff = cfg.fossil.on ? rng.normal(1, cfg.fossil.sd) : 1;
    const inv = cfg.inv.on ? rng.normal(cfg.inv.mean, cfg.inv.sd) : cfg.inv.mean;
    const age = cfg.age.on ? rng.uniform(cfg.age.min, cfg.age.max) : 0.5 * (cfg.age.min + cfg.age.max);
    const land = cfg.land.on ? rng.normal(cfg.land.mean, cfg.land.sd) : cfg.land.mean;
    const tcre = cfg.tcre.on ? rng.uniform(cfg.tcre.min, cfg.tcre.max) : 1.8;
    const beta = Math.min(Math.max(betaFor(land, CAL), CAL.B[0]), CAL.B[CAL.B.length - 1]);
    let pr;
    if (cfg.layers === 3) { const k2 = 3 / age, c = k1For(inv, k2, CAL); pr = { layers: 3, k1m: c.k1, k2m: k2 }; rates.k1.push(c.k1); rates.k2.push(k2); rates.atLimit.push(c.atLimit); }
    else { const kd = kdFor(inv, CAL); pr = { layers: 2, k_d: kd }; rates.k_d.push(kd); rates.atLimit.push(inv >= CAL.inv2[CAL.inv2.length - 1]); }
    const p = makeParams(Object.assign(pr, { beta, feedback: !!cfg.feedback, tcre, k_a: cfg.k_a ?? 0.2, T: cfg.T ?? 10 }));
    const Ef = new Float64Array(N), El = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      const k = i < 175 ? 1 : scenarioFactor(cfg.scenario, i - 174, cfg.custom || 0);
      const base = i < 175 ? i : 174;
      Ef[i] = ff * D.fossil[base] * k; El[i] = fl * D.luc[base] * k;
    }
    const r = runModel(Ef, El, p, hist[0], null, cfg.dt || 0.02);
    let pk = 0;
    for (let i = 0; i < N; i++) { ppm[j * N + i] = r.ppm[i]; if (r.ppm[i] > r.ppm[pk]) pk = i; }
    peakYear[j] = 1850 + pk; peakVal[j] = r.ppm[pk]; dT[j] = r.dT[N - 1];
    let eT = 0, oT = 0, lT = 0; for (let i = 0; i < Math.min(175, N); i++) { eT += Ef[i] + El[i]; oT += r.ocean[i]; lT += r.land[i]; }
    share.ocean[j] = oT / eT; share.land[j] = lT / eT; share.air[j] = 1 - share.ocean[j] - share.land[j];
    inp.f_luc.push(fl); inp.f_fossil.push(ff); inp.inventory.push(inv); inp.deep_age.push(age); inp.land_target.push(land); inp.tcre.push(tcre); rates.beta.push(beta);
    if (onProgress && (j % 25 === 24 || j === n - 1)) onProgress(j + 1);
  }
  return { n, N, end, ppm, peakYear, peakVal, dT, inp, rates, share };
}
if (typeof module !== "undefined") module.exports = { runMC, makeRng, interp, k1For, kdFor, betaFor };
