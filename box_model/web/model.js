/* Carbon box model (BEAM ocean + land). Same equations as box_model.py and ocean3.py.
   p.layers = 2: surface + deep ocean (k_d).  p.layers = 3: surface + thermocline + deep (k1m, k2m).
   p.feedback: ocean chemistry re-evaluated at T + dT, dT = TCRE x cumulative emissions / 1000. */
const GTC_PER_PPM = 2.124, V_TH = 12.5, V_LO = 37.5;
function chemistryAt(Tc, S = 35) {
  const T = Tc + 273.15, lg = Math.log10;
  const k0 = Math.exp(9345.17 / T - 60.2409 + 23.3585 * Math.log(T / 100) + S * (0.023517 - 0.00023656 * T + 0.0047036 * (T / 100) ** 2));
  const pK1 = -13.721 + 0.031334 * T + 3235.76 / T + 1.3e-5 * S * T - 0.1031 * S ** 0.5;
  const pK2 = 5371.96 + 1.671221 * T + 0.22913 * S + 18.3802 * lg(S) - 128375.28 / T - 2194.30 * lg(T) - 8.0944e-4 * S * T - 5617.11 * lg(S) / T + 2.136 * S / T;
  return { k1: 10 ** -pK1, k2: 10 ** -pK2, kH: 55.57 / (1.027 * k0) };
}
const CHEM10 = chemistryAt(10);
function makeParams(o = {}) {
  const p = Object.assign({ k_a: 0.2, k_d: 0.05, delta: 50, k1: 8.00e-7, k2: 4.53e-10, kH: 1.23e3, AM: 1.77e20, OM: 7.8e22,
    MUP_ref: 713, pH_ref: 8.29, NPP0: 60, beta: 0.4, tau_L: 35, land: true, layers: 2, k1m: 1.0, k2m: 0.005,
    feedback: false, tcre: 1.8, T: 10 }, o);
  const c = chemistryAt(p.T);                      // BEAM constants scaled from 10 C to T
  p.k1 *= c.k1 / CHEM10.k1; p.k2 *= c.k2 / CHEM10.k2; p.kH *= c.kH / CHEM10.kH;
  p.k1b = p.k1; p.k2b = p.k2; p.kHb = p.kH; p.chemT = c;
  p.A = p.kH * p.AM / (p.OM / (p.delta + 1));
  const h = 10 ** -p.pH_ref;
  p.alk = (p.k1 / h + 2 * p.k1 * p.k2 / h ** 2) * p.MUP_ref * Bfac(h, p);
  return p;
}
function setWarming(p, dT) {                         // feedback: re-evaluate chemistry at T + dT
  const c = chemistryAt(p.T + dT);
  p.k1 = p.k1b * c.k1 / p.chemT.k1; p.k2 = p.k2b * c.k2 / p.chemT.k2; p.kH = p.kHb * c.kH / p.chemT.kH;
  p.A = p.kH * p.AM / (p.OM / (p.delta + 1));
}
function Bfac(h, p) { return 1 / (1 + p.k1 / h + p.k1 * p.k2 / (h * h)); }
function hplus(MUP, p) {
  const a = p.alk, b = (p.alk - MUP) * p.k1, c = (p.alk - 2 * MUP) * p.k1 * p.k2;
  return (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
}
function oceanPco2(MUP, p) { return p.A * Bfac(hplus(MUP, p), p) * MUP; }
function equilibrium(ppm, p) {
  const MAT = ppm * GTC_PER_PPM; let lo = 1, hi = 5000;
  for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2; if (oceanPco2(m, p) < MAT) lo = m; else hi = m; }
  const MUP = (lo + hi) / 2;
  return p.layers === 3 ? [MAT, MUP, V_TH * MUP, V_LO * MUP, p.NPP0 * p.tau_L] : [MAT, MUP, 0, p.delta * MUP, p.NPP0 * p.tau_L];
}
/* state s = [atmosphere, surface, thermocline, deep, land]; returns [F_ocean, F_1, F_2, F_land] */
function fluxes(s, p, MAT0) {
  const Fo = p.k_a * (s[0] - oceanPco2(s[1], p));
  let F1, F2;
  if (p.layers === 3) { F1 = p.k1m * (s[1] - s[2] / V_TH); F2 = p.k2m * (s[2] - s[3] * V_TH / V_LO); }
  else { F1 = p.k_d * (s[1] - s[3] / p.delta); F2 = 0; }
  const Fl = p.land ? p.NPP0 * (1 + p.beta * Math.log(s[0] / MAT0)) - s[4] / p.tau_L : 0;
  return [Fo, F1, F2, Fl];
}
function step(s, p, MAT0, E, dt) {
  const [Fo, F1, F2, Fl] = fluxes(s, p, MAT0);
  const t = p.layers === 3 ? [s[2] + dt * (F1 - F2), s[3] + dt * F2] : [0, s[3] + dt * F1];
  return { s: [s[0] + dt * (E - Fo - Fl), s[1] + dt * (Fo - F1), t[0], t[1], s[4] + dt * Fl], Fo, F1, Fl };
}
/* Annual emissions in GtC/yr. prescribed: optional ppm array (concentration-driven). */
function runModel(Ef, El, p0, startPpm, prescribed = null, dt = 0.01) {
  const p = Object.assign({}, p0);
  let s = equilibrium(startPpm, p); const MAT0 = s[0], n = Math.round(1 / dt), N = Ef.length;
  const keys = ["ppm", "ocean", "land", "mix", "MAT", "MUP", "MTH", "MLO", "ML", "pH", "dT"];
  const o = {}; keys.forEach(k => o[k] = new Float64Array(N));
  let cum = 0;
  for (let i = 0; i < N; i++) {
    let so = 0, sl = 0, sm = 0; const E = Ef[i] + El[i];
    for (let k = 0; k < n; k++) {
      if (prescribed) s[0] = prescribed[i] * GTC_PER_PPM;
      if (p.feedback) setWarming(p, p.tcre * cum / 1000);
      const r = step(s, p, MAT0, E, dt);
      s = r.s; so += r.Fo * dt; sl += r.Fl * dt; sm += r.F1 * dt; cum += E * dt;
    }
    o.MAT[i] = s[0]; o.MUP[i] = s[1]; o.MTH[i] = s[2]; o.MLO[i] = s[3]; o.ML[i] = s[4]; o.ppm[i] = s[0] / GTC_PER_PPM;
    o.ocean[i] = so; o.land[i] = sl; o.mix[i] = sm; o.pH[i] = -Math.log10(hplus(s[1], p)); o.dT[i] = p.tcre * cum / 1000;
  }
  o.state = s; o.MAT0 = MAT0; return o;
}
/* Joos et al. (2013) style: background state, atmosphere held at bg ppm, 100 GtC over the first year (no feedback). */
function pulseTest(p0, bgState, MAT0, bg = 389, pulse = 100, years = 1000, dt = 0.01) {
  const p = Object.assign({}, p0, { feedback: false });
  let c = bgState.slice(); c[0] = bg * GTC_PER_PPM; let q = c.slice(); const n = Math.round(1 / dt);
  const atm = new Float64Array(years), ocean = new Float64Array(years), land = new Float64Array(years);
  for (let y = 0; y < years; y++) {
    for (let k = 0; k < n; k++) {
      const rc = step(c, p, MAT0, 0, dt); const E = rc.Fo + rc.Fl;          // holds the control atmosphere fixed
      c = rc.s; c[0] = bg * GTC_PER_PPM;
      q = step(q, p, MAT0, E + (y === 0 ? pulse : 0), dt).s;
    }
    atm[y] = (q[0] - c[0]) / pulse; ocean[y] = (q[1] + q[2] + q[3] - c[1] - c[2] - c[3]) / pulse; land[y] = (q[4] - c[4]) / pulse;
  }
  return { atm, ocean, land };
}
if (typeof module !== "undefined") module.exports = { makeParams, runModel, pulseTest, equilibrium, chemistryAt, GTC_PER_PPM };
