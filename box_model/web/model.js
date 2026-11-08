/* Four-box carbon model (BEAM ocean + land box). Same equations as box_model.py. */
const GTC_PER_PPM = 2.124;
function chemistryAt(Tc, S = 35) {
  const T = Tc + 273.15, lg = Math.log10;
  const k0 = Math.exp(9345.17 / T - 60.2409 + 23.3585 * Math.log(T / 100) + S * (0.023517 - 0.00023656 * T + 0.0047036 * (T / 100) ** 2));
  const pK1 = -13.721 + 0.031334 * T + 3235.76 / T + 1.3e-5 * S * T - 0.1031 * S ** 0.5;
  const pK2 = 5371.96 + 1.671221 * T + 0.22913 * S + 18.3802 * lg(S) - 128375.28 / T - 2194.30 * lg(T) - 8.0944e-4 * S * T - 5617.11 * lg(S) / T + 2.136 * S / T;
  return { k1: 10 ** -pK1, k2: 10 ** -pK2, kH: 55.57 / (1.027 * k0) };
}
function makeParams(o = {}) {
  const p = Object.assign({ k_a: 0.2, k_d: 0.05, delta: 50, k1: 8.00e-7, k2: 4.53e-10, kH: 1.23e3, AM: 1.77e20, OM: 7.8e22,
    MUP_ref: 713, pH_ref: 8.29, NPP0: 60, beta: 0.4, tau_L: 35, land: true }, o);
  if (o.T != null) Object.assign(p, chemistryAt(o.T));
  p.A = p.kH * p.AM / (p.OM / (p.delta + 1));
  const h = 10 ** -p.pH_ref;
  p.alk = (p.k1 / h + 2 * p.k1 * p.k2 / h ** 2) * p.MUP_ref * Bfac(h, p);
  return p;
}
function Bfac(h, p) { return 1 / (1 + p.k1 / h + p.k1 * p.k2 / (h * h)); }
function hplus(MUP, p) {
  const a = p.alk, b = (p.alk - MUP) * p.k1, c = (p.alk - 2 * MUP) * p.k1 * p.k2;
  return (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
}
function oceanPco2(MUP, p) { return p.A * Bfac(hplus(MUP, p), p) * MUP; }
function revelle(MUP, p) { const e = 1e-4, u = MUP * (1 + e), d = MUP * (1 - e); return (Math.log(oceanPco2(u, p)) - Math.log(oceanPco2(d, p))) / (Math.log(u) - Math.log(d)); }
function equilibrium(ppm, p) {
  const MAT = ppm * GTC_PER_PPM; let lo = 1, hi = 5000;
  for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2; if (oceanPco2(m, p) < MAT) lo = m; else hi = m; }
  const MUP = (lo + hi) / 2;
  return [MAT, MUP, p.delta * MUP, p.NPP0 * p.tau_L];
}
function fluxes(s, p, MAT0) {
  const Fo = p.k_a * (s[0] - oceanPco2(s[1], p));
  const Fm = p.k_d * (s[1] - s[2] / p.delta);
  const Fl = p.land ? p.NPP0 * (1 + p.beta * Math.log(s[0] / MAT0)) - s[3] / p.tau_L : 0;
  return [Fo, Fm, Fl];
}
/* Annual emissions in GtC/yr. prescribed: optional ppm array (concentration-driven). */
function runModel(Ef, El, p, startPpm, prescribed = null, dt = 0.01) {
  let s = equilibrium(startPpm, p); const MAT0 = s[0], n = Math.round(1 / dt), N = Ef.length;
  const o = { ppm: new Float64Array(N), ocean: new Float64Array(N), land: new Float64Array(N), mix: new Float64Array(N),
    MAT: new Float64Array(N), MUP: new Float64Array(N), MLO: new Float64Array(N), ML: new Float64Array(N), pH: new Float64Array(N) };
  for (let i = 0; i < N; i++) {
    let so = 0, sl = 0, sm = 0;
    for (let k = 0; k < n; k++) {
      if (prescribed) s[0] = prescribed[i] * GTC_PER_PPM;
      const [Fo, Fm, Fl] = fluxes(s, p, MAT0);
      s = [s[0] + dt * (Ef[i] + El[i] - Fo - Fl), s[1] + dt * (Fo - Fm), s[2] + dt * Fm, s[3] + dt * Fl];
      so += Fo * dt; sl += Fl * dt; sm += Fm * dt;
    }
    o.MAT[i] = s[0]; o.MUP[i] = s[1]; o.MLO[i] = s[2]; o.ML[i] = s[3]; o.ppm[i] = s[0] / GTC_PER_PPM;
    o.ocean[i] = so; o.land[i] = sl; o.mix[i] = sm; o.pH[i] = -Math.log10(hplus(s[1], p));
  }
  o.state = s; o.MAT0 = MAT0; return o;
}
/* Joos et al. (2013) style: background state, atmosphere held at bg ppm, 100 GtC over the first year. */
function pulseTest(p, bgState, MAT0, bg = 389, pulse = 100, years = 1000, dt = 0.01) {
  let c = bgState.slice(); c[0] = bg * GTC_PER_PPM; let q = c.slice(); const n = Math.round(1 / dt);
  const atm = new Float64Array(years), ocean = new Float64Array(years), land = new Float64Array(years);
  for (let y = 0; y < years; y++) {
    for (let k = 0; k < n; k++) {
      const [Fo, Fm, Fl] = fluxes(c, p, MAT0);
      c = [c[0], c[1] + dt * (Fo - Fm), c[2] + dt * Fm, c[3] + dt * Fl];
      const Ep = Fo + Fl + (y === 0 ? pulse : 0);
      const [Fo2, Fm2, Fl2] = fluxes(q, p, MAT0);
      q = [q[0] + dt * (Ep - Fo2 - Fl2), q[1] + dt * (Fo2 - Fm2), q[2] + dt * Fm2, q[3] + dt * Fl2];
    }
    atm[y] = (q[0] - c[0]) / pulse; ocean[y] = (q[1] + q[2] - c[1] - c[2]) / pulse; land[y] = (q[3] - c[3]) / pulse;
  }
  return { atm, ocean, land };
}
if (typeof module !== "undefined") module.exports = { makeParams, runModel, pulseTest, equilibrium, revelle, chemistryAt, GTC_PER_PPM };
