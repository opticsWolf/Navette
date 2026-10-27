import { heatmap, lineChart } from "./charts.js";

const $ = (s) => document.querySelector(s);
const pct = (v) => (v * 100).toFixed(v < 0.01 ? 3 : 2) + " %";
const nm = (v) => (+v).toFixed(v % 1 ? 1 : 0) + " nm";
const deg = (v) => (+v).toFixed(1) + "°";

// Theme toggle: remembers the choice; the OS preference is the default.
$("#theme-toggle").addEventListener("click", () => {
  const root = document.documentElement;
  const dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
  root.dataset.theme = dark ? "light" : "dark";
  try { localStorage.setItem("navette-theme", root.dataset.theme); } catch (e) {}
  document.dispatchEvent(new Event("themechange"));
});

// Copy buttons.
async function copy(btn, text) {
  try { await navigator.clipboard.writeText(text); } catch (e) { return; }
  btn.classList.add("done");
  setTimeout(() => btn.classList.remove("done"), 1200);
}
document.querySelectorAll("[data-copy]").forEach((b) => b.addEventListener("click", () => copy(b, b.dataset.copy)));
$("#copy-snippet").addEventListener("click", (e) => {
  copy(e.currentTarget, $("#snippet").innerText);
  e.currentTarget.textContent = "Copied";
  setTimeout(() => (e.currentTarget.textContent = "Copy"), 1200);
});

// A segmented control: one button "on", calls back with its data value.
function seg(selector, key, onChange) {
  const group = document.querySelectorAll(selector);
  group.forEach((b) => b.addEventListener("click", () => {
    group.forEach((o) => o.classList.toggle("on", o === b));
    onChange(b.dataset[key]);
  }));
}

const data = await fetch("data/figures.json").then((r) => r.json());
const v = data.meta.navette;
$("#version-pill").textContent = "v" + v;
$("#data-version").textContent = v;

// --- hero: the mirror, with an angle slider ----------------------------------
{
  const B = data.bragg;
  let pol = "u", ai = 0;
  const curve = () => B.wl.map((_, i) => pol === "s" ? B.Rs[ai][i] : pol === "p" ? B.Rp[ai][i] : 0.5 * (B.Rs[ai][i] + B.Rp[ai][i]));
  const opts = {
    x: { min: 380, max: 780, label: "Wavelength (nm)" },
    y: { min: 0, max: 1, label: "Reflectance", fmt: (t) => t.toFixed(1) },
    series: [{ name: "", x: B.wl, y: curve(), spectral: true, width: 2.6, area: "var(--hi-film)" }],
    legend: false, spectrumStrip: true, tipX: nm, tipY: pct,
    aria: "Reflectance of a quarter-wave mirror versus wavelength",
  };
  const chart = lineChart($("#hero-chart"), opts);
  const update = () => { opts.series[0].y = curve(); chart.redraw(); };
  $("#angle").addEventListener("input", (e) => {
    ai = +e.target.value;
    $("#angle-out").textContent = deg(B.angles[ai]);
    update();
  });
  seg(".hero-fig [data-pol]", "pol", (p) => { pol = p; update(); });
}

// --- stability: silicon wafer --------------------------------------------------
{
  const S = data.stability;
  $("#wafer-first").textContent = Math.round(S.first_finite_nm) + " nm";
  $("#wafer-d").textContent = Math.round(S.d_um) + " µm";
  lineChart($("#stability-chart"), {
    x: { min: 300, max: 1200, label: "Wavelength (nm)" },
    y: { min: 0, max: 0.8, label: "Reflectance", fmt: (t) => t.toFixed(1) },
    bands: [
      { x0: 300, x1: S.first_finite_nm, label: "naive TMM: NaN" },
      { x0: 1000, x1: 1200, label: "aliased", fill: "var(--hi-film)" },
    ],
    // Naive first, so where the two agree the Navette line stays on top.
    series: [
      { name: "Naive transfer matrix", x: S.wl, y: S.naive, color: "var(--s2)", width: 0, points: 1.7 },
      { name: "Front-surface Fresnel", x: S.wl, y: S.fresnel, color: "var(--s5)", width: 1.4, dash: "5 4" },
      { name: "Navette (S-matrix, incoherent wafer)", x: S.wl, y: S.R, color: "var(--s1)", width: 2.4 },
    ],
    tipX: nm, tipY: pct,
  });
}

// --- R(λ, θ) map ---------------------------------------------------------------
{
  const B = data.bragg;
  const draw = (key) => heatmap($("#bragg-map"), {
    x: B.wl, y: B.angles, z: B[key], zmin: 0, zmax: 1, map: "magma",
    xLabel: "Wavelength (nm)", yLabel: "Angle of incidence (°)", zLabel: key === "Rs" ? "Rs" : "Rp",
    tipFmt: (x, y, z) => `<div class="tip-x">${x.toFixed(0)} nm · ${y.toFixed(1)}°</div><div>${key} <span>${pct(z)}</span></div>`,
  });
  draw("Rs");
  seg("[data-map]", "map", draw);
}

// --- coherence -----------------------------------------------------------------
{
  const C = data.coherence;
  $("#fringe-period").textContent = C.period_nm.toFixed(2) + " nm";
  lineChart($("#coherence-chart"), {
    x: { min: 600, max: 601, label: "Wavelength (nm)", fmt: (t) => t.toFixed(1) },
    y: { min: 0.1, max: 0.5, label: "Reflectance", fmt: (t) => t.toFixed(1) },
    series: [
      { name: "Fully coherent", x: C.wl, y: C.coherent, color: "var(--s5)", width: 0.9 },
      { name: "Slide flagged incoherent", x: C.wl, y: C.incoherent, color: "var(--s1)", width: 2.6 },
    ],
    tipX: (t) => t.toFixed(3) + " nm", tipY: pct,
  });
}

// --- group delay and GDD -------------------------------------------------------
{
  const G = data.group_delay;
  const opts = { x: { min: 470, max: 690, label: "Wavelength (nm)" }, y: {}, series: [], vlines: [{ x: G.lam0, label: "λ₀" }], tipX: nm };
  const set = (k) => {
    const gd = k === "GD";
    opts.y = gd ? { min: 0, max: 12, label: "Group delay (fs)" } : { min: -60, max: 60, label: "GDD (fs²)" };
    opts.tipY = gd ? (t) => t.toFixed(2) + " fs" : (t) => t.toFixed(1) + " fs²";
    opts.series = [{ name: k, x: G.wl, y: G[k], color: gd ? "var(--s3)" : "var(--s4)", width: 2.2 }];
    opts.legend = false;
  };
  set("GD");
  const chart = lineChart($("#gd-chart"), opts);
  seg("[data-gd]", "gd", (k) => { set(k); chart.redraw(); });
}

// --- SPR -----------------------------------------------------------------------
{
  const P = data.spr;
  const colors = ["var(--s5)", "var(--s3)", "var(--s2)", "var(--s4)"];
  lineChart($("#spr-chart"), {
    x: { min: P.angles[0], max: P.angles[P.angles.length - 1], label: "Angle in the prism (°)" },
    y: { min: 0, max: 1, label: "Rp", fmt: (t) => t.toFixed(1) },
    series: Object.entries(P.curves).map(([d, y], i) => ({ name: `${d} nm gold`, x: P.angles, y, color: colors[i], width: d === "50" ? 2.6 : 1.6 })),
    tipX: (t) => t.toFixed(2) + "°", tipY: pct,
  });
}

// --- ellipsometry --------------------------------------------------------------
{
  const E = data.ellipsometry;
  const opts = {
    x: { min: 0, max: 90, label: "Ψ (°)" },
    y: { min: 60, max: 300, label: "Δ (°)" },
    series: [
      { name: "trajectory", x: E.psi, y: E.delta, color: "var(--s3)", width: 2 },
      { name: "", x: [], y: [], color: "var(--s2)", width: 0, points: 6 },
    ],
    legend: false, tipX: deg, tipY: deg,
  };
  const chart = lineChart($("#ellips-chart"), opts);
  const set = (i) => {
    opts.series[1].x = [E.psi[i]];
    opts.series[1].y = [E.delta[i]];
    opts.series[1].name = `${E.d[i]} nm`;
    $("#oxide-out").textContent = nm(E.d[i]);
    chart.redraw();
  };
  $("#oxide").max = E.d.length - 1;
  $("#oxide").addEventListener("input", (e) => set(+e.target.value));
  set(+$("#oxide").value);
}

// --- dispersion ----------------------------------------------------------------
{
  const D = data.dispersion;
  const colors = ["var(--s1)", "var(--s2)", "var(--s3)"];
  const opts = { x: { min: 300, max: 1000, label: "Wavelength (nm)" }, y: {}, series: [], tipX: nm, tipY: (t) => t.toFixed(3) };
  const set = (c) => {
    opts.y = c === "n" ? { min: 0, max: 5, label: "n" } : { min: 0, max: 7, label: "k" };
    opts.series = Object.entries(D.materials).map(([name, m], i) => ({ name, x: D.wl, y: m[c], color: colors[i], width: 2 }));
  };
  set("n");
  const chart = lineChart($("#disp-chart"), opts);
  seg("[data-nk]", "nk", (c) => { set(c); chart.redraw(); });
}

// --- needle synthesis ----------------------------------------------------------
{
  const N = data.needle;
  // Log axis: the design sits two decades under bare glass. Zeros of R are
  // floored to the axis minimum rather than dropped.
  const toPct = (a) => a.map((x) => Math.max(1e-3, x * 100));
  lineChart($("#needle-spectrum"), {
    x: { min: 380, max: 780, label: "Wavelength (nm)" },
    y: { min: 1e-3, max: 50, log: true, label: "Reflectance (%)", fmt: (t) => (t >= 1 ? t.toFixed(0) : +t.toPrecision(1) + "") },
    bands: [{ x0: N.band[0], x1: N.band[1], fill: "var(--hi-film)" }],
    series: [
      { name: "Bare glass, n = 1.52", x: N.wl, y: toPct(N.R_bare), color: "var(--s5)", width: 1.6, dash: "6 4" },
      { name: `Needle design, ${N.layers.length} films`, x: N.wl, y: toPct(N.R), spectral: true, width: 2.8 },
    ],
    spectrumStrip: true, tipX: nm, tipY: (t) => t.toFixed(3) + " %",
  });
  // Late needle passes can end worse than they began (an insertion the
  // re-optimization cannot use), and cleanup's merge usually takes it back.
  // The run ends with one more optimization of every thickness; that final
  // design sits one step past the last cycle, so the chart ends where the
  // spectrum and stack come from.
  const H = N.history;
  const last = H[H.length - 1].cycle, fin = last + 1;
  const cyc = H.map((h) => h.cycle);
  lineChart($("#needle-merit"), {
    x: {
      min: 0, max: fin, label: "Macro cycle",
      ticks: [...cyc.filter((c) => c % 2 === 0 && c < last), last, fin],
      fmt: (c) => (c === fin ? "final" : String(c)),
    },
    y: { min: 1, max: 1e6, log: true, label: "Merit" },
    series: [
      { name: "after needle pass", x: cyc, y: H.map((h) => h.needle), color: "var(--s5)", width: 1.5, dash: "4 3" },
      { name: "after cleanup", x: cyc, y: H.map((h) => h.merit), color: "var(--s1)", width: 2.2, points: 3.2 },
      { name: "final design", x: [last, fin], y: [H[H.length - 1].merit, N.final.merit], color: "var(--s2)", width: 1.5, dash: "2 3", points: 4.5 },
    ],
    tipX: (c) => (c === fin ? `final optimization · ${N.final.layers} films` : `cycle ${c} · ${H[c].layers} films`),
    tipY: (t) => t.toPrecision(4),
  });
  $("#needle-rms").textContent = (N.rms_in_band * 100).toFixed(2);
  const box = $("#needle-stack");
  const total = N.layers.reduce((s, l) => s + l.d, 0);
  for (const l of N.layers) {
    const div = document.createElement("div");
    div.className = l.n > 2 ? "hi" : "lo";
    div.style.flex = `${l.d / total} 0 0`;
    div.title = `${l.n > 2 ? "TiO₂" : "SiO₂"} · ${l.d.toFixed(1)} nm`;
    box.appendChild(div);
  }
  const sub = document.createElement("div");
  sub.className = "sub";
  sub.title = "glass substrate";
  box.appendChild(sub);
  $("#stack-sub").textContent = `${N.layers.length} films · ${total.toFixed(0)} nm · ambient on the left`;
}

// --- rugate --------------------------------------------------------------------
{
  const G = data.rugate;
  $("#rugate-rows").textContent = G.rows_graded;
  $("#rugate-rows-s").textContent = G.rows_stepped;
  lineChart($("#rugate-chart"), {
    x: { min: 450, max: 780, label: "Wavelength (nm)" },
    y: { min: 0, max: 1, label: "Reflectance", fmt: (t) => t.toFixed(1) },
    series: [
      { name: `Stepped, ${G.rows_stepped} films`, x: G.wl, y: G.stepped, color: "var(--s2)", width: 1.6 },
      { name: `Graded, ${G.rows_graded} rows`, x: G.wl, y: G.graded, color: "var(--s1)", width: 2.6 },
    ],
    spectrumStrip: true, tipX: nm, tipY: pct,
  });
}
