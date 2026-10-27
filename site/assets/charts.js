// Small dependency-free SVG/canvas charts for the Navette site.
// Every chart reads CSS custom properties, so light and dark themes come free.

const NS = "http://www.w3.org/2000/svg";

function el(tag, attrs = {}, parent) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}

function css(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

// SVG presentation attributes do not take var() everywhere; resolve it here.
function rc(c) {
  if (typeof c !== "string") return c;
  return c.replace(/var\((--[\w-]+)\)/g, (_, v) => css(v));
}

// "Nice" tick values covering [lo, hi].
function ticks(lo, hi, count = 5) {
  const span = hi - lo;
  const step0 = span / count;
  const mag = 10 ** Math.floor(Math.log10(step0));
  const err = step0 / mag;
  const step = (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1) * mag;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) {
    out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  }
  return out;
}

function logTicks(lo, hi) {
  const out = [];
  for (let e = Math.ceil(Math.log10(lo)); e <= Math.floor(Math.log10(hi)); e++) out.push(10 ** e);
  return out;
}

function fmt(v) {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e4 || a < 1e-3) {
    const e = Math.floor(Math.log10(a));
    const m = v / 10 ** e;
    return (Math.abs(m - 1) < 1e-9 ? "" : m.toFixed(0) + "·") + "10" + sup(e);
  }
  return +v.toPrecision(4) + "";
}

function sup(n) {
  const map = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
  return String(n).split("").map((c) => map[c]).join("");
}

// Approximate visible colour for a wavelength in nm (Bruton's piecewise fit).
export function wavelengthColor(wl, alpha = 1) {
  let r = 0, g = 0, b = 0;
  if (wl >= 380 && wl < 440) { r = -(wl - 440) / 60; b = 1; }
  else if (wl < 490) { g = (wl - 440) / 50; b = 1; }
  else if (wl < 510) { g = 1; b = -(wl - 510) / 20; }
  else if (wl < 580) { r = (wl - 510) / 70; g = 1; }
  else if (wl < 645) { r = 1; g = -(wl - 645) / 65; }
  else if (wl <= 780) { r = 1; }
  let f = 1;
  if (wl < 420) f = 0.3 + 0.7 * (wl - 380) / 40;
  else if (wl > 700) f = 0.3 + 0.7 * (780 - wl) / 80;
  if (wl < 380 || wl > 780) f = 0.25, r = g = b = 0.6;
  const c = (x) => Math.round(255 * (x * f) ** 0.8);
  return `rgba(${c(r)},${c(g)},${c(b)},${alpha})`;
}

// Perceptual colour maps for heatmaps (stops from matplotlib).
const MAPS = {
  magma: ["#000004", "#1c1044", "#4f127b", "#812581", "#b5367a", "#e55064", "#fb8761", "#fec287", "#fcfdbf"],
  viridis: ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"],
};

function hexRgb(h) {
  const n = parseInt(h.slice(1), 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}

export function colormap(name, t) {
  const stops = MAPS[name].map(hexRgb);
  const x = Math.min(1, Math.max(0, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(x));
  const f = x - i;
  return stops[i].map((c, k) => Math.round(c + f * (stops[i + 1][k] - c)));
}

/**
 * Line chart.
 * opts: { x: {min,max,label,log}, y: {min,max,label,log,fmt}, series: [{name,x,y,color,width,dash,area,points}],
 *         bands: [{x0,x1,fill,label}], vlines: [{x,label}], height, legend, gradientX: fn(x)->color }
 */
export function lineChart(container, opts) {
  const draw = () => {
    container.innerHTML = "";
    const W = Math.max(280, container.clientWidth);
    const H = opts.height || Math.round(Math.min(380, Math.max(240, W * 0.52)));
    const narrow = W < 520;
    const m = { l: narrow ? 46 : 58, r: 16, t: 14, b: narrow ? 40 : 46 };
    const iw = W - m.l - m.r, ih = H - m.t - m.b;
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, class: "chart", role: "img" });
    if (opts.aria) svg.setAttribute("aria-label", opts.aria);
    container.appendChild(svg);

    const X = opts.x, Y = opts.y;
    const sx = X.log
      ? (v) => m.l + (Math.log10(v) - Math.log10(X.min)) / (Math.log10(X.max) - Math.log10(X.min)) * iw
      : (v) => m.l + (v - X.min) / (X.max - X.min) * iw;
    const sy = Y.log
      ? (v) => m.t + ih - (Math.log10(v) - Math.log10(Y.min)) / (Math.log10(Y.max) - Math.log10(Y.min)) * ih
      : (v) => m.t + ih - (v - Y.min) / (Y.max - Y.min) * ih;

    const defs = el("defs", {}, svg);
    const clipId = "clip" + Math.random().toString(36).slice(2);
    el("rect", { x: m.l, y: m.t, width: iw, height: ih }, el("clipPath", { id: clipId }, defs));

    // bands
    for (const b of opts.bands || []) {
      const x0 = sx(Math.max(b.x0, X.min)), x1 = sx(Math.min(b.x1, X.max));
      el("rect", { x: x0, y: m.t, width: Math.max(0, x1 - x0), height: ih, fill: rc(b.fill || "var(--band)"), class: b.cls || "" }, svg);
      if (b.label) {
        const t = el("text", { x: (x0 + x1) / 2, y: m.t + 16, class: "band-label", "text-anchor": "middle" }, svg);
        t.textContent = b.label;
      }
    }

    // grid + axes
    const g = el("g", { class: "grid" }, svg);
    const xt = X.ticks || (X.log ? logTicks(X.min, X.max) : ticks(X.min, X.max, narrow ? 4 : 7));
    const yt = Y.ticks || (Y.log ? logTicks(Y.min, Y.max) : ticks(Y.min, Y.max, 5));
    for (const v of xt) {
      el("line", { x1: sx(v), x2: sx(v), y1: m.t, y2: m.t + ih }, g);
      const t = el("text", { x: sx(v), y: m.t + ih + 18, "text-anchor": "middle", class: "tick" }, svg);
      t.textContent = (X.fmt || fmt)(v);
    }
    for (const v of yt) {
      el("line", { x1: m.l, x2: m.l + iw, y1: sy(v), y2: sy(v) }, g);
      const t = el("text", { x: m.l - 8, y: sy(v) + 4, "text-anchor": "end", class: "tick" }, svg);
      t.textContent = (Y.fmt || fmt)(v);
    }
    el("rect", { x: m.l, y: m.t, width: iw, height: ih, class: "frame" }, svg);
    if (X.label) {
      const t = el("text", { x: m.l + iw / 2, y: H - 6, "text-anchor": "middle", class: "axis-label" }, svg);
      t.textContent = X.label;
    }
    if (Y.label) {
      const t = el("text", { x: 0, y: 0, "text-anchor": "middle", class: "axis-label", transform: `translate(13 ${m.t + ih / 2}) rotate(-90)` }, svg);
      t.textContent = Y.label;
    }

    // optional spectral strip under the x axis
    if (opts.spectrumStrip) {
      const gid = "spec" + Math.random().toString(36).slice(2);
      const lg = el("linearGradient", { id: gid, x1: 0, x2: 1, y1: 0, y2: 0 }, defs);
      for (let i = 0; i <= 40; i++) {
        const v = X.min + (X.max - X.min) * i / 40;
        el("stop", { offset: i / 40, "stop-color": wavelengthColor(v) }, lg);
      }
      el("rect", { x: m.l, y: m.t + ih - 4, width: iw, height: 4, fill: `url(#${gid})` }, svg);
    }

    const plot = el("g", { "clip-path": `url(#${clipId})` }, svg);
    const paths = [];
    for (const s of opts.series) {
      let d = "", pen = false;
      const pts = [];
      for (let i = 0; i < s.x.length; i++) {
        const xv = s.x[i], yv = s.y[i];
        const ok = yv !== null && Number.isFinite(yv) && (!Y.log || yv > 0) && (!X.log || xv > 0);
        if (!ok) { pen = false; continue; }
        const px = sx(xv), py = sy(Y.log ? yv : yv);
        d += (pen ? "L" : "M") + px.toFixed(2) + " " + py.toFixed(2);
        pen = true;
        pts.push([px, py, xv, yv]);
      }
      if (s.area && pts.length) {
        const base = sy(Y.log ? Y.min : Math.max(Y.min, 0));
        const a = `M${pts[0][0]} ${base}` + pts.map((p) => `L${p[0]} ${p[1]}`).join("") + `L${pts[pts.length - 1][0]} ${base}Z`;
        el("path", { d: a, fill: rc(s.area), stroke: "none" }, plot);
      }
      let stroke = rc(s.color);
      if (s.spectral) {
        const gid = "sg" + Math.random().toString(36).slice(2);
        const lg = el("linearGradient", { id: gid, gradientUnits: "userSpaceOnUse", x1: sx(X.min), x2: sx(X.max), y1: 0, y2: 0 }, defs);
        for (let i = 0; i <= 40; i++) {
          const v = X.min + (X.max - X.min) * i / 40;
          el("stop", { offset: i / 40, "stop-color": wavelengthColor(v) }, lg);
        }
        stroke = `url(#${gid})`;
      }
      const p = el("path", { d, fill: "none", stroke, "stroke-width": s.width ?? 2, "stroke-dasharray": s.dash || "", "stroke-linejoin": "round", "stroke-linecap": "round", class: "series" }, plot);
      if (s.points) for (const q of pts) el("circle", { cx: q[0], cy: q[1], r: s.points, fill: rc(s.color) }, plot);
      paths.push({ s, pts, p });
    }

    for (const v of opts.vlines || []) {
      el("line", { x1: sx(v.x), x2: sx(v.x), y1: m.t, y2: m.t + ih, class: "vline" }, svg);
      if (v.label) {
        const t = el("text", { x: sx(v.x) + 6, y: m.t + 14, class: "vline-label" }, svg);
        t.textContent = v.label;
      }
    }

    // legend
    if (opts.legend !== false && opts.series.some((s) => s.name)) {
      const lg = document.createElement("div");
      lg.className = "legend";
      for (const s of opts.series) {
        if (!s.name) continue;
        const it = document.createElement("span");
        const sw = document.createElement("i");
        sw.style.background = s.spectral ? "linear-gradient(90deg,#7a3cff,#2fb8ff,#3ddc84,#ffd23f,#ff4d4d)" : s.color;
        if (s.dash) sw.className = "dash";
        it.append(sw, document.createTextNode(s.name));
        lg.appendChild(it);
      }
      container.appendChild(lg);
    }

    // hover readout
    const cross = el("line", { y1: m.t, y2: m.t + ih, class: "cross", visibility: "hidden" }, svg);
    const dots = paths.map((q) => el("circle", { r: 3.5, fill: rc(q.s.spectral ? "var(--fg)" : q.s.color), visibility: "hidden", class: "hover-dot" }, svg));
    const tip = document.createElement("div");
    tip.className = "tip";
    container.appendChild(tip);
    const hide = () => { cross.setAttribute("visibility", "hidden"); dots.forEach((d) => d.setAttribute("visibility", "hidden")); tip.style.opacity = 0; };
    svg.addEventListener("pointerleave", hide);
    svg.addEventListener("pointermove", (ev) => {
      const r = svg.getBoundingClientRect();
      const px = (ev.clientX - r.left) * (W / r.width);
      if (px < m.l || px > m.l + iw) return hide();
      cross.setAttribute("x1", px); cross.setAttribute("x2", px); cross.setAttribute("visibility", "visible");
      let html = "", xval = null;
      paths.forEach((q, k) => {
        if (!q.pts.length) { dots[k].setAttribute("visibility", "hidden"); return; }
        let best = q.pts[0];
        for (const p of q.pts) if (Math.abs(p[0] - px) < Math.abs(best[0] - px)) best = p;
        if (Math.abs(best[0] - px) > iw / 25) { dots[k].setAttribute("visibility", "hidden"); return; }
        dots[k].setAttribute("cx", best[0]); dots[k].setAttribute("cy", best[1]); dots[k].setAttribute("visibility", "visible");
        xval = best[2];
        html += `<div><b style="color:${q.s.spectral ? "var(--fg)" : q.s.color}">●</b> ${q.s.name || ""} <span>${(opts.tipY || fmt)(best[3])}</span></div>`;
      });
      if (xval === null) { tip.style.opacity = 0; return; }
      tip.innerHTML = `<div class="tip-x">${(opts.tipX || fmt)(xval)}</div>` + html;
      tip.style.opacity = 1;
      const cx = (px / W) * r.width;
      tip.style.left = Math.min(r.width - tip.offsetWidth - 4, Math.max(4, cx + 12)) + "px";
      tip.style.top = "8px";
    });
  };
  draw();
  observe(container, draw);
  return { redraw: draw, set(o) { Object.assign(opts, o); draw(); } };
}

/** Heatmap on canvas: z[row=y][col=x]. */
export function heatmap(container, { x, y, z, zmin = 0, zmax = 1, map = "magma", xLabel, yLabel, zLabel, tipFmt }) {
  const draw = () => {
    container.innerHTML = "";
    const W = Math.max(280, container.clientWidth);
    const narrow = W < 520;
    const H = Math.round(Math.min(400, Math.max(240, W * 0.5)));
    const m = { l: narrow ? 46 : 58, r: narrow ? 44 : 64, t: 14, b: narrow ? 40 : 46 };
    const iw = W - m.l - m.r, ih = H - m.t - m.b;
    const wrap = document.createElement("div");
    wrap.className = "heat-wrap";
    wrap.style.height = H + "px";
    container.appendChild(wrap);

    const cv = document.createElement("canvas");
    const dpr = window.devicePixelRatio || 1;
    cv.width = Math.round(iw * dpr); cv.height = Math.round(ih * dpr);
    Object.assign(cv.style, { position: "absolute", left: m.l + "px", top: m.t + "px", width: iw + "px", height: ih + "px", borderRadius: "4px" });
    wrap.appendChild(cv);
    const ctx = cv.getContext("2d");
    const img = ctx.createImageData(cv.width, cv.height);
    const nx = x.length, ny = y.length;
    for (let py = 0; py < cv.height; py++) {
      const fy = (1 - py / (cv.height - 1)) * (ny - 1);
      const j = Math.min(ny - 2, Math.floor(fy)), ty = fy - j;
      for (let px = 0; px < cv.width; px++) {
        const fx = (px / (cv.width - 1)) * (nx - 1);
        const i = Math.min(nx - 2, Math.floor(fx)), tx = fx - i;
        const v = (z[j][i] * (1 - tx) + z[j][i + 1] * tx) * (1 - ty) + (z[j + 1][i] * (1 - tx) + z[j + 1][i + 1] * tx) * ty;
        const c = colormap(map, (v - zmin) / (zmax - zmin));
        const o = (py * cv.width + px) * 4;
        img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);

    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, class: "chart overlay" });
    wrap.appendChild(svg);
    const sx = (v) => m.l + (v - x[0]) / (x[nx - 1] - x[0]) * iw;
    const sy = (v) => m.t + ih - (v - y[0]) / (y[ny - 1] - y[0]) * ih;
    for (const v of ticks(x[0], x[nx - 1], narrow ? 4 : 7)) {
      const t = el("text", { x: sx(v), y: m.t + ih + 18, "text-anchor": "middle", class: "tick" }, svg);
      t.textContent = fmt(v);
    }
    for (const v of ticks(y[0], y[ny - 1], 5)) {
      const t = el("text", { x: m.l - 8, y: sy(v) + 4, "text-anchor": "end", class: "tick" }, svg);
      t.textContent = fmt(v);
    }
    const xl = el("text", { x: m.l + iw / 2, y: H - 6, "text-anchor": "middle", class: "axis-label" }, svg);
    xl.textContent = xLabel;
    const yl = el("text", { "text-anchor": "middle", class: "axis-label", transform: `translate(13 ${m.t + ih / 2}) rotate(-90)` }, svg);
    yl.textContent = yLabel;

    // colour bar
    const defs = el("defs", {}, svg);
    const gid = "cb" + Math.random().toString(36).slice(2);
    const lg = el("linearGradient", { id: gid, x1: 0, x2: 0, y1: 1, y2: 0 }, defs);
    for (let i = 0; i <= 20; i++) el("stop", { offset: i / 20, "stop-color": `rgb(${colormap(map, i / 20).join(",")})` }, lg);
    const bx = m.l + iw + (narrow ? 10 : 14);
    el("rect", { x: bx, y: m.t, width: 10, height: ih, fill: `url(#${gid})`, rx: 2 }, svg);
    for (const v of ticks(zmin, zmax, 4)) {
      const t = el("text", { x: bx + 14, y: m.t + ih - (v - zmin) / (zmax - zmin) * ih + 4, class: "tick" }, svg);
      t.textContent = fmt(v);
    }
    if (zLabel && !narrow) {
      const t = el("text", { "text-anchor": "middle", class: "axis-label", transform: `translate(${W - 4} ${m.t + ih / 2}) rotate(90)` }, svg);
      t.textContent = zLabel;
    }

    const tip = document.createElement("div");
    tip.className = "tip";
    wrap.appendChild(tip);
    cv.addEventListener("pointerleave", () => (tip.style.opacity = 0));
    cv.addEventListener("pointermove", (ev) => {
      const r = cv.getBoundingClientRect();
      const fx = (ev.clientX - r.left) / r.width * (nx - 1);
      const fy = (1 - (ev.clientY - r.top) / r.height) * (ny - 1);
      const i = Math.round(fx), j = Math.round(fy);
      if (i < 0 || j < 0 || i >= nx || j >= ny) return;
      tip.innerHTML = tipFmt(x[i], y[j], z[j][i]);
      tip.style.opacity = 1;
      tip.style.left = Math.min(W - tip.offsetWidth - 4, m.l + (ev.clientX - r.left) + 12) + "px";
      tip.style.top = Math.max(4, m.t + (ev.clientY - r.top) - 40) + "px";
    });
  };
  draw();
  observe(container, draw);
  return { redraw: draw };
}

function observe(container, draw) {
  let w = container.clientWidth, raf = 0;
  new ResizeObserver(() => {
    if (container.clientWidth === w) return;
    w = container.clientWidth;
    cancelAnimationFrame(raf);
    raf = requestAnimationFrame(draw);
  }).observe(container);
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", draw);
  document.addEventListener("themechange", draw);
}

export { css, el, rc };
