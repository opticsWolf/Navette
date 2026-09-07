# -*- coding: utf-8 -*-
"""Color demands end to end: bound ``compile_merit_spec`` -> merit/residuals
vs independent oracles (Rust ``color_merit`` + ``merit`` twins, plan R3)."""
import json

import numpy as np
import pytest

from navette._smatrix import compile_merit_spec
from navette.synthesis import sim_curves_from_arrays
import navette._color as C
from navette.data import load_cie_table

CMF = load_cie_table("CIE", "cmf", "CIE_xyz_1931_2deg.json")
D65 = load_cie_table("CIE", "sds", "CIE_std_illum_D65_S_D65.json")
LAM = np.asarray(CMF["lambda"], dtype=np.float64)
CX = np.asarray(CMF["x_bar(lambda)"], dtype=np.float64)
CY = np.asarray(CMF["y_bar(lambda)"], dtype=np.float64)
CZ = np.asarray(CMF["z_bar(lambda)"], dtype=np.float64)
DLAM = np.asarray(D65["lambda"], dtype=np.float64)
E65 = np.asarray(D65["S_D65(lambda)"], dtype=np.float64)


def window(a, b):
  """Interior slice [a, b] of the 1 nm tables (node-exact resample)."""
  m = (LAM >= a) & (LAM <= b)
  me = (DLAM >= a) & (DLAM <= b)
  return (LAM[m], CX[m], CY[m], CZ[m], DLAM[me], E65[me])


def demand(curve="Rs", angle=0.0, quantity="Lab", reference=(60.0, 10.0, -20.0),
           distance="DeltaE2000", weight=1.0, tables=None, names=False,
           wr=None):
  if names:
    illum, cmf = "D65", "1931_2deg"
  else:
    wl, x, y, z, el, ev = tables
    illum = {"wavelengths": el.tolist(), "values": ev.tolist()}
    cmf = {"wavelengths": wl.tolist(),
           "xyz": np.stack([x, y, z], axis=1).tolist()}
  doc = {"spectral": [], "angular": [],
         "color": [{"curve": curve, "angle": angle,
                    "illuminant": illum, "observer": cmf,
                    "quantity": quantity,
                    "reference": (list(reference) if isinstance(reference, (tuple, list)) else reference),
                    "distance": distance, "weight": weight,
                    "wavelength_range": wr}],
         "cache_size": 128, "tolerance_floor": 1e-12}
  return compile_merit_spec(json.dumps(doc))


def xyz_oracle(row, wl, x, y, z, ev):
  """Pure-Python-loop XYZ (mirrors the Rust summation op-for-op).
  Aligned grids only (no resample) — callers slice interior nodes."""
  n = len(wl)
  dw = [wl[i + 1] - wl[i] for i in range(n - 1)] + [wl[-1] - wl[-2]]
  den = 0.0
  for i in range(n):
    den += ev[i] * y[i] * dw[i]
  k = 1.0 / den
  out = [0.0, 0.0, 0.0]
  for c, cc in enumerate((x, y, z)):
    s = 0.0
    for i in range(n):
      w = row[i] * ev[i] * k * dw[i]
      s += w * cc[i]
    out[c] = s
  return out


def xyy_of(xyz):
  s = xyz[0] + xyz[1] + xyz[2]
  inv = 1.0 / s
  return [xyz[0] * inv, xyz[1] * inv, xyz[1]]


def test_xyy_channels_hex_against_oracle():
  # Strictest test in the batch: fully analytic path, bitwise.
  wl, x, y, z, el, ev = window(500.0, 519.0)
  row = 0.2 + 0.6 * (wl - wl[0]) / (wl[-1] - wl[0])
  ref = (0.35, 0.36, 0.55)
  spec = demand(quantity="XyY", reference=ref, distance="Channels",
                weight=2.0, tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  assert spec.n_residuals() == 1 and len(got) == 1
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  c = xyy_of(xyz)
  f = 2.0 * sum(((a - b) / 1.0) ** 2 for a, b in zip(c, ref))
  assert got[0].hex() == (f ** 0.5).hex()
  assert spec.merit(sim, 1e6).hex() == (got[0] * got[0]).hex()


def test_lab_de00_against_color_pipeline():
  wl, x, y, z, el, ev = window(400.0, 700.0)
  assert len(wl) == 301
  row = 0.5 + 0.3 * np.sin((wl - 400.0) / 300.0 * np.pi)
  ref = (62.0, 18.0, -34.0)
  spec = demand(quantity="Lab", reference=ref, distance="DeltaE2000",
                tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  lab = C.XYZ_to_Lab(np.array([xyz]), illuminant=white)[0]
  d = C.delta_E_CIE2000(np.array([lab]), np.array([ref]), 1.0, 1.0, 1.0, False)[0]
  assert got[0] == pytest.approx(float(d), rel=1e-12)


def test_missing_curve_penalty_parity():
  wl, x, y, z, el, ev = window(500.0, 519.0)
  spec = demand(quantity="XyY", reference=(0.3, 0.3, 0.5),
                distance="Channels",
                tables=(wl, x, y, z, el, ev))
  row = np.full((1, len(wl)), 0.5)
  full = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row})
  assert np.isfinite(spec.merit(full, 1e6))
  empty = sim_curves_from_arrays(np.array([0.0]), wl, {})
  assert spec.merit(empty, 1e6) == 1e6  # one penalty (shared key group)
  with pytest.raises(ValueError, match="Rs"):
    spec.residuals(empty)


def test_n_residuals_counts_color_as_one():
  wl, x, y, z, el, ev = window(500.0, 519.0)
  spec = demand(quantity="XyY", reference=(0.3, 0.3, 0.5),
                distance="Channels", tables=(wl, x, y, z, el, ev))
  assert spec.n_residuals() == 1
  row = np.full((1, len(wl)), 0.5)
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row})
  assert len(spec.residuals(sim)) == 1


def test_named_tables_resolve():
  wl, x, y, z, el, ev = window(500.0, 519.0)
  spec = demand(quantity="Lab", reference=(60.0, 10.0, -20.0),
                distance="DeltaE76", names=True)
  row = np.full((1, len(wl)), 0.5)
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row})
  assert np.isfinite(spec.merit(sim, 1e6))


def _lch_of(lab):
  l, a, b = lab
  import math
  h = math.degrees(math.atan2(b, a))
  return [l, math.hypot(a, b), h if h >= 0 else h + 360.0]


def test_lch_channels_against_color_pipeline():
  wl, x, y, z, el, ev = window(400.0, 700.0)
  row = 0.5 + 0.3 * np.sin((wl - 400.0) / 300.0 * np.pi)
  ref = (62.0, 18.0, 100.0)
  spec = demand(quantity="LCh", reference=ref, distance="DeltaE76",
                tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  lab = C.XYZ_to_Lab(np.array([xyz]), illuminant=white)[0]
  lch = C.Lab_to_LCHab(np.array([lab]))[0]
  # DeltaE76 path converts the LCh ref to Lab; Channels compares in LCh.
  lab_ref = C.LCHab_to_Lab(np.array([ref]))[0]
  d = C.delta_E_CIE1976(np.array([lab]), np.array([lab_ref]), )[0]
  assert got[0] == pytest.approx(float(d), rel=1e-12)


def test_lch_channels_hue_cut_twin():
  # Ref hue 181 deg past the op hue: wrapped gap is 179, raw would be 181.
  wl, x, y, z, el, ev = window(400.0, 700.0)
  row = 0.5 + 0.3 * np.sin((wl - 400.0) / 300.0 * np.pi)
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  lab = C.XYZ_to_Lab(np.array([xyz]), illuminant=white)[0]
  h0 = _lch_of(lab.tolist())[2]
  ref = (60.0, 20.0, h0 + 181.0)
  spec = demand(quantity="LCh", reference=ref, distance="Channels",
                tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)[0]
  lch = C.Lab_to_LCHab(np.array([lab]))[0]
  wrap = lambda d: d - 360.0 * round(d / 360.0)
  expect = sum(((a - b) if i < 2 else wrap(a - b)) ** 2
               for i, (a, b) in enumerate(zip(lch, ref))) ** 0.5
  raw = sum((a - b) ** 2 for a, b in zip(lch, ref)) ** 0.5
  assert got == pytest.approx(expect, rel=1e-12)
  assert abs(got - raw) > 1.0  # wrapped, provably not raw


def test_oklab_channels_hex_against_pipeline():
  wl, x, y, z, el, ev = window(400.0, 700.0)
  row = 0.5 + 0.3 * np.sin((wl - 400.0) / 300.0 * np.pi)
  ref = (0.55, 0.02, -0.03)
  spec = demand(quantity="Oklab", reference=ref, distance="Channels",
                tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  adapted = C.chromatic_adaptation_VonKries(
    np.array([xyz]), white, list(C.REF_WHITE_D65), False)[0]
  ok = C.XYZ_to_Oklab(np.array([adapted]))[0]
  f = sum((a - b) ** 2 for a, b in zip(ok, ref))
  assert got[0].hex() == (f ** 0.5).hex()


def test_y_scalar_hex_against_hand_oracle():
  wl, x, y, z, el, ev = window(500.0, 519.0)
  row = 0.2 + 0.6 * (wl - wl[0]) / (wl[-1] - wl[0])
  spec = demand(quantity="Y", reference=0.45, distance="Channels",
                weight=2.0, tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  assert spec.n_residuals() == 1
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  f = 2.0 * ((xyz[1] - 0.45) / 1.0) ** 2
  assert got[0].hex() == (f ** 0.5).hex()


def _domwl_oracle(xyz, white, locus):
  """Independent ray-vs-locus oracle (same math, numpy scalar ops)."""
  s = sum(xyz)
  ws = sum(white)
  w = [white[0] / ws, white[1] / ws]
  p = [xyz[0] / s, xyz[1] / s]
  d = [p[0] - w[0], p[1] - w[1]]
  if d[0] * d[0] + d[1] * d[1] < 1e-18:
    return (0.0, 0.0)
  fwd, bwd = None, None
  for (ax, ay, la), (bx, by, lb) in zip(locus[:-1], locus[1:]):
    ex, ey = bx - ax, by - ay
    den = d[0] * ey - d[1] * ex
    if abs(den) < 1e-300:
      continue
    awx, awy = ax - w[0], ay - w[1]
    t = (awx * ey - awy * ex) / den
    sg = (awx * d[1] - awy * d[0]) / den
    if not (0.0 <= sg <= 1.0):
      continue
    lam = la + sg * (lb - la)
    if t > 1.0 - 1e-9 and (fwd is None or t < fwd[0]):
      fwd = (t, lam)
    u = -t
    if u > 1e-9 and (bwd is None or u < bwd[0]):
      bwd = (u, lam)
  if fwd is not None:
    return (fwd[1], 1.0 / fwd[0])
  if bwd is not None:
    return (bwd[1], -1.0 / bwd[0])
  raise AssertionError("oracle missed locus")


def test_domwl_against_independent_oracle():
  wl, x, y, z, el, ev = window(400.0, 700.0)
  row = np.exp(-((wl - 550.0) / 15.0) ** 2)
  ref = (550.0, 0.9)
  spec = demand(quantity="DomWl", reference=ref, distance="Channels",
                tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  locus = [(a / (a + b + c), b / (a + b + c), w)
           for a, b, c, w in zip(x.tolist(), y.tolist(), z.tolist(), wl.tolist())
           if a + b + c > 0]
  lam, purity = _domwl_oracle(xyz, white, locus)
  assert lam == pytest.approx(550.0, abs=3.0) and purity > 0.8
  f = ((lam - ref[0]) / 1.0) ** 2 + ((purity - ref[1]) / 1.0) ** 2
  assert got[0] == pytest.approx(f ** 0.5, rel=1e-9)


def test_srgb_channels_against_pipeline():
  # In-gamut gray ramp: kernel == bound pipeline (adapt + unclipped map).
  wl, x, y, z, el, ev = window(400.0, 700.0)
  row = np.linspace(0.05, 0.9, len(wl))
  ref = (0.5, 0.5, 0.5)
  spec = demand(quantity="sRGB", reference=ref, distance="Channels",
                tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  adapted = C.chromatic_adaptation_VonKries(
    np.array([xyz]), white, list(C.REF_WHITE_D65), False)
  rgb = C.XYZ_to_sRGB(adapted, clip=False)[0]
  f = sum((a - b) ** 2 for a, b in zip(rgb, ref))
  assert got[0] == pytest.approx(f ** 0.5, rel=1e-12)


def test_luv_channels_against_pipeline():
  wl, x, y, z, el, ev = window(400.0, 700.0)
  row = 0.5 + 0.3 * np.sin((wl - 400.0) / 300.0 * np.pi)
  ref = (50.0, 10.0, -20.0)
  spec = demand(quantity="Luv", reference=ref, distance="Channels",
                tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  luv = C.XYZ_to_Luv(np.array([xyz]), illuminant=white)[0]
  f = sum((a - b) ** 2 for a, b in zip(luv, ref))
  assert got[0] == pytest.approx(f ** 0.5, rel=1e-12)


def test_xyz_channels_hex_against_hand_oracle():
  wl, x, y, z, el, ev = window(500.0, 519.0)
  row = 0.2 + 0.6 * (wl - wl[0]) / (wl[-1] - wl[0])
  ref = (0.3, 0.4, 0.2)
  spec = demand(quantity="XYZ", reference=ref, distance="Channels",
                weight=2.0, tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  f = 2.0 * sum((a - b) ** 2 for a, b in zip(xyz, ref))
  assert got[0].hex() == (f ** 0.5).hex()


def test_din99_channels_against_numpy():
  # Independent oracle: Lab (via _color) -> DIN99 closed form in numpy.
  import math
  wl, x, y, z, el, ev = window(400.0, 700.0)
  row = 0.5 + 0.3 * np.sin((wl - 400.0) / 300.0 * np.pi)
  ref = (50.0, 5.0, 5.0)
  spec = demand(quantity="Din99", reference=ref, distance="Channels",
                tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  lab = C.XYZ_to_Lab(np.array([xyz]), illuminant=white)[0]
  def din99(lab):
    l, a, b = lab
    l99 = 105.509 * math.log1p(0.0158 * l)
    e = a * math.cos(math.radians(16.0)) + b * math.sin(math.radians(16.0))
    f = 0.7 * (-a * math.sin(math.radians(16.0)) + b * math.cos(math.radians(16.0)))
    g = math.hypot(e, f)
    if g < 1e-12:
      return (l99, 0.0, 0.0)
    c = math.log1p(0.045 * g) / 0.045
    h = math.atan2(f, e)
    return (l99, c * math.cos(h), c * math.sin(h))
  d = din99(lab.tolist())
  f = sum((a - b) ** 2 for a, b in zip(d, ref))
  assert got[0] == pytest.approx(f ** 0.5, rel=1e-12)


def test_white_yellow_against_hand_oracles():
  wl, x, y, z, el, ev = window(500.0, 519.0)
  row = 0.2 + 0.6 * (wl - wl[0]) / (wl[-1] - wl[0])
  xyz = xyz_oracle(row.tolist(), wl.tolist(), x.tolist(), y.tolist(),
                   z.tolist(), ev.tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  # Whiteness pair vs closed form (1e-12; same fractions, both sides).
  wspec = demand(quantity="White", reference=(90.0, 2.0), distance="Channels",
                 tables=(wl, x, y, z, el, ev))
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = wspec.residuals(sim)
  s, ws = sum(xyz), sum(white)
  w = 100.0 * xyz[1] + 800.0 * (white[0] / ws - xyz[0] / s) + 1700.0 * (white[1] / ws - xyz[1] / s)
  tw = 1000.0 * (white[0] / ws - xyz[0] / s) - 650.0 * (white[1] / ws - xyz[1] / s)
  assert got[0] == pytest.approx(((w - 90.0) ** 2 + (tw - 2.0) ** 2) ** 0.5, rel=1e-12)
  # Yellowness scalar vs E313 hand formula (HEX: identical op order).
  yspec = demand(quantity="Yellow", reference=5.0, distance="Channels",
                 tables=(wl, x, y, z, el, ev))
  goty = yspec.residuals(sim)
  yi = 100.0 * (1.3013 * xyz[0] - 1.1498 * xyz[2]) / xyz[1]
  # sqrt(r^2) rounds once vs abs(): 1e-12, not hex.
  assert goty[0] == pytest.approx(abs(yi - 5.0), rel=1e-12)


def test_wavelength_range_windows_the_sample_only():
  # Full tables + wr=[500, 519] vs oracle integrating ONLY that slice.
  # White stays full-range: oracle Lab uses the full-table white.
  wl, x, y, z, el, ev = window(400.0, 700.0)
  row = 0.5 + 0.3 * np.sin((wl - 400.0) / 300.0 * np.pi)
  ref = (60.0, 10.0, -20.0)
  spec = demand(quantity="Lab", reference=ref, distance="DeltaE2000",
                tables=(wl, x, y, z, el, ev), wr=[500.0, 519.0])
  sim = sim_curves_from_arrays(np.array([0.0]), wl, {"Rs": row.reshape(1, -1)})
  got = spec.residuals(sim)
  assert spec.n_residuals() == 1 and len(got) == 1
  m = (wl >= 500.0) & (wl <= 519.0)
  xyz = xyz_oracle(row[m].tolist(), wl[m].tolist(), x[m].tolist(),
                   y[m].tolist(), z[m].tolist(), ev[m].tolist())
  white = xyz_oracle([1.0] * len(wl), wl.tolist(), x.tolist(), y.tolist(),
                     z.tolist(), ev.tolist())
  lab = C.XYZ_to_Lab(np.array([xyz]), illuminant=white)[0]
  tlab = np.asarray(ref)
  d = C.delta_E_CIE2000(np.array([lab]), np.array([ref]), 1.0, 1.0, 1.0, False)[0]
  assert got[0] == pytest.approx(float(d), rel=1e-12)
  # Reversed range refuses LOUD at compile (named, not silent).
  with pytest.raises(Exception, match="lo < hi"):
    demand(quantity="Lab", reference=ref, distance="DeltaE2000",
           tables=(wl, x, y, z, el, ev), wr=[519.0, 500.0])
