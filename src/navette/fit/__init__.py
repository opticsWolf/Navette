# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""navette.fit: fit any scalar of a model against measured data.

A fit is a :class:`FitProblem`: a program's model (its materials, groups,
structures and architect), the samples it was measured on
(:class:`Environment`), the data (:class:`Dataset`, :class:`NkDataset`),
and the scalars that may move (:class:`Parameter`). :func:`fit` runs it,
in one stage or by a :class:`Strategy`, and returns a :class:`FitReport`::

    from navette.fit import Dataset, Environment, FitProblem, Parameter, fit

    problem = FitProblem(
        model=sections,                     # a program's sections mapping
        environments=[Environment("sample", ambient="air", substrate="glass",
                                  stack=[{"structure": "ar"}])],
        data=[Dataset("R", "R_avg", wl, measured, sigma=0.002, angle_deg=8)],
        parameters=[Parameter(select="kind=thickness")],
    )
    report = fit(problem, progress=print)
    print(report.summary())

The problem is a program document with a ``fit`` section
(:meth:`FitProblem.to_document`), so it saves and loads as a program file
(:meth:`FitProblem.save`, :meth:`FitProblem.load`). The engine is Rust
(``navette::fit``); this module builds the request and reads the reply.

Strategies, generated strategies (:func:`orderings`, :func:`groupings`,
:func:`cumulative`), the :func:`tournament` and the :func:`benchmark`
compare ways of fitting one problem. :func:`to_solver` reads instrument
data in a named convention (:func:`conventions`). :func:`fit_dispersion`
and :func:`table_to_model` fit a material model to optical constants.

Addresses name the scalars: ``layer[structure:row].thickness`` (rows count
films from 0), ``group[G].n_factor``, ``material[TiO2].A``. Selectors pick
them: ``*``, ``kind=thickness``, ``tag=H``, ``set:name``, or a glob.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from navette.config.program import PROGRAM_SCHEMA_VERSION

__all__ = [
    "Parameter",
    "Tie",
    "Environment",
    "Dataset",
    "NkDataset",
    "Stage",
    "Strategy",
    "FitProblem",
    "FittedParameter",
    "Curve",
    "FitReport",
    "TournamentEntry",
    "BenchmarkEntry",
    "BenchmarkReport",
    "Truth",
    "fit",
    "evaluate",
    "orderings",
    "groupings",
    "cumulative",
    "generate",
    "tournament",
    "benchmark",
    "draw_truth",
    "conventions",
    "to_solver",
    "from_solver",
    "fit_dispersion",
    "table_to_model",
    "plot_fit",
]

ArrayLike = Union[Sequence[float], np.ndarray]


def _native():
    from navette import _fit
    return _fit


def _floats(x: ArrayLike) -> List[float]:
    return [float(v) for v in np.asarray(x, dtype=np.float64).ravel()]


def _drop_none(d: Dict[str, Any]) -> Dict[str, Any]:
    return {k: v for k, v in d.items() if v is not None}


# ---------------------------------------------------------------------------
# The request
# ---------------------------------------------------------------------------


@dataclass
class Parameter:
    """One scalar (``addr``) or every model scalar a selector picks (``select``).

    Unset bounds and scale take the defaults of the scalar's kind; an unset
    ``value`` reads the model. ``log`` fits ``ln x``. ``prior`` is
    ``(value, sigma)``: a residual that holds the parameter near ``value``.
    A ``select`` expands over film thicknesses, group factors and summands,
    and material numbers, skipping addresses declared elsewhere.
    """

    addr: Optional[str] = None
    select: Optional[str] = None
    value: Optional[float] = None
    lb: Optional[float] = None
    ub: Optional[float] = None
    scale: Optional[float] = None
    log: bool = False
    fixed: bool = False
    tags: Sequence[str] = ()
    prior: Optional[Tuple[float, float]] = None

    def to_json(self) -> Dict[str, Any]:
        if (self.addr is None) == (self.select is None):
            raise ValueError("Parameter: give exactly one of addr and select")
        d = _drop_none({
            "addr": self.addr, "select": self.select, "value": self.value,
            "lb": self.lb, "ub": self.ub, "scale": self.scale,
            "prior": list(self.prior) if self.prior is not None else None,
        })
        if self.log:
            d["log"] = True
        if self.fixed:
            d["fixed"] = True
        if self.tags:
            d["tags"] = list(self.tags)
        return d

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> "Parameter":
        d = dict(d)
        if d.get("prior") is not None:
            d["prior"] = tuple(d["prior"])
        d["tags"] = tuple(d.get("tags", ()))
        return cls(**d)


@dataclass
class Tie:
    """``target = factor * driver + offset`` for each ``(target, factor, offset)``."""

    driver: str
    targets: Sequence[Union[str, Tuple[str, float, float]]]

    def to_json(self) -> Dict[str, Any]:
        out = []
        for t in self.targets:
            if isinstance(t, str):
                out.append({"target": t})
            else:
                target, factor, offset = t
                out.append({"target": target, "factor": factor, "offset": offset})
        return {"driver": self.driver, "targets": out}

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> "Tie":
        return cls(d["driver"], [
            (t["target"], t.get("factor", 1.0), t.get("offset", 0.0)) for t in d["targets"]
        ])


@dataclass
class Environment:
    """A sample: a chain of the model's structures and held films.

    ``stack`` items are ``{"structure": name}`` (with optional
    ``"inverted"`` and ``"repeat"``) or ``{"layer": material,
    "thickness": nm}``, between an optional ambient and substrate.
    """

    name: str
    ambient: Optional[str] = None
    substrate: Optional[str] = None
    stack: Sequence[Mapping[str, Any]] = ()

    def to_json(self) -> Dict[str, Any]:
        return _drop_none({"name": self.name, "ambient": self.ambient,
                           "substrate": self.substrate,
                           "stack": [dict(s) for s in self.stack]})

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> "Environment":
        return cls(d["name"], d.get("ambient"), d.get("substrate"), list(d.get("stack", [])))


@dataclass
class Dataset:
    """One measured spectrum in solver conventions (fractions, radians for
    Psi/Delta, the solver's Delta sign).

    ``channel`` is a solver curve (``Rs``, ``Tp``, ``R_avg``, ``Psi_R``,
    ``Delta_R``, ...); ``sigma`` one number or one per point. ``env``
    defaults to the problem's first environment. Instrument data goes
    through :func:`to_solver` first.
    """

    id: str
    channel: str
    wavelengths: ArrayLike
    values: ArrayLike
    sigma: Union[float, ArrayLike] = 1.0
    angle_deg: float = 0.0
    env: Optional[str] = None
    weight: float = 1.0
    window: Optional[Tuple[float, float]] = None
    source: Optional[str] = None

    def to_json(self) -> Dict[str, Any]:
        sigma = float(self.sigma) if np.ndim(self.sigma) == 0 else _floats(self.sigma)
        d = _drop_none({
            "id": self.id, "env": self.env, "channel": self.channel,
            "angle_deg": float(self.angle_deg),
            "wavelengths": _floats(self.wavelengths), "values": _floats(self.values),
            "sigma": sigma,
            "window": list(self.window) if self.window is not None else None,
            "source": self.source,
        })
        if self.weight != 1.0:
            d["weight"] = float(self.weight)
        return d

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> "Dataset":
        if "angle_rad" in d and d["angle_rad"] is not None:
            angle = float(np.degrees(d["angle_rad"]))
        else:
            angle = float(d.get("angle_deg", 0.0))
        sigma = d["sigma"]
        return cls(
            id=d["id"], channel=d["channel"],
            wavelengths=np.asarray(d["wavelengths"], dtype=np.float64),
            values=np.asarray(d["values"], dtype=np.float64),
            sigma=sigma if np.ndim(sigma) == 0 else np.asarray(sigma, dtype=np.float64),
            angle_deg=angle, env=d.get("env"), weight=d.get("weight", 1.0),
            window=tuple(d["window"]) if d.get("window") else None,
            source=d.get("source"),
        )


@dataclass
class NkDataset:
    """Measured optical constants of one material of the model.

    ``columns`` maps ``"n"``, ``"k"``, ``"eps1"`` or ``"eps2"`` to
    ``(values, sigma)``.
    """

    id: str
    material: str
    wavelengths: ArrayLike
    columns: Mapping[str, Tuple[ArrayLike, Union[float, ArrayLike]]]
    weight: float = 1.0
    window: Optional[Tuple[float, float]] = None

    def to_json(self) -> Dict[str, Any]:
        cols = []
        for quantity, (values, sigma) in self.columns.items():
            s = float(sigma) if np.ndim(sigma) == 0 else _floats(sigma)
            cols.append({"quantity": quantity, "values": _floats(values), "sigma": s})
        d = _drop_none({
            "id": self.id, "material": self.material,
            "wavelengths": _floats(self.wavelengths), "columns": cols,
            "window": list(self.window) if self.window is not None else None,
        })
        if self.weight != 1.0:
            d["weight"] = float(self.weight)
        return d

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> "NkDataset":
        return cls(
            d["id"], d["material"], np.asarray(d["wavelengths"], dtype=np.float64),
            {c["quantity"]: (np.asarray(c["values"], dtype=np.float64), c["sigma"])
             for c in d["columns"]},
            d.get("weight", 1.0), tuple(d["window"]) if d.get("window") else None,
        )


@dataclass
class Stage:
    """One stage: the selectors it frees, what it scores, how it solves.

    ``kind`` is ``None`` (least squares), ``{"scan": {"param": sel,
    "points": n}}``, ``{"multi_start": {"starts": n, "sampler":
    "lhs", "seed": 0}}`` or ``{"global": {"method": "de", "budget": n,
    "seed": 0}}`` (``"de"``, ``"de_inject"`` or ``"cma_es"``; needs the
    ``opt-basin`` build). ``repeat`` is ``{"until_no_gain": {"rel": 1e-3,
    "max": 5}}``; ``accept`` ``"if_better"`` rolls a stage back that did
    not lower its chi2. ``lsq`` overrides solver settings (``backend``,
    ``max_iterations``, ``max_evals``, ``ftol``, ``xtol``, ``gtol``,
    ``jacobian``: ``"analytic"`` or ``"fd"``).
    """

    name: str
    free: Sequence[str]
    datasets: Optional[Sequence[str]] = None
    weights: Optional[Mapping[str, float]] = None
    window: Optional[Tuple[float, float]] = None
    kind: Optional[Mapping[str, Any]] = None
    lsq: Optional[Mapping[str, Any]] = None
    repeat: Optional[Mapping[str, Any]] = None
    accept: Optional[str] = None

    def to_json(self) -> Dict[str, Any]:
        return _drop_none({
            "name": self.name, "free": list(self.free),
            "datasets": list(self.datasets) if self.datasets is not None else None,
            "weights": [[k, float(v)] for k, v in self.weights.items()] if self.weights else None,
            "window": list(self.window) if self.window is not None else None,
            "kind": dict(self.kind) if self.kind is not None else None,
            "lsq": dict(self.lsq) if self.lsq is not None else None,
            "repeat": dict(self.repeat) if self.repeat is not None else None,
            "accept": self.accept,
        })

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> "Stage":
        kind = d.get("kind")
        repeat = d.get("repeat")
        return cls(
            d["name"], list(d["free"]), d.get("datasets"),
            {k: v for k, v in d["weights"]} if d.get("weights") else None,
            tuple(d["window"]) if d.get("window") else None,
            None if kind in (None, "lsq") else kind,
            d.get("lsq") or None,
            None if repeat in (None, "once") else repeat,
            None if d.get("accept") in (None, "always") else d.get("accept"),
        )


@dataclass
class Strategy:
    """An ordered list of stages."""

    name: str
    stages: Sequence[Stage]

    def to_json(self) -> Dict[str, Any]:
        return {"name": self.name, "stages": [s.to_json() for s in self.stages]}

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> "Strategy":
        return cls(d["name"], [Stage.from_json(s) for s in d["stages"]])


@dataclass
class FitProblem:
    """What to fit: a model, its samples, the data, the scalars that move.

    ``model`` is a program's ``sections`` mapping (``materials``,
    ``groups``, ``structures``, ``architect``) without a ``fit`` entry.
    ``sets`` names unions of selectors for stages (``{"H": ["tag=H"]}``).
    ``lsq`` holds the solver settings every stage starts from.
    """

    model: Mapping[str, Any]
    data: Sequence[Dataset] = ()
    parameters: Sequence[Parameter] = ()
    environments: Sequence[Environment] = ()
    nk_data: Sequence[NkDataset] = ()
    ties: Sequence[Tie] = ()
    sets: Mapping[str, Sequence[str]] = field(default_factory=dict)
    strategy: Optional[Strategy] = None
    coherence_mode: int = 0
    lsq: Mapping[str, Any] = field(default_factory=dict)
    name: Optional[str] = None

    def fit_section(self, strategy: Optional[Strategy] = None) -> Dict[str, Any]:
        """The program's ``fit`` section; ``strategy`` replaces the problem's."""
        strategy = strategy if strategy is not None else self.strategy
        section: Dict[str, Any] = {
            "parameters": [p.to_json() for p in self.parameters],
            "coherence_mode": int(self.coherence_mode),
            "lsq": dict(self.lsq),
        }
        if self.environments:
            section["environments"] = [e.to_json() for e in self.environments]
        if self.data:
            section["data"] = [d.to_json() for d in self.data]
        if self.nk_data:
            section["nk_data"] = [d.to_json() for d in self.nk_data]
        if self.ties:
            section["ties"] = [t.to_json() for t in self.ties]
        if self.sets:
            section["sets"] = [{"name": k, "members": list(v)} for k, v in self.sets.items()]
        if strategy is not None:
            section["strategy"] = strategy.to_json()
        return section

    def to_document(self, strategy: Optional[Strategy] = None) -> Dict[str, Any]:
        """The program document (schema v4) this problem is."""
        if "fit" in self.model:
            raise ValueError("FitProblem.model holds a 'fit' section; pass the model only")
        doc: Dict[str, Any] = {
            "schema_version": PROGRAM_SCHEMA_VERSION,
            "kind": "program",
            "sections": {**dict(self.model), "fit": self.fit_section(strategy)},
        }
        if self.name is not None:
            doc["name"] = self.name
        return doc

    def to_json(self, strategy: Optional[Strategy] = None) -> str:
        return json.dumps(self.to_document(strategy))

    def save(self, path: Union[str, Path]) -> None:
        """Write the problem as a JSON program file."""
        Path(path).write_text(json.dumps(self.to_document(), indent=1), encoding="utf-8")

    @classmethod
    def from_document(cls, doc: Mapping[str, Any]) -> "FitProblem":
        """The problem a program document's ``fit`` section describes."""
        sections = dict(doc["sections"])
        section = sections.pop("fit", None)
        if section is None:
            raise ValueError("program document has no 'fit' section")
        strategy = section.get("strategy")
        return cls(
            model=sections,
            data=[Dataset.from_json(d) for d in section.get("data", [])],
            parameters=[Parameter.from_json(p) for p in section.get("parameters", [])],
            environments=[Environment.from_json(e) for e in section.get("environments", [])],
            nk_data=[NkDataset.from_json(d) for d in section.get("nk_data", [])],
            ties=[Tie.from_json(t) for t in section.get("ties", [])],
            sets={s["name"]: list(s["members"]) for s in section.get("sets", [])},
            strategy=Strategy.from_json(strategy) if strategy else None,
            coherence_mode=section.get("coherence_mode", 0),
            lsq=dict(section.get("lsq") or {}),
            name=doc.get("name"),
        )

    @classmethod
    def load(cls, path: Union[str, Path]) -> "FitProblem":
        """Read a program file (JSON or YAML) with a ``fit`` section."""
        from navette.config.program import load_document
        kind, name, payload = load_document(path)
        if kind != "program":
            raise ValueError(f"{path}: a fit needs a program document, got kind {kind!r}")
        return cls.from_document({"name": name, "sections": payload})


# ---------------------------------------------------------------------------
# The reply
# ---------------------------------------------------------------------------


@dataclass
class FittedParameter:
    addr: str
    start: Optional[float]
    value: Optional[float]
    std_error: Optional[float]
    lb: Optional[float]
    ub: Optional[float]
    fixed: bool
    at_bound: bool


@dataclass
class Curve:
    """A dataset at the solution: what a plot of the fit draws.

    ``label`` is the channel (spectra) or the quantity (nk data).
    """

    id: str
    kind: str
    label: str
    wavelengths: np.ndarray
    measured: np.ndarray
    simulated: np.ndarray
    sigma: np.ndarray
    angle_deg: Optional[float] = None
    env: Optional[str] = None
    material: Optional[str] = None

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> "Curve":
        arr = lambda k: np.asarray(d[k], dtype=np.float64)  # noqa: E731
        return cls(
            d["id"], d["kind"], d.get("channel") or d.get("quantity"),
            arr("wavelengths"), arr("measured"), arr("simulated"), arr("sigma"),
            d.get("angle_deg"), d.get("env"), d.get("material"),
        )


def _nan(x: Optional[float]) -> float:
    return float("nan") if x is None else float(x)


def _square(v: Optional[Sequence[Optional[float]]], n: int) -> Optional[np.ndarray]:
    if v is None:
        return None
    return np.asarray([_nan(x) for x in v], dtype=np.float64).reshape(n, n)


@dataclass
class FitReport:
    """A fit's result.

    ``chi2`` sums ``(weight * (sim - value) / sigma)^2``; ``covariance`` and
    ``correlation`` are over the free parameters in problem order, in their
    own units. ``stages`` lists what each stage did when a strategy ran.
    ``curves`` hold every dataset at the solution.
    """

    parameters: List[FittedParameter]
    datasets: List[Dict[str, Any]]
    chi2: float
    chi2_reduced: float
    aic: float
    bic: float
    n_residuals: int
    n_free: int
    covariance: Optional[np.ndarray]
    correlation: Optional[np.ndarray]
    covariance_note: Optional[str]
    termination: str
    iterations: int
    evals: int
    backend: str
    curves: List[Curve] = field(default_factory=list)
    strategy: Optional[str] = None
    stages: List[Dict[str, Any]] = field(default_factory=list)

    @classmethod
    def _from_report(cls, r: Mapping[str, Any]) -> "FitReport":
        n = r["n_free"]
        return cls(
            parameters=[FittedParameter(**p) for p in r["parameters"]],
            datasets=list(r["datasets"]),
            chi2=_nan(r["chi2"]), chi2_reduced=_nan(r["chi2_reduced"]),
            aic=_nan(r["aic"]), bic=_nan(r["bic"]),
            n_residuals=r["n_residuals"], n_free=n,
            covariance=_square(r["covariance"], n), correlation=_square(r["correlation"], n),
            covariance_note=r["covariance_note"], termination=r["termination"],
            iterations=r["iterations"], evals=r["evals"], backend=r["backend"],
        )

    @classmethod
    def from_json(cls, d: Mapping[str, Any]) -> "FitReport":
        rep = cls._from_report(d["report"])
        rep.curves = [Curve.from_json(c) for c in d.get("curves", [])]
        s = d.get("strategy")
        if s is not None:
            rep._adopt_strategy(s)
        return rep

    def _adopt_strategy(self, s: Mapping[str, Any]) -> None:
        # The end report is an evaluation over everything any stage freed;
        # how the solver stopped is the last stage's.
        self.strategy = s["name"]
        self.stages = list(s["stages"])
        self.evals = s["evals"]
        if self.stages:
            self.termination = self.stages[-1]["termination"]
            self.iterations = sum(st["iterations"] for st in self.stages)

    def parameter(self, addr: str) -> FittedParameter:
        for p in self.parameters:
            if p.addr == addr:
                return p
        raise KeyError(addr)

    @property
    def values(self) -> Dict[str, float]:
        """``{addr: value}`` for every parameter."""
        return {p.addr: _nan(p.value) for p in self.parameters}

    def summary(self) -> str:
        """A plain-text table of the parameters and the fit statistics."""
        lines = [f"chi2 {self.chi2:.6g}  chi2_red {self.chi2_reduced:.4g}  "
                 f"n_free {self.n_free}  m {self.n_residuals}  evals {self.evals}  "
                 f"({self.termination})"]
        if self.strategy is not None:
            for s in self.stages:
                mark = "" if s["accepted"] else "  (rolled back)"
                lines.append(f"  stage {s['name']!r}: {len(s['free'])} free, chi2 "
                             f"{s['chi2_start']:.4g} -> {s['chi2']:.4g}, "
                             f"{s['evals']} evals{mark}")
        width = max((len(p.addr) for p in self.parameters), default=10)
        for p in self.parameters:
            err = "" if p.std_error is None else f" +- {p.std_error:.3g}"
            flag = "  fixed" if p.fixed else ("  AT BOUND" if p.at_bound else "")
            lines.append(f"  {p.addr:<{width}}  {_nan(p.value):12.6g}{err}{flag}")
        if self.covariance_note:
            lines.append(f"  no covariance: {self.covariance_note}")
        return "\n".join(lines)


ProgressFn = Callable[[Dict[str, Any]], Any]


def fit(problem: FitProblem, strategy: Optional[Strategy] = None,
        progress: Optional[ProgressFn] = None) -> FitReport:
    """Fit ``problem``, by ``strategy`` (default: the problem's own; with
    none, every free parameter in one stage named ``fit``).

    ``progress(event)`` gets a dict per stage start (``event ==
    "stage_start"``, ``stage``, ``name``, ``free``), per lower chi2 in a
    stage (``"improved"``, ``stage``, ``evals``, ``chi2``) and per stage end
    (``"stage_end"``, ``chi2``, ``accepted``). It may be called from
    several threads during a multi-start stage. Raising in it stops the
    fit; the error comes back as a ``ValueError`` naming it.
    """
    out = _native().fit_program_json(problem.to_json(strategy), progress)
    return FitReport.from_json(json.loads(out))


def evaluate(problem: FitProblem) -> FitReport:
    """The report and curves at the start values, nothing moved."""
    return FitReport.from_json(json.loads(_native().evaluate_program_json(problem.to_json())))


# ---------------------------------------------------------------------------
# Strategies compared
# ---------------------------------------------------------------------------

StrategyLike = Union[Strategy, Mapping[str, Any]]


def orderings(sets: Sequence[str], sample: Optional[Tuple[int, int]] = None,
              polish: bool = True) -> Dict[str, Any]:
    """Every order of ``sets``, one stage per set, then all of them free
    (``polish``); ``sample=(count, seed)`` takes a seeded sample of orders."""
    return {"orderings": _drop_none({"sets": list(sets), "polish": polish,
                                     "sample": list(sample) if sample else None})}


def groupings(variants: Sequence[Sequence[Sequence[str]]], polish: bool = True) -> Dict[str, Any]:
    """One strategy per variant; a variant is a list of blocks of sets freed together."""
    return {"groupings": {"variants": [[list(b) for b in v] for v in variants],
                          "polish": polish}}


def cumulative(sets: Sequence[str], sample: Optional[Tuple[int, int]] = None) -> Dict[str, Any]:
    """In each order of ``sets``: free the first, then the first two, and so on."""
    return {"cumulative": _drop_none({"sets": list(sets),
                                      "sample": list(sample) if sample else None})}


def _strategies_json(strategies: Sequence[StrategyLike]) -> str:
    return json.dumps([s.to_json() if isinstance(s, Strategy) else dict(s) for s in strategies])


def generate(strategies: Sequence[StrategyLike], cap: Optional[int] = None
             ) -> Tuple[List[Strategy], List[str]]:
    """Expand strategies and generators into strategies, with notes of what
    was sampled (orders above ``cap``, default 24, are sampled)."""
    n = _native()
    cap = n.DEFAULT_ORDERING_CAP if cap is None else cap
    out = json.loads(n.generate_json_text(_strategies_json(strategies), cap))
    return [Strategy.from_json(s) for s in out["strategies"]], list(out["notes"])


@dataclass
class TournamentEntry:
    name: str
    rank: int
    chi2: float
    chi2_reduced: float
    bic: float
    rel_uncertainty: float
    evals: int
    time_s: float
    error: Optional[str]
    stages: List[Dict[str, Any]]
    report: Optional[FitReport]


def tournament(problem: FitProblem, strategies: Sequence[StrategyLike],
               cap: Optional[int] = None) -> Tuple[List[TournamentEntry], List[str]]:
    """Run every strategy on ``problem`` (in parallel) and rank by final
    chi2, then evaluations, then name; failed strategies last."""
    n = _native()
    cap = n.DEFAULT_ORDERING_CAP if cap is None else cap
    out = json.loads(n.tournament_json_text(problem.to_json(), _strategies_json(strategies), cap))
    entries = []
    for e in out["entries"]:
        r = e["report"]
        rep = None
        if r is not None:
            rep = FitReport._from_report(r["report"])
            rep._adopt_strategy(r)
        entries.append(TournamentEntry(
            e["name"], e["rank"], _nan(e["chi2"]), _nan(e["chi2_reduced"]), _nan(e["bic"]),
            _nan(e["rel_uncertainty"]), e["evals"], e["time_s"], e["error"],
            list(r["stages"]) if r else [], rep))
    return entries, list(out["notes"])


@dataclass
class BenchmarkEntry:
    name: str
    rank: int
    mean_recovery: float
    worst_recovery: float
    mean_chi2_reduced: float
    mean_evals: float
    successes: int
    failures: int
    runs: List[Dict[str, Any]]


@dataclass
class BenchmarkReport:
    seeds: List[int]
    n_free: int
    entries: List[BenchmarkEntry]
    notes: List[str]

    def summary(self) -> str:
        lines = [f"{'rank':>4}  {'strategy':<30} {'recovery':>10} {'worst':>10} "
                 f"{'chi2_red':>9} {'in 2sd':>6} {'failed':>6} {'evals':>8}"]
        for e in self.entries:
            lines.append(f"{e.rank:>4}  {e.name:<30} {e.mean_recovery:>10.3e} "
                         f"{e.worst_recovery:>10.3e} {e.mean_chi2_reduced:>9.3f} "
                         f"{e.successes:>3}/{len(self.seeds):<2} {e.failures:>6} "
                         f"{e.mean_evals:>8.0f}")
        return "\n".join(lines)


def benchmark(problem: FitProblem, strategies: Sequence[StrategyLike], seeds: Sequence[int],
              cap: Optional[int] = None) -> BenchmarkReport:
    """Every strategy on a truth drawn per seed from the groups' error
    laws (:func:`draw_truth`), from the problem's starts; ranked by
    failures, then recovery of the truth, then evaluations. Deterministic
    per seed."""
    n = _native()
    cap = n.DEFAULT_ORDERING_CAP if cap is None else cap
    out = json.loads(n.benchmark_json_text(problem.to_json(), _strategies_json(strategies),
                                           [int(s) for s in seeds], cap))
    entries = [BenchmarkEntry(
        e["name"], e["rank"], _nan(e["mean_recovery"]), _nan(e["worst_recovery"]),
        _nan(e["mean_chi2_reduced"]), _nan(e["mean_evals"]), e["successes"],
        e["failures"], list(e["runs"])) for e in out["entries"]]
    return BenchmarkReport(list(out["seeds"]), out["n_free"], entries, list(out["notes"]))


@dataclass
class Truth:
    """A truth drawn from the groups' error laws and its simulated data."""

    seed: int
    values: Dict[str, float]
    data: List[Dataset]


def draw_truth(problem: FitProblem, seed: int) -> Truth:
    """Draw every free parameter from its group's error law, then simulate
    every dataset on that truth with noise at its sigma. Deterministic per
    seed; ``data`` replaces the problem's to fit it."""
    out = json.loads(_native().truth_json_text(problem.to_json(), int(seed)))
    return Truth(out["seed"], {a: v for a, v in out["values"]},
                 [Dataset.from_json(d) for d in out["data"]])


# ---------------------------------------------------------------------------
# Instrument data and optical constants
# ---------------------------------------------------------------------------


def conventions() -> Dict[str, Dict[str, Any]]:
    """The preset instrument conventions, by name."""
    return json.loads(_native().conventions_json_text())


def _convention(convention: Union[str, Mapping[str, Any]]) -> str:
    return convention if isinstance(convention, str) else json.dumps(dict(convention))


def to_solver(raw: Mapping[str, Any], convention: Union[str, Mapping[str, Any]]) -> List[Dataset]:
    """Read one instrument record in ``convention`` (a preset name or a
    convention mapping) into solver-convention datasets: one for a
    spectrum, ``<id>.psi`` and ``<id>.delta`` for an ellipsometric record.

    ``raw`` holds ``id``, optionally ``env``, ``quantity`` (``{"spectrum": "R_avg"}``
    or ``{"ellipsometry": "R"}``), ``angle`` and ``abscissa`` in the
    convention's units, ``columns`` and, for a sigma column, ``sigma``.
    """
    raw = dict(raw)
    raw["abscissa"] = _floats(raw["abscissa"])
    raw["columns"] = [_floats(c) for c in raw["columns"]]
    raw["sigma"] = [_floats(c) for c in raw.get("sigma", [])]
    out = _native().to_solver_json_text(json.dumps(raw), _convention(convention))
    return [Dataset.from_json(d) for d in json.loads(out)]


def from_solver(id: str, data: Sequence[Dataset],
                convention: Union[str, Mapping[str, Any]]) -> Dict[str, Any]:
    """Write solver-convention data (one spectrum, or a Psi and a Delta) as
    the instrument would: the inverse of :func:`to_solver`, sigma aside."""
    text = json.dumps([d.to_json() for d in data])
    return json.loads(_native().from_solver_json_text(id, text, _convention(convention)))


def fit_dispersion(spec: Mapping[str, Any], nk_data: Sequence[NkDataset],
                   parameters: Optional[Sequence[Parameter]] = None,
                   lsq: Optional[Mapping[str, Any]] = None
                   ) -> Tuple[Dict[str, Any], FitReport]:
    """Fit one material spec (``{"model": ..., "params": ...}``) to measured
    optical constants; every dataset names the same material, which is the
    spec's name in addresses (``material[<name>].A``). ``parameters=None``
    frees every number of the spec. Returns the fitted spec and the report.
    """
    params = None if parameters is None else json.dumps([p.to_json() for p in parameters])
    out = json.loads(_native().fit_dispersion_json_text(
        json.dumps(dict(spec)), json.dumps([d.to_json() for d in nk_data]), params,
        json.dumps(dict(lsq or {}))))
    return out["spec"], FitReport.from_json(out)


def table_to_model(table: Mapping[str, Any], candidates: Sequence[Mapping[str, Any]],
                   lsq: Optional[Mapping[str, Any]] = None) -> List[Dict[str, Any]]:
    """Fit each candidate spec to a ``Table`` spec's n (and k) and rank them
    by BIC, best first; failed candidates last with their ``error``."""
    out = json.loads(_native().table_to_model_json_text(
        json.dumps(dict(table)), json.dumps([dict(c) for c in candidates]),
        json.dumps(dict(lsq or {}))))
    for c in out:
        if c["report"] is not None:
            c["report"] = FitReport._from_report(c["report"])
    return out


def plot_fit(report: FitReport, axes=None):
    """Measured (with sigma) and simulated curves, one panel per dataset.

    Needs matplotlib; returns the figure.
    """
    import matplotlib.pyplot as plt
    n = len(report.curves)
    if axes is None:
        fig, axes = plt.subplots(n, 1, figsize=(7, 2.6 * n), squeeze=False)
        axes = axes[:, 0]
    else:
        fig = axes[0].figure
    for ax, c in zip(axes, report.curves):
        ax.errorbar(c.wavelengths, c.measured, yerr=c.sigma, fmt=".", ms=3,
                    label="measured", alpha=0.6)
        ax.plot(c.wavelengths, c.simulated, "-", label="fit")
        where = f" @ {c.angle_deg:g} deg" if c.angle_deg is not None else ""
        ax.set_title(f"{c.id}: {c.label}{where}", fontsize=9)
        ax.set_xlabel("wavelength (nm)")
        ax.legend(fontsize=8)
    fig.tight_layout()
    return fig
