"""Deterministic figure rendering — the ONE definition of the settings.

Every figure a run writes must be byte-reproducible, so two runs at the same seeds
can be compared file by file and the figures can take part in the re-run
determinism assertion instead of being excluded from it.

Two matplotlib rcParams carry that guarantee:

``svg.hashsalt``
    The SVG backend derives element ids from this salt. When it is None —
    matplotlib's default — the ids come from ``uuid.uuid4()`` instead, so every
    render of the same figure produces different bytes. PNG output has no ids and
    is unaffected, which is why the symptom is always "the PNG matches and the SVG
    does not".

``path.simplify``
    Vertex simplification is a rendering-time transformation of the plotted data.
    Off, so a vector figure carries the values it was given.

**Why this module exists rather than each figure module setting its own.**
rcParams are process-global mutable state, so these are not four independent
settings — they are one setting that any importer can silently undo. That is not
hypothetical: ``data/figures.py`` calls ``plt.rcdefaults()`` to get a clean
baseline, which reset ``svg.hashsalt`` to None for the whole process. Stage A
draws before Stage B, C and D in every real run, so every vector figure written
after it was non-reproducible — and only Stage B asserted on it, so Stages C and D
carried the same defect unnoticed. The full suite caught it and targeted runs did
not, purely because targeted runs happened to execute Stage B first.

So the values live here once, in the Lead-owned package, and a module that resets
rcParams re-applies them immediately afterwards. Adding a figure module means
calling :func:`apply_deterministic_rcparams` at import; a module that forgets is
caught by ``tests/integration/test_figure_determinism.py``.

This is a RENDERING guarantee only. No scientific value, threshold or computation
is affected by either parameter.
"""
from __future__ import annotations

#: The frozen rendering settings. Not scientific constants — they govern how a
#: figure is serialized, never what it plots.
DETERMINISTIC_RCPARAMS: dict[str, object] = {
    "svg.hashsalt": "hotspot3d",
    "path.simplify": False,
}

#: SVG metadata that must be suppressed: matplotlib writes a ``<dc:date>`` element
#: with the wall-clock time unless the key is explicitly set to None, and a
#: timestamp makes two otherwise identical renders differ.
SVG_METADATA: dict[str, None] = {"Date": None}


def apply_deterministic_rcparams() -> None:
    """Apply the settings to the current matplotlib process state.

    Call at import in every module that writes figures, and again immediately
    after any ``rcdefaults()``/``rc_context`` reset, which clears them.

    Imports matplotlib lazily so importing :mod:`hotspot3d.utils` never requires
    it; callers that treat matplotlib as optional already guard the import and
    should call this only once it has succeeded.
    """
    import matplotlib

    matplotlib.rcParams.update(DETERMINISTIC_RCPARAMS)
