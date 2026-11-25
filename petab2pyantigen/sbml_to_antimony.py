"""SBML -> Antimony through Tellurium, without the clutter.

Tellurium's Antimony export is faithful but noisy.  Two kinds of noise matter:

* ``X is "display name"`` lines -- one per reaction, plus ``time_unit is "time"`` and
  ``<model> is "<model>"`` -- and the ``model isDescribedBy / created / creatorN ...``
  annotation block.  They name things; they define nothing.
* ``X is Y`` between identifiers, which is what SBML ``comp`` modules turn into: every
  port of every submodel is tied to a parameter or species of its parent.  A PyAntiGen
  module writes each reaction with exactly the names it should have, so none of this is
  needed.  The clean fix is at the source: flatten the ``comp`` hierarchy with libsbml
  before exporting, so no module and no ``is`` is ever created.  If an identifier ``is``
  line survives anyway it is reported, never silently rewritten.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import libsbml

_DISPLAY_NAME = re.compile(r'^\s*[\w.]+\s+is\s+"[^"]*"\s*;?\s*$')
_IDENT_IS = re.compile(r'^\s*([\w.]+)\s+is\s+([\w.]+)\s*;?\s*$')
_ANNOTATION = re.compile(r'^\s*model\s+(isDescribedBy|created|modified|creator\d*\.\w+|notes|is|isVersionOf|hasPart)\b')


@dataclass
class AntimonyResult:
    text: str                                              # cleaned Antimony
    raw: str                                               # what Tellurium returned
    removed: dict = field(default_factory=dict)            # category -> number of lines removed
    leftover_is: list = field(default_factory=list)        # identifier `is` lines that remain


def flatten_comp(sbml_path: str | Path) -> str:
    """Return the SBML text with any ``comp`` submodels flattened into one model."""
    doc = libsbml.readSBMLFromFile(str(sbml_path))
    if doc.isPackageEnabled("comp"):
        props = libsbml.ConversionProperties()
        props.addOption("flatten comp", True)
        if doc.convert(props) != libsbml.LIBSBML_OPERATION_SUCCESS:
            raise RuntimeError("libsbml could not flatten the comp package: "
                               + doc.getErrorLog().toString())
    return libsbml.writeSBMLToString(doc)


def clean_antimony(text: str) -> tuple[str, dict, list]:
    """Drop display names and annotations; return (text, counts removed, identifier-`is` leftovers)."""
    out, removed, leftover = [], {"display names": 0, "annotations": 0, "notes": 0}, []
    in_notes = False
    for line in text.splitlines():
        if in_notes:                                      # model notes ``` ... ```
            removed["notes"] += 1
            if "```" in line:
                in_notes = False
            continue
        if _ANNOTATION.match(line) and "notes" in line:
            removed["notes"] += 1
            in_notes = line.count("```") == 1
            continue
        if _ANNOTATION.match(line):
            removed["annotations"] += 1
            continue
        if _DISPLAY_NAME.match(line):
            removed["display names"] += 1
            continue
        if _IDENT_IS.match(line):
            leftover.append(line.strip())
        if line.strip() in ("// Notes:", "// Display Names:"):
            continue
        out.append(line)
    cleaned = re.sub(r"\n{3,}", "\n\n", "\n".join(out)).strip() + "\n"
    return cleaned, {k: v for k, v in removed.items() if v}, leftover


def sbml_to_antimony(sbml_path: str | Path, flatten: bool = True) -> AntimonyResult:
    """Convert an SBML file to clean Antimony with Tellurium."""
    import tellurium as te
    source = flatten_comp(sbml_path) if flatten else Path(sbml_path).read_text(encoding="utf-8")
    raw = te.sbmlToAntimony(source)
    text, removed, leftover = clean_antimony(raw)
    return AntimonyResult(text=text, raw=raw, removed=removed, leftover_is=leftover)
