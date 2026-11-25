"""Measure how much text a definition takes, in ways that cannot flatter it.

Four numbers per file (or per line range of a file):

  bytes       what is on disk
  lines       non-blank lines
  code        non-blank lines that are not comment-only or docstring (Python); for data files, the same as ``lines``
  gzip        the size after gzip -9: how much of the text is repetition a general-purpose compressor already removes
              (for Python, of the code lines only, so prose in comments and docstrings is not counted against the module)

``bytes`` and ``lines`` reward terseness of formatting (a pretty-printed XML has a line per tag); ``code`` is the honest
count for Python, which a docstring-rich module would otherwise inflate; ``gzip`` is the baseline any representation has to
beat to claim it removed *information* and not just repetition.
"""
from __future__ import annotations

import ast
import gzip
import io
import tokenize
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Size:
    bytes: int = 0
    lines: int = 0
    code: int = 0
    gzip: int = 0

    def __add__(self, other: "Size") -> "Size":
        return Size(self.bytes + other.bytes, self.lines + other.lines, self.code + other.code, self.gzip + other.gzip)

    def ratio(self, other: "Size", field: str = "bytes") -> float:
        """How many times bigger ``other`` is than ``self``."""
        return getattr(other, field) / max(getattr(self, field), 1)


def _gz(data: bytes) -> int:
    return len(gzip.compress(data, compresslevel=9, mtime=0))


def _code_lines(text: str) -> set[int]:
    """1-based numbers of the lines of Python source that hold code: not blank, not comment-only, not a docstring."""
    doc = set()
    try:
        for node in ast.walk(ast.parse(text)):
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                first = node.body[0] if node.body else None
                if isinstance(first, ast.Expr) and isinstance(getattr(first, "value", None), ast.Constant) \
                        and isinstance(first.value.value, str):
                    doc.update(range(first.lineno, first.end_lineno + 1))
    except SyntaxError:
        pass
    skip = {tokenize.COMMENT, tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT, tokenize.ENDMARKER}
    code = set()
    for tok in tokenize.generate_tokens(io.StringIO(text).readline):
        if tok.type not in skip:
            code.update(range(tok.start[0], tok.end[0] + 1))
    return code - doc


def measure(path: str | Path, first: int | None = None, last: int | None = None) -> Size:
    """Size of a file, or of lines ``first..last`` (1-based, inclusive) of it."""
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    all_lines = text.splitlines(keepends=True)
    lo, hi = (first or 1), (last or len(all_lines))
    chosen = all_lines[lo - 1:hi]
    data = "".join(chosen).encode("utf-8")
    nonblank = sum(1 for ln in chosen if ln.strip())
    if path.suffix == ".py":
        keep = sorted(n for n in _code_lines(text) if lo <= n <= hi)
        return Size(len(data), nonblank, len(keep), _gz("".join(all_lines[n - 1] for n in keep).encode("utf-8")))
    return Size(len(data), nonblank, nonblank, _gz(data))


def measure_text(text: str, python: bool = False) -> Size:
    """Size of a string (generated Antimony, exported SBML)."""
    data = text.encode("utf-8")
    nonblank = sum(1 for ln in text.splitlines() if ln.strip())
    if python:
        lines = text.splitlines(keepends=True)
        keep = sorted(_code_lines(text))
        return Size(len(data), nonblank, len(keep), _gz("".join(lines[n - 1] for n in keep).encode("utf-8")))
    return Size(len(data), nonblank, nonblank, _gz(data))


def class_range(path: str | Path, name: str) -> tuple[int, int]:
    """First and last line of ``class <name>`` in a Python file."""
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == name:
            return node.lineno, node.end_lineno
    raise KeyError(name)


def function_range(path: str | Path, name: str) -> tuple[int, int]:
    """First and last line of the top-level function ``name`` in a Python file."""
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node.lineno, node.end_lineno
    raise KeyError(name)


def total(sizes) -> Size:
    out = Size()
    for s in sizes:
        out = out + s
    return out
