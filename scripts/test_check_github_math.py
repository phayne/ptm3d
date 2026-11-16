#!/usr/bin/env python3
"""Self-test for the GitHub math checkers. No pytest, no deps.

    python3 scripts/test_check_github_math.py

Covers scripts/check_github_math.py (run in CI), tools/check_markdown_math.py and the
shared macro list tools/github_math_macros.{py,txt}. The cases that matter:

  * Both checkers passed docs/thermal-model.md while GitHub showed 15 of its equations
    broken: \\thickspace was printed in red and every \\tag{...} equation collapsed into a
    column of symbols. They only modeled the Markdown stage. The blocked-macro cases below
    pin that the renderer stage is checked too, in $...$, $$...$$, $`...`$ and ```math.
  * The fixed forms (operator alone, \\thinspace before a unit, \\qquad \\text{(B1)}) must
    pass, or the checkers would push authors back to the broken ones.
  * A macro not yet seen working on GitHub is a warning, not an error: exit status 0.
  * The allow list is derived from a saved GitHub page; --derive must report red nodes and
    labeled rows (mlabeledtr) as failures, not as allowed macros.
"""

import importlib.util
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cgm = load("cgm", os.path.join(HERE, "check_github_math.py"))
cmm = load("cmm", os.path.join(ROOT, "tools", "check_markdown_math.py"))
ghm = load("ghm", os.path.join(ROOT, "tools", "github_math_macros.py"))

_fail = 0


def check(cond, what):
    global _fail
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        _fail += 1


def write(tmp, name, text):
    path = os.path.join(tmp, name)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


BROKEN = {
    "thickspace in $$": "Text.\n\n$$a \\thickspace=\\thickspace b$$\n\nMore.\n",
    "tag in $$": "Text.\n\n$$x = y \\tag{B1}$$\n\nMore.\n",
    "thickspace inline": "Where $S_0 = 1361\\thickspace\\text{W m}^{-2}$ is used.\n",
    "thickspace in $`...`$": "Where $`a \\thickspace b`$ is used.\n",
    "tag in ```math": "Text.\n\n```math\nx = y \\tag{2}\n```\n\nMore.\n",
    # The first display equation of docs/thermal-model.md as it was on main.
    "thermal-model B1": (
        "$$\\rho(z)\\thinspace c_p(T)\\thinspace \\frac{\\partial T}{\\partial t} "
        "\\thickspace=\\thickspace \\frac{\\partial}{\\partial z}\\negthinspace\\left[\\thinspace "
        "k(T, z)\\thinspace \\frac{\\partial T}{\\partial z} \\thinspace\\right] \\tag{B1}$$\n"),
}
FIXED = {
    "operator alone, numbered": (
        "$$\\rho(z)\\thinspace c_p(T)\\thinspace \\frac{\\partial T}{\\partial t} = "
        "\\frac{\\partial}{\\partial z}\\negthinspace\\left[\\thinspace k(T, z)\\thinspace "
        "\\frac{\\partial T}{\\partial z} \\thinspace\\right] \\qquad \\text{(B1)}$$\n"),
    "thin space before unit": "Where $S_0 = 1361\\thinspace\\text{W m}^{-2}$ is used.\n",
    "math fence": "Text.\n\n```math\nx = \\frac{a}{b} \\qquad \\text{(2)}\n```\n",
}

with tempfile.TemporaryDirectory() as tmp:
    for label, text in BROKEN.items():
        path = write(tmp, "broken.md", text)
        findings, _ = cgm.scan(path)
        check(any(tok == "macro" for _, tok, _ in findings),
              f"check_github_math flags {label}")
        issues = cmm.check_markdown(path)
        check(any(rule == "macro" for _, rule, _ in issues),
              f"check_markdown_math flags {label}")
    for label, text in FIXED.items():
        path = write(tmp, "fixed.md", text)
        findings, warnings = cgm.scan(path)
        check(not findings and not warnings, f"check_github_math passes {label}")
        issues = cmm.check_markdown(path)
        check(not issues, f"check_markdown_math passes {label}")

    # Unverified macro: warning only, exit status 0, macro named in the output.
    path = write(tmp, "unverified.md", "Angle $\\varphi_0$ here.\n")
    for argv in ([sys.executable, os.path.join(HERE, "check_github_math.py"), path],
                 [sys.executable, os.path.join(ROOT, "tools", "check_markdown_math.py"), path]):
        r = subprocess.run(argv, capture_output=True, text=True)
        check(r.returncode == 0 and "\\varphi" in r.stdout and "unverified" in r.stdout,
              f"{os.path.basename(argv[1])}: unverified macro warns with exit 0")
    # Blocked macro: exit status 1.
    path = write(tmp, "blocked.md", BROKEN["tag in $$"])
    for argv in ([sys.executable, os.path.join(HERE, "check_github_math.py"), path],
                 [sys.executable, os.path.join(ROOT, "tools", "check_markdown_math.py"), path]):
        r = subprocess.run(argv, capture_output=True, text=True)
        check(r.returncode == 1, f"{os.path.basename(argv[1])}: blocked macro exits 1")

    # Derivation from a saved page: red node and labeled row are failures.
    saved = write(tmp, "saved.html", (
        '<math-renderer><math data-latex="a \\thickspace = b" display="block">'
        '<mi data-latex="a">a</mi><mtext mathcolor="red" data-latex="\\thickspace">'
        '\\thickspace</mtext><mo data-latex="=">=</mo></math></math-renderer>'
        '<math-renderer><math data-latex="x \\tag{1}"><mtable><mlabeledtr><mtd>'
        '<mtext data-latex="\\text{(}">(</mtext></mtd><mtd><mi data-latex="x">x</mi></mtd>'
        '</mlabeledtr></mtable></math></math-renderer>'
        '<math-renderer><math data-latex="\\alpha + \\beta"><mi data-latex="\\alpha">a</mi>'
        '</math></math-renderer>'))
    ok, red = ghm.derive(saved)
    check(red == {"\\thickspace", "\\tag"}, f"derive: red = {sorted(red)}")
    check("\\alpha" in ok and "\\beta" in ok and "\\tag" not in ok and "\\thickspace" not in ok,
          f"derive: allowed = {sorted(ok)}")

# The data file must keep the two measured failures blocked.
allow, block = ghm.load()
check({"\\thickspace", "\\tag"} <= set(block), "data file blocks \\thickspace and \\tag")
check(not (allow & set(block)), "no macro is both allowed and blocked")

print(f"\n{'FAILED' if _fail else 'passed'}: {_fail} failure(s)")
sys.exit(1 if _fail else 0)
