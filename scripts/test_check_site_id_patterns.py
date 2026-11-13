#!/usr/bin/env python3
"""Self-test for ``scripts/check_site_id_patterns.py``.

A hygiene gate that cannot fail is a hygiene gate nobody should trust, and this
project has already paid for that lesson twice: the deny-list check skips on a
machine with no local list (which is every CI runner), and a pull request
merged 61 site-identifier hits through it. So the shape-only check that
replaces it in CI gets its own gate, run as a CI step BEFORE the check itself,
exactly like ``scripts/test_check_american_english.py``.

What is pinned here:

* the checker finds a planted identifier shape, and does not find one in a
  clean file (so it can fail, and it does not fire on everything);
* every shape is caught, not just the first, including the short form
  (PTM-820) in a file, in a file NAME, and in lower case;
* the ``ALLOW_LINES`` table exempts exactly its lines -- same path and same
  content -- and nothing else, and its exemptions are counted, not silent;
* commit MESSAGES are checked (``--message-file``, ``--commits``,
  ``--github-event``), only new commits are scanned for an event, comment
  lines of a message being edited are ignored, and a range that cannot be
  listed is exit 2;
* the exit codes: 0 clean, 1 on a hit, **2 when git cannot answer** -- the last
  one is the whole point, because reporting clean on a broken git is how a
  gate learns to say yes;
* the inline escape hatch and the ALLOW_PATHS table both actually exempt;
* the checker does not print the matching text, only the line and the shape;
* the checker's own files (and the commit-msg hook) are clean under the
  checker, with ALLOW_PATHS empty and no ALLOW_LINES entry for them -- so the
  gate is not exempting itself to pass.

**Every shaped string in this file is assembled at run time from fragments.**
Written out whole, they would be literals in a tracked file, and the checker
would flag its own self-test -- correctly. Assembling them is not a trick to
get around the gate; it is the demonstration that the gate works on the file
that plants them.

Run directly::

    python3 scripts/test_check_site_id_patterns.py
"""

from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import check_site_id_patterns as C                            # noqa: E402

REPO_ROOT = HERE.parent

#: A catalog-shaped identifier, assembled so this file contains no literal of
#: that shape. It is not a real identifier: the region code and number are
#: invented, and the gate cannot tell the difference -- which is correct, since
#: shape is all it knows.
PLANTED_ID = "A03" + "ZZ" + "9999"

#: A site-filename-shaped coordinate field, likewise assembled. 89.999 S is not
#: a place this project works on; it is a shape.
PLANTED_COORD = "_8" + "9.999S" + "_"

#: A short-form identifier (PTM-820), likewise assembled. The number is
#: invented; the gate cannot tell, and should not.
PLANTED_SHORT = "MP0" + "99"

#: A tag of the kind the tests are supposed to use instead -- invented, and
#: deliberately NOT of any of the shapes.
SAFE_TAG = "TESTSITE01"


class Failure(AssertionError):
    pass


def _git():
    for git in ("/Library/Developer/CommandLineTools/usr/bin/git", "git"):
        try:
            subprocess.run([git, "--version"], capture_output=True, check=True)
            return git
        except (OSError, subprocess.CalledProcessError):
            continue
    return None


def _commit(root: Path, message: str, name: str = "f.txt") -> str:
    """Commit one file change with ``message``; returns the new sha."""
    git = _git()
    target = root / name
    target.write_text((target.read_text() if target.exists() else "") + "x\n")
    ident = ["-c", "user.name=t", "-c", "user.email=t@example.invalid",
             "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null"]
    subprocess.run([git, "-C", str(root)] + ident + ["add", "-A"],
                   capture_output=True, check=True)
    subprocess.run([git, "-C", str(root)] + ident +
                   ["commit", "-q", "-m", message],
                   capture_output=True, check=True)
    return subprocess.run([git, "-C", str(root), "rev-parse", "HEAD"],
                          capture_output=True, check=True,
                          text=True).stdout.strip()


def _repo(tmp: Path, files: dict) -> Path:
    """A throwaway git repository with ``files`` written and staged."""
    git = _git()
    assert git, "no working git binary to build the fixture repository with"
    subprocess.run([git, "init", "-q", str(tmp)], capture_output=True,
                   check=True)
    for name, body in files.items():
        target = tmp / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(body)
    subprocess.run([git, "-C", str(tmp), "add", "-A"], capture_output=True,
                   check=True)
    return tmp


def _run(argv, cwd=None):
    """``(exit code, stdout, stderr)`` from the checker's own main().

    ``SystemExit`` counts as an exit code, because to a CI step it IS one:
    argparse raises it for an unrecognized flag, and a harness that let it
    escape would report the strictness arm as a crash rather than as the
    exit 2 the shell actually sees.
    """
    out, err = io.StringIO(), io.StringIO()
    previous = os.getcwd()
    if cwd:
        os.chdir(cwd)
    try:
        with redirect_stdout(out), redirect_stderr(err):
            try:
                code = C.main(["check_site_id_patterns.py"] + argv)
            except SystemExit as exc:
                code = exc.code if isinstance(exc.code, int) else 1
    finally:
        os.chdir(previous)
    return code, out.getvalue(), err.getvalue()


# --------------------------------------------------------------------------

def check_a_planted_catalog_identifier_is_found():
    """The load-bearing arm: it can fail."""
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {
            "dirty.py": f'MESH = "{PLANTED_ID}_R3km.ply"\n',
            "clean.py": f'MESH = "{SAFE_TAG}_R3km.ply"\n',
        })
        hits = C.scan(C.tracked_files(str(root)), str(root))
        names = {p for p, _n, _w in hits}
        if "dirty.py" not in names:
            raise Failure("a planted catalog identifier was NOT found; the "
                          "checker cannot fail, so its passing says nothing")
        if "clean.py" in names:
            raise Failure("a file with an invented, non-matching tag was "
                          "flagged; the pattern is too broad to keep")


def check_a_planted_coordinate_field_is_found():
    """The second shape is live too, not decoration next to the first."""
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {
            "paths.md": f"the mesh is at core{PLANTED_COORD}R3km_linear.ply\n",
        })
        hits = C.scan(C.tracked_files(str(root)), str(root))
        if not any(p == "paths.md" for p, _n, _w in hits):
            raise Failure("the coordinate-field shape was not caught, so only "
                          "one of the two patterns is doing anything")


def check_the_exit_code_is_1_on_a_hit_and_0_when_clean():
    """A checker that prints hits and returns 0 cannot gate anything."""
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"dirty.py": f"# {PLANTED_ID}\n"})
        code, out, _err = _run(["--root", str(root)], cwd=str(root))
        if code != 1:
            raise Failure(f"expected exit 1 on a hit, got {code}")
        if "dirty.py:1" not in out:
            raise Failure(f"the hit did not name the file and line: {out!r}")

    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"clean.py": f"# {SAFE_TAG}\n"})
        code, out, _err = _run(["--root", str(root)], cwd=str(root))
        if code != 0:
            raise Failure(f"expected exit 0 on a clean tree, got {code}: {out}")
        if "clean" not in out:
            raise Failure(f"a clean tree did not say so: {out!r}")


def check_git_failure_exits_2_and_never_reports_clean():
    """The failure mode this whole file exists for.

    The gate it replaces reports a SKIP when it cannot run, which in CI is
    indistinguishable from a pass. This one must be unable to do that.
    """
    with tempfile.TemporaryDirectory() as td:
        # A directory that is not a repository: git exits non-zero.
        code, out, err = _run(["--root", td], cwd=td)
        if code != 2:
            raise Failure(
                f"expected exit 2 when the tracked list cannot be obtained, "
                f"got {code}; stdout={out!r}")
        if "clean" in out:
            raise Failure("it reported clean with nothing to scan")
        if "will NOT report clean" not in err:
            raise Failure(f"the reason was not stated on stderr: {err!r}")

    # And the same through the API, so a future caller cannot get [] instead.
    with tempfile.TemporaryDirectory() as td:
        try:
            C.tracked_files(td)
        except C.GitUnavailable:
            pass
        else:
            raise Failure("tracked_files returned a list for a non-repository "
                          "instead of raising; an empty list reads as clean")


def check_a_site_shaped_FILENAME_is_found():
    """The likeliest leak of all: a file named after the run it came from.

    A first draft scanned contents only, and a tracked mesh whose filename
    carried both a catalog identifier and a coordinate field reported clean --
    doubly so, because a .ply is never opened. The contents here are innocent
    on purpose: the name is the whole hit.
    """
    mesh = f"{PLANTED_ID}{PLANTED_COORD}R3km_linear.ply"
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {
            mesh: "not really a mesh, and deliberately clean inside\n",
            f"notes_{PLANTED_ID}.md": "also clean inside\n",
            "fine.md": "clean name, clean contents\n",
        })
        hits = C.scan(C.tracked_files(str(root)), str(root))
        names = {p for p, _n, _w in hits}
        if mesh not in names:
            raise Failure("a .ply NAMED with an identifier shape was not "
                          "caught; a binary suffix must skip the CONTENTS, "
                          "never the name")
        if f"notes_{PLANTED_ID}.md" not in names:
            raise Failure("a text file named with an identifier shape was not "
                          "caught by its name")
        if "fine.md" in names:
            raise Failure("a clean file was flagged")
        for path, number, why in hits:
            if number == C.PATH_LINE and "NAME" not in why:
                raise Failure(f"a name hit did not say so: {why!r}")


def check_the_patterns_are_case_insensitive():
    """The deny list matches case-insensitively; this half must agree.

    Where they disagree is exactly where this half is the only gate that runs,
    so a lowercased filename or a ``str.lower()``-ed fixture key would ship.
    """
    lowered = PLANTED_ID.lower() + PLANTED_COORD.lower() + "r3km.ply"
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"paths.md": f"see {lowered}\n"})
        hits = C.scan(C.tracked_files(str(root)), str(root))
        if not hits:
            raise Failure("a lowercased identifier shape was not caught; this "
                          "half and the deny list disagree on case")


def check_an_unreadable_named_path_exits_2():
    """Naming a path that cannot be read must not report clean.

    Explicit-path mode is how a pre-commit hook or a changed-files CI step
    would call this. A typo'd or renamed path used to be swallowed and counted
    as one clean file -- the same "learned to say yes" failure as a broken git,
    on a different road.
    """
    with tempfile.TemporaryDirectory() as td:
        code, out, err = _run(["--root", td, "does/not/exist.py"], cwd=td)
        if code != 2:
            raise Failure(f"expected exit 2 for an unreadable named path, "
                          f"got {code}; stdout={out!r}")
        if "clean" in out:
            raise Failure("it reported clean on a file it could not open")
        if "could not open" not in err and "could not be read" not in err:
            raise Failure(f"the reason was not stated: {err!r}")

    # And through the API, so a future caller cannot get [] instead.
    with tempfile.TemporaryDirectory() as td:
        try:
            C.check("does/not/exist.py", td)
        except C.Unreadable:
            pass
        else:
            raise Failure("check() returned a list for an unopenable path")


def check_one_line_with_both_shapes_is_reported_once():
    """The count has to mean what it says.

    A site-named path carries both shapes, which is the common case, and both
    patterns firing on it used to report the same line twice -- "2 lines in 1
    file" about one line.
    """
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {
            "notes.md": f"the mesh is {PLANTED_ID}{PLANTED_COORD}R3km.ply\n"})
        hits = [h for h in C.scan(C.tracked_files(str(root)), str(root))
                if h[1] != C.PATH_LINE]
        if len(hits) != 1:
            raise Failure(f"one line matching both shapes produced {len(hits)} "
                          f"reports: {hits}")


def check_spice_text_kernels_are_read_as_text():
    """.tf/.tpc are TEXT kernels, and a frame kernel can name a site frame.

    The deny-list check's own suffix table calls them binary; this one does
    not, and that difference is deliberate rather than an oversight.
    """
    for suffix in (".tf", ".tpc"):
        if suffix in C.BINARY_SUFFIXES:
            raise Failure(f"{suffix} is a SPICE TEXT kernel and must not be "
                          f"in BINARY_SUFFIXES")
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {
            "frames.tf": f"FRAME_{PLANTED_ID}_TOPO = -1000\n"})
        if not C.scan(C.tracked_files(str(root)), str(root)):
            raise Failure("an identifier shape inside a .tf text kernel was "
                          "not caught")


def check_the_inline_escape_hatch_exempts_one_line():
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {
            "note.md": f"a coincidence {PLANTED_ID}  <!-- {C.INLINE_ALLOW} -->\n"
                       f"and a real one {PLANTED_ID}\n",
        })
        hits = C.scan(C.tracked_files(str(root)), str(root))
        lines = sorted(n for _p, n, _w in hits)
        if lines != [2]:
            raise Failure(f"expected only line 2 to be reported, got {lines}")


def check_allow_paths_exempts_a_named_file():
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"vendor.py": f"# {PLANTED_ID}\n"})
        saved = dict(C.ALLOW_PATHS)
        try:
            C.ALLOW_PATHS["vendor.py"] = "test"
            hits = C.scan(C.tracked_files(str(root)), str(root))
        finally:
            C.ALLOW_PATHS.clear()
            C.ALLOW_PATHS.update(saved)
        if hits:
            raise Failure(f"ALLOW_PATHS did not exempt the file it names: {hits}")


def check_the_report_never_prints_the_matching_text():
    """A gate against publishing an identifier must not publish it in the log."""
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"dirty.py": f'SITE = "{PLANTED_ID}"\n',
                                f"mesh_{PLANTED_ID}.ply": "clean\n"})
        _code, out, err = _run(["--root", str(root)], cwd=str(root))
        if "<site-id>.ply (name)" not in out:
            raise Failure(f"a NAME hit was not reported with the identifier "
                          f"redacted: {out!r}")
        if PLANTED_ID in out or PLANTED_ID in err:
            raise Failure("the report quoted the matching text; a CI log is "
                          "readable by everyone who can read the build")


def check_the_checker_and_this_file_are_clean_with_no_exemptions():
    """The gate is not passing by exempting itself.

    ``ALLOW_PATHS`` is empty on purpose: a regex for a shape does not have that
    shape, and strings assembled at run time are not literals. If either file
    ever needs an exemption to pass, that is a design change and this arm is
    the one that will say so.
    """
    if C.ALLOW_PATHS:
        raise Failure(f"ALLOW_PATHS is not empty: {sorted(C.ALLOW_PATHS)}; "
                      f"every entry is a hole, and neither of the gate's own "
                      f"files should need one")
    for name in ("scripts/check_site_id_patterns.py",
                 "scripts/test_check_site_id_patterns.py",
                 ".githooks/commit-msg"):
        if any(path == name for path, _d in C.ALLOW_LINES):
            raise Failure(f"ALLOW_LINES exempts a line of {name}")
        hits = C.check(name, str(REPO_ROOT))
        if hits:
            raise Failure(f"{name} matches its own patterns at lines "
                          f"{[n for n, _w in hits]}")


def check_ordinary_text_does_not_match():
    """The patterns are narrow enough to survive contact with the repository.

    A gate that fires on normal prose gets switched off, and the switching-off
    is quiet. These are the near-misses worth pinning: a version-like token, a
    northern coordinate field, and the catalog prefix followed by punctuation
    rather than a region code.
    """
    benign = [
        "A03 is the catalog, not an identifier",
        "the prefix A03-something appears in prose",
        # Northern hemisphere, and so not the shape. The degree value is
        # invented: writing a REAL one here would pass this shape-only gate and
        # fail the deny-list one, which is what happened on the first draft.
        "core_89.999N_R3km_linear.ply",
        "a duration of 8.5 days",
        "version A0300 of the format",       # digits with no region code
        f"{SAFE_TAG}_00.000S_000.000E_R3km_linear.ply",
    ]
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"prose.md": "\n".join(benign) + "\n"})
        hits = C.scan(C.tracked_files(str(root)), str(root))
        if hits:
            raise Failure(f"ordinary text was flagged at lines "
                          f"{[n for _p, n, _w in hits]}: too broad")


def check_the_short_form_is_found_in_contents_names_and_lower_case():
    """PTM-820: the short form, in each of the places it actually lands."""
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {
            "comment.cpp": f"// profile of the {PLANTED_SHORT} window\n",
            "lower.md": f"see runs/{PLANTED_SHORT.lower()}_window3/\n",
            f"configs/{PLANTED_SHORT}_w6h.yaml": "clean: inside\n",
            "fine.md": f"{SAFE_TAG} only\n",
        })
        hits = C.scan(C.tracked_files(str(root)), str(root))
        names = {p for p, _n, _w in hits}
        for want in ("comment.cpp", "lower.md",
                     f"configs/{PLANTED_SHORT}_w6h.yaml"):
            if want not in names:
                raise Failure(f"the short form was not caught in {want}")
        if "fine.md" in names:
            raise Failure("a clean file was flagged")
        code, out, err = _run(["--root", str(root)], cwd=str(root))
        if code != 1:
            raise Failure(f"expected exit 1 on a short-form hit, got {code}")
        if PLANTED_SHORT in out + err or PLANTED_SHORT.lower() in out + err:
            raise Failure("the report quoted the short-form text")


def check_the_short_form_has_no_false_positives_on_near_misses():
    """Bounded on both sides, so a longer word or number is not a hit."""
    benign = [
        "TE" + PLANTED_SHORT + " is a longer word",      # letter before it
        "a" + PLANTED_SHORT.lower() + " token",           # letter before it
        PLANTED_SHORT + "1 is a five-digit number",       # digit after it
        "an MP3 file, and MP04 has one digit too few",
        "the MP prefix alone",
    ]
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"prose.md": "\n".join(benign) + "\n"})
        hits = C.scan(C.tracked_files(str(root)), str(root))
        if hits:
            raise Failure(f"near-misses of the short form were flagged at "
                          f"lines {[n for _p, n, _w in hits]}: too broad")


def check_allow_lines_exempts_exactly_its_lines_and_counts_them():
    """Same path AND same content, or it is a hit. Never silent."""
    line = f"    // measured on {PLANTED_SHORT} (2026-11-04)"
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {
            "src/a.cpp": f"int x;\n{line}\n",
            "src/b.cpp": f"{line}\n",                       # other path
            "src/c.cpp": f"{line} edited\n",                # edited line
        })
        saved = dict(C.ALLOW_LINES)
        try:
            C.ALLOW_LINES.clear()
            C.ALLOW_LINES[C.line_key("src/a.cpp", line)] = "test"
            seen = []
            hits = C.scan(C.tracked_files(str(root)), str(root), seen)
            code, out, _err = _run(["--root", str(root), "src/a.cpp"],
                                   cwd=str(root))
        finally:
            C.ALLOW_LINES.clear()
            C.ALLOW_LINES.update(saved)
        names = sorted({p for p, _n, _w in hits})
        if names != ["src/b.cpp", "src/c.cpp"]:
            raise Failure(f"ALLOW_LINES exempted the wrong lines: hits in "
                          f"{names}, expected only b.cpp and c.cpp")
        if len(seen) != 1:
            raise Failure(f"the exemption was not counted: {seen}")
        if code != 0 or "exempted by ALLOW_LINES" not in out:
            raise Failure(f"an exempted line passed silently or failed: "
                          f"exit {code}, {out!r}")


def check_allow_lines_entries_are_the_grandfathered_ones_only():
    """The table is the PTM-820 ruling and nothing else."""
    for (path, digest), why in C.ALLOW_LINES.items():
        if "PTM-820" not in why:
            raise Failure(f"ALLOW_LINES entry for {path} does not cite the "
                          f"ruling it rests on: {why!r}")
        if len(digest) != 64 or not path.startswith("src/cpp/"):
            raise Failure(f"malformed or unexpected ALLOW_LINES entry {path}")
    if len(C.ALLOW_LINES) > 5:
        raise Failure(f"ALLOW_LINES grew to {len(C.ALLOW_LINES)} entries; the "
                      f"ruling grandfathered five existing lines, not new ones")


def check_a_commit_message_file_is_checked():
    """The commit-msg hook's mode: a synthetic message with the short form."""
    with tempfile.TemporaryDirectory() as td:
        dirty = Path(td) / "dirty.txt"
        dirty.write_text(f"perf: faster store loop\n\nMeasured on the "
                         f"{PLANTED_SHORT} production window.\n")
        clean = Path(td) / "clean.txt"
        # The comment block git appends names the branch and the files; a
        # site-shaped NAME there is the tree check's job, not the message's.
        clean.write_text(f"perf: faster store loop\n\nMeasured on a "
                         f"production window.\n# On branch {PLANTED_SHORT}\n"
                         f"# ------------------------ >8 ------------------\n"
                         f"diff --git a/{PLANTED_ID}.md b/{PLANTED_ID}.md\n")
        code, out, err = _run(["--message-file", str(dirty)])
        if code != 1:
            raise Failure(f"a message with the short form exited {code}")
        if "message line 3" not in out:
            raise Failure(f"the hit did not name the message line: {out!r}")
        if PLANTED_SHORT in out + err:
            raise Failure("the report quoted the message text")
        code, out, _err = _run(["--message-file", str(clean)])
        if code != 0:
            raise Failure(f"a clean message (with comment lines) exited "
                          f"{code}: {out!r}")
        code, _out, _err = _run(["--message-file", str(Path(td) / "nope")])
        if code != 2:
            raise Failure(f"an unreadable message file exited {code}, not 2")


def check_a_commit_range_is_checked_and_only_that_range():
    """Old history outside the range can never fail; the range always can."""
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"README.md": "x\n"})
        old = _commit(root, f"old: carries {PLANTED_SHORT} (already on main)")
        base = _commit(root, "clean: a later commit")
        dirty = _commit(root, f"perf: measured on {PLANTED_SHORT.lower()}")
        head = _commit(root, f"Merge pull request #1 from o/{PLANTED_ID}-x")
        n, hits = C.scan_commits(f"{base}..{head}", str(root))
        if n != 2 or sorted(h[0] for h in hits) != sorted([dirty, head]):
            raise Failure(f"range scan: {n} commits, hits {hits}")
        n, hits = C.scan_commits(f"{old}..{base}", str(root))
        if n != 1 or hits:
            raise Failure(f"a clean range (old history excluded) reported "
                          f"{hits}")
        code, out, err = _run(["--root", str(root), "--commits",
                               f"{base}..{dirty}"], cwd=str(root))
        if code != 1 or PLANTED_SHORT.lower() in out + err:
            raise Failure(f"--commits on a dirty range: exit {code}, {out!r}")
        code, out, _err = _run(["--root", str(root), "--commits",
                                f"{dirty}..{dirty}"], cwd=str(root))
        if code != 0 or "clean (0 commit" not in out:
            raise Failure(f"an empty range: exit {code}, {out!r}")
        code, out, err = _run(["--root", str(root), "--commits",
                               "no-such-ref..HEAD"], cwd=str(root))
        if code != 2 or "clean" in out:
            raise Failure(f"an unresolvable range exited {code}, not 2")


def check_the_github_event_range_is_new_commits_only():
    """PR -> base..head; push -> before..after; no range -> tip, warned."""
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"README.md": "x\n"})
        old = _commit(root, f"old: carries {PLANTED_SHORT}")
        base = _commit(root, "clean")
        head = _commit(root, "clean too")
        pr = {"pull_request": {"base": {"sha": base}, "head": {"sha": head}}}
        if C.range_from_github_event("pull_request", pr, str(root))[0] != \
                f"{base}..{head}":
            raise Failure("pull_request range is not base..head")
        rng, why = C.range_from_github_event(
            "push", {"before": base, "after": head}, str(root))
        if rng != f"{base}..{head}" or "WARNING" in why:
            raise Failure(f"push range: {rng} ({why})")
        for before in (C.NULL_SHA, "1" * 40):
            rng, why = C.range_from_github_event(
                "push", {"before": before, "after": head}, str(root))
            if rng != head or not why.startswith("WARNING"):
                raise Failure(f"a push with no usable before: {rng} ({why})")
        rng, why = C.range_from_github_event("workflow_dispatch", {}, str(root))
        if rng != "HEAD" or not why.startswith("WARNING"):
            raise Failure(f"dispatch: {rng} ({why})")

        event = Path(td) / "event.json"
        saved = {k: os.environ.get(k) for k in
                 ("GITHUB_EVENT_NAME", "GITHUB_EVENT_PATH")}
        try:
            os.environ["GITHUB_EVENT_NAME"] = "push"
            os.environ["GITHUB_EVENT_PATH"] = str(event)
            event.write_text(json.dumps({"before": base, "after": head}))
            code, out, _err = _run(["--root", str(root), "--github-event"],
                                   cwd=str(root))
            if code != 0:
                raise Failure(f"a clean push range with dirty OLD history "
                              f"exited {code}: {out!r}")
            event.write_text(json.dumps({"before": old + "~1",
                                         "after": head}))
            code, out, _err = _run(["--root", str(root), "--github-event"],
                                   cwd=str(root))
            if code != 1:
                raise Failure(f"a push range containing the dirty commit "
                              f"exited {code}")
            event.write_text("{not json")
            code, _out, _err = _run(["--root", str(root), "--github-event"],
                                    cwd=str(root))
            if code != 2:
                raise Failure(f"an unreadable event payload exited {code}")
            os.environ.pop("GITHUB_EVENT_PATH")
            code, _out, _err = _run(["--root", str(root), "--github-event"],
                                    cwd=str(root))
            if code != 2:
                raise Failure(f"no event payload exited {code}, not 2")
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v


def check_the_commit_msg_hook_refuses_and_accepts():
    """The hook end to end, through a real ``git commit``."""
    hook = REPO_ROOT / ".githooks" / "commit-msg"
    if not os.access(hook, os.X_OK):
        raise Failure(f"{hook} is missing or not executable")
    git = _git()
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"README.md": "x\n"})
        (root / "scripts").mkdir()
        for name in ("check_site_id_patterns.py",):
            (root / "scripts" / name).write_text(
                (HERE / name).read_text())
        (root / ".githooks").mkdir()
        (root / ".githooks" / "commit-msg").write_text(hook.read_text())
        os.chmod(root / ".githooks" / "commit-msg", 0o755)
        ident = ["-c", "user.name=t", "-c", "user.email=t@example.invalid",
                 "-c", "commit.gpgsign=false", "-c", "core.hooksPath=.githooks"]
        subprocess.run([git, "-C", str(root), "add", "-A"], check=True,
                       capture_output=True)
        bad = subprocess.run([git, "-C", str(root)] + ident +
                             ["commit", "-q", "-m",
                              f"perf: measured on {PLANTED_SHORT}"],
                             capture_output=True, text=True)
        if bad.returncode == 0:
            raise Failure("the hook let a short-form message through")
        good = subprocess.run([git, "-C", str(root)] + ident +
                              ["commit", "-q", "-m", "perf: measured"],
                              capture_output=True, text=True)
        if good.returncode != 0:
            raise Failure(f"the hook refused a clean message: {good.stdout}"
                          f"{good.stderr}")


def _png(chunks):
    """A minimal valid 1x1 PNG with extra ``(type, data)`` chunks after IHDR."""
    import struct
    import zlib

    def chunk(kind, data):
        return (struct.pack(">I", len(data)) + kind + data +
                struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF))

    out = b"\x89PNG\r\n\x1a\n" + chunk(
        b"IHDR", struct.pack(">IIBBBBB", 1, 1, 8, 0, 0, 0, 0))
    for kind, data in chunks:
        out += chunk(kind, data)
    return (out + chunk(b"IDAT", zlib.compress(b"\x00\x00")) +
            chunk(b"IEND", b""))


def _png_with_shape_in_base64():
    """Base64 of a valid PNG whose ENCODED text contains a planted shape.

    The shape is base64 characters, so it is placed in the encoding directly:
    the bytes it decodes to go in a private chunk, aligned so the encoder
    reproduces it. The decoded image carries no text chunk at all.
    """
    import base64
    encoded = PLANTED_ID + "AAA"                     # 12 base64 characters
    raw = base64.b64decode(encoded)
    prefix = _png([])[:33]                           # signature + IHDR
    # Chunk header (8 bytes) + padding must bring the raw bytes to a multiple
    # of 3 bytes from the start of the file, where base64 groups begin.
    pad = (-(len(prefix) + 8)) % 3
    png = _png([(b"prVt", b"\x00" * pad + raw)])
    text = base64.b64encode(png).decode("ascii")
    assert PLANTED_ID in text, "fixture: the shape did not land in the base64"
    return text


def _notebook(cells):
    return json.dumps({"cells": cells, "metadata": {}, "nbformat": 4,
                       "nbformat_minor": 5}, indent=1)


def _code_cell(source="x = 1\n", outputs=()):
    return {"cell_type": "code", "execution_count": 1, "metadata": {},
            "source": [source], "outputs": list(outputs)}


def _image_output(b64, mime="image/png"):
    return {"output_type": "display_data", "metadata": {},
            "data": {mime: b64, "text/plain": ["<Figure>"]}}


def check_a_notebook_image_payload_is_not_scanned_as_text():
    """The base64 of a plotted figure is data, not text.

    The same payload scanned as raw lines IS flagged -- that is what failed the
    tutorial notebooks -- so this arm is discriminating, not vacuous.
    """
    b64 = _png_with_shape_in_base64()
    body = _notebook([_code_cell(outputs=[_image_output(b64)])])
    if not any(rx.search(line) for rx, _why in C.COMPILED
               for line in body.splitlines()):
        raise Failure("fixture: the raw lines do not match, so this arm "
                      "would pass whether or not images are skipped")
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), {"nb.ipynb": body,
                                "attached.ipynb": _notebook([{
                                    "cell_type": "markdown", "metadata": {},
                                    "source": ["![x](attachment:fig.png)"],
                                    "attachments": {"fig.png": {
                                        "image/png": [b64]}}}])})
        code, out, _err = _run(["--root", str(root)], cwd=str(root))
        if code != 0:
            raise Failure(f"a shape occurring by chance inside base64 image "
                          f"data was flagged (exit {code}): {out!r}")


def check_a_notebook_shape_in_text_is_still_found():
    """Source, text outputs, PNG text chunks: each still fails, by place."""
    import base64
    import zlib
    clean_png = base64.b64encode(_png([])).decode("ascii")

    def png_b64(chunks):
        return base64.b64encode(_png(chunks)).decode("ascii")

    shape = PLANTED_ID.encode("ascii")
    cases = {
        "source.ipynb": [_code_cell(source=f"mesh = '{PLANTED_ID}.ply'\n")],
        "markdown.ipynb": [{"cell_type": "markdown", "metadata": {},
                            "source": ["intro\n", f"see {PLANTED_SHORT}\n"]}],
        "stream.ipynb": [_code_cell(outputs=[{
            "output_type": "stream", "name": "stdout",
            "text": ["loading\n", f"loaded {PLANTED_ID}\n"]}])],
        "plain.ipynb": [_code_cell(outputs=[{
            "output_type": "execute_result", "execution_count": 1,
            "metadata": {}, "data": {"text/plain": [f"'{PLANTED_ID}'"]}}])],
        "tEXt.ipynb": [_code_cell(outputs=[_image_output(
            png_b64([(b"tEXt", b"Title\0" + shape)]))])],
        "zTXt.ipynb": [_code_cell(outputs=[_image_output(
            png_b64([(b"zTXt", b"Title\0\0" + zlib.compress(shape))]))])],
        "iTXt.ipynb": [_code_cell(outputs=[_image_output(
            png_b64([(b"iTXt", b"Title\0\0\0en\0\0" + shape)]))])],
        "notpng.ipynb": [_code_cell(outputs=[_image_output(
            base64.b64encode(b"not an image " + shape).decode("ascii")
            + " " + PLANTED_ID)])],
        "clean.ipynb": [_code_cell(outputs=[_image_output(clean_png)])],
    }
    files = {name: _notebook(cells) for name, cells in cases.items()}
    files["broken.ipynb"] = "{ not json\n" + PLANTED_ID + "\n"
    with tempfile.TemporaryDirectory() as td:
        root = _repo(Path(td), files)
        hits = C.scan(C.tracked_files(str(root)), str(root))
        found = {p: n for p, n, _w in hits}
        for name in files:
            if name == "clean.ipynb":
                if name in found:
                    raise Failure(f"a clean notebook was flagged at "
                                  f"{found[name]}")
            elif name not in found:
                raise Failure(f"a shape in {name} was not caught")
        where = {"source.ipynb": "cells[0].source",
                 "stream.ipynb": "cells[0].outputs[0].text",
                 "tEXt.ipynb": '"image/png"] tEXt',
                 "iTXt.ipynb": '"image/png"] iTXt',
                 "broken.ipynb": 2}
        for name, want in where.items():
            got = found.get(name)
            if (got != want) if isinstance(want, int) else (want not in str(got)):
                raise Failure(f"{name}: hit reported at {got!r}, expected "
                              f"{want!r}")
        code, out, err = _run(["--root", str(root)], cwd=str(root))
        if code != 1:
            raise Failure(f"notebook hits exited {code}, not 1")
        if PLANTED_ID in out + err or PLANTED_SHORT in out + err:
            raise Failure("the report quoted the matching notebook text")


def check_a_usage_error_exits_2():
    code, _out, _err = _run(["--nonsense"])
    if code != 2:
        raise Failure(f"expected exit 2 on an unrecognized flag, got {code}")
    code, _out, _err = _run(["--commits", "A..B", "--message-file", "m"])
    if code != 2:
        raise Failure(f"two message modes at once exited {code}, not 2")
    code, _out, _err = _run(["--commits", "A..B", "some/path.py"])
    if code != 2:
        raise Failure(f"a PATH with a message mode exited {code}, not 2")


CHECKS = [
    check_a_planted_catalog_identifier_is_found,
    check_a_planted_coordinate_field_is_found,
    check_a_site_shaped_FILENAME_is_found,
    check_the_patterns_are_case_insensitive,
    check_the_exit_code_is_1_on_a_hit_and_0_when_clean,
    check_git_failure_exits_2_and_never_reports_clean,
    check_an_unreadable_named_path_exits_2,
    check_one_line_with_both_shapes_is_reported_once,
    check_spice_text_kernels_are_read_as_text,
    check_the_inline_escape_hatch_exempts_one_line,
    check_allow_paths_exempts_a_named_file,
    check_the_report_never_prints_the_matching_text,
    check_the_checker_and_this_file_are_clean_with_no_exemptions,
    check_ordinary_text_does_not_match,
    check_the_short_form_is_found_in_contents_names_and_lower_case,
    check_the_short_form_has_no_false_positives_on_near_misses,
    check_allow_lines_exempts_exactly_its_lines_and_counts_them,
    check_allow_lines_entries_are_the_grandfathered_ones_only,
    check_a_commit_message_file_is_checked,
    check_a_commit_range_is_checked_and_only_that_range,
    check_the_github_event_range_is_new_commits_only,
    check_the_commit_msg_hook_refuses_and_accepts,
    check_a_notebook_image_payload_is_not_scanned_as_text,
    check_a_notebook_shape_in_text_is_still_found,
    check_a_usage_error_exits_2,
]


def main() -> int:
    failed = 0
    for fn in CHECKS:
        try:
            fn()
        except SystemExit as exc:
            # argparse exits on --nonsense; _run catches the code, so a
            # SystemExit escaping here is a real defect, not the usage arm.
            failed += 1
            print(f"FAIL {fn.__name__}: escaped SystemExit({exc.code})")
        except AssertionError as exc:
            failed += 1
            print(f"FAIL {fn.__name__}:\n  {exc}")
        else:
            print(f"PASS {fn.__name__}")
    print(f"\n[test_check_site_id_patterns] {len(CHECKS) - failed}/"
          f"{len(CHECKS)} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
