#!/usr/bin/env python3
"""Refuse the SHAPE of a landing-site identifier in the tracked tree.

PTM3D works on lunar landing-site candidates whose identifiers, coordinates and
region names are not published with the model. They are kept out of the tracked
tree, and ``src/python/tests/test_tracked_tree_hygiene.py`` is the check that
keeps them out -- against a deny list of the actual terms, which deliberately
lives OUTSIDE the repository (``~/.ptm3d/sensitive_identifiers.json``).

That design has one hole, and it is not hypothetical. **A machine with no local
list cannot run that check**, so the check SKIPS -- honestly, by name, but it
skips. Continuous-integration runners have no local list and never will, so the
one gate that would have caught the leak does not run in the place where merges
are decided. A pull request reintroduced 61 site-identifier hits across three
files in September 2026 and merged green for exactly this reason.

This file is the half of the check that CAN run anywhere, because it carries no
identifiers at all -- only their *shape*:

* the catalog identifier form: the three-character catalog prefix, one to three
  letters of region code, two to four digits;
* the site-filename form these meshes and matrices are named with: an
  underscore, a high-southern-latitude degree value to three decimals, ``S``,
  an underscore;
* the SHORT site form used in run names and in conversation: the two-letter
  ``MP`` prefix, a zero, two digits, standing alone (not inside a longer word
  or number). Added 2026-11-04 (PTM-820) after that form reached commit
  messages and code comments; it carries no location, but it is still a site
  identifier, and this project does not publish those.

The same shapes are also checked in COMMIT MESSAGES, which the tracked-tree
scan cannot see (PTM-820: the short form was found in seven commit messages
on main while every tracked file was clean):

* ``--commits RANGE`` scans the message of every commit in a git revision
  range (``BASE..HEAD``);
* ``--github-event`` derives that range from the GitHub Actions event (the PR
  range, or the push ``before..after``), so CI scans only NEW commits and the
  history already on main cannot fail it. A push that creates a ref (a
  release tag, a new branch) scans the commits it adds that are not already
  on the default branch -- none, for a tag on a commit already there -- and
  the message of an annotated tag; if the clone cannot say what is already
  on the default branch, it scans the tagged commit alone and warns. It never
  scans the whole history;
* ``--message-file FILE`` scans one message before it is committed, which is
  what ``.githooks/commit-msg`` runs.

Jupyter notebooks are read as notebooks, not as lines of text. A plotted figure
is stored in the notebook as base64-encoded image data, often hundreds of
kilobytes on a single line, and in that much base64 the short shapes above
occur by chance. So for a ``.ipynb`` file the check parses the JSON and scans
every piece of text in it -- cell sources, text and HTML outputs, tracebacks,
metadata, attachment names -- but not the encoded pixels of an ``image/png``,
``image/jpeg`` or ``image/gif`` entry. What it scans in an image instead is the
only place in one where readable text can live: its text metadata (PNG text
chunks, JPEG comment and application segments, GIF comment and application
blocks), after decoding. A notebook that is not valid JSON, or an image entry
that does not decode as the format it claims, is scanned as plain text, so
nothing is skipped because it could not be parsed. A notebook hit is reported
by its place in the notebook (cell, output, field, line) rather than a file
line number.

Neither pattern is or contains a real identifier, a real coordinate or a real
place name, so both are safe in a tracked file -- which is the entire reason
this check exists as a separate file rather than as more entries on the local
list. It is strictly WEAKER than the deny-list check: it cannot see a bare
region code used as a Python identifier, an unprefixed site number, or a
northern-hemisphere coordinate. It is not a replacement. It is the part that
runs in CI, and it is sized to catch the way these identifiers actually come
back -- pasted in whole, as a path or a fixture name.

Run the local deny-list check too, whenever the machine has the list::

    python3 src/python/tests/test_tracked_tree_hygiene.py

Usage::

    python3 scripts/check_site_id_patterns.py            # whole tracked tree
    python3 scripts/check_site_id_patterns.py PATH ...   # specific paths
    python3 scripts/check_site_id_patterns.py --commits origin/main..HEAD
    python3 scripts/check_site_id_patterns.py --github-event   # in CI
    python3 scripts/check_site_id_patterns.py --message-file .git/COMMIT_EDITMSG

Exits 0 clean, 1 on a hit, and 2 on a usage error OR when the tracked file list
(or the commit range) cannot be obtained. That last one matters: a checker that reports "clean"
because git failed is a checker that has learned to say yes, and this one
replaced a gate that was already failing by not running.

``scripts/test_check_site_id_patterns.py`` pins all of that, including the exit
codes and the fact that this file's own patterns do not match this file. Run it
after touching this one.
"""

from __future__ import annotations

import argparse
import base64
import binascii
import hashlib
import json
import os
import re
import subprocess
import sys
import zlib

# --- the patterns -----------------------------------------------------------
# Shapes only. Nothing here names a site, a coordinate or a region: read this
# table aloud and it tells you the FORM of the identifiers this project keeps
# private, not one of them. Anything that would name one belongs on the local
# deny list instead, which is not in this repository.
#
# Kept deliberately narrow. A pattern that also matched ordinary text would be
# switched off by the first person it inconvenienced, and a gate that gets
# switched off is worse than no gate, because the switching-off is quiet.
PATTERNS = [
    (
        r"A03[A-Z]{1,3}[0-9]{2,4}",
        "a landing-site catalog identifier (prefix + region code + number)",
    ),
    (
        r"_8[0-9]\.[0-9]{3}S_",
        "a site-mesh filename's high-southern-latitude coordinate field",
    ),
    # Bounded on both sides so it cannot fire inside a longer token: a word
    # ending in "...mp" followed by digits, or a run of base64 text (an
    # unbounded case-insensitive form matches base64 easily). An underscore is NOT a
    # boundary character here on purpose: a "<tag>_" or "_<tag>" filename
    # field is the likeliest way the short form is committed.
    (
        r"(?<![A-Za-z0-9])MP0[0-9]{2}(?![0-9])",
        "a short landing-site identifier (MP prefix + zero + two digits)",
    ),
]

#: Specific LINES exempt by content, keyed on ``(path, sha256 of the line with
#: surrounding whitespace stripped)``. Hashes, not text, so this table does not
#: itself carry the identifier, and content rather than line number, so an
#: entry survives edits elsewhere in the file and dies the moment its own line
#: is edited (an edited line is a new line and must not carry the form).
#:
#: The entries are grandfathered under Paul's ruling of 2026-11-04 (PTM-820):
#: the short form carries no location information, so lines that ALREADY
#: carried it when the check was added stay as they are; the check is for new
#: text. They are the five code comments that stack/r2026-11-v2 brings to
#: main (they are not on main when this table is written, so on main every
#: entry is reported as "not present" -- a notice, not a failure). Nothing may
#: be ADDED here except by the same kind of ruling: fix new text instead.
#: Remove an entry once its line is gone; the whole-tree scan names stale ones.
ALLOW_LINES: dict[tuple[str, str], str] = {
    ("src/cpp/illumination/illumination_model.cpp",
     "b32ff919fb00d01ebab744d0cc4c10af948024885753906a37aa5b236486b7d7"):
        "PTM-820 grandfathered: production-profile comment (store loop)",
    ("src/cpp/illumination/illumination_model.cpp",
     "84130230b715e395bbd1cc886f9e68fe4a23f1f89655d4fba52973006ec274c9"):
        "PTM-820 grandfathered: parallel-build FP-order comment",
    ("src/cpp/illumination/illumination_model.h",
     "d011b47272fc27e910e07f4d8c5fe9805dc06759cd44d2c2e5a63292ad388e07"):
        "PTM-820 grandfathered: per-frame array cost comment",
    ("src/cpp/illumination/illumination_model.h",
     "fc91269cb88ac8b21d3f4bcf6dd8cf74e38419a70059792bd690ef102ae6bef9"):
        "PTM-820 grandfathered: prep-state scan comment",
    ("src/cpp/main/ptm_thermal_sim.cpp",
     "ad9dbead8394d659fb6b6c7434a85d23fbd6c30bc953efb3c39b546b4516b759"):
        "PTM-820 grandfathered: single-threaded wall-time comment",
}


def line_key(path, line):
    """The ``ALLOW_LINES`` key for one line of one tracked file."""
    return (path.replace(os.sep, "/"),
            hashlib.sha256(line.strip().encode("utf-8")).hexdigest())

#: Files exempt by path, each with the reason it is exempt. It is EMPTY, and
#: keeping it empty is the goal: every entry is a hole in the gate.
#:
#: In particular this gate does not exempt itself, and does not need to. A
#: regex for the shape does not have the shape -- ``A03`` followed by ``[`` is
#: not ``A03`` followed by a letter -- so the pattern table above does not match
#: the pattern table above, and the self-test pins that rather than leaving it
#: to luck. The self-test likewise plants its shaped strings by assembling them
#: at run time, so it carries no literal either. Both files are scanned like
#: any other, which is what makes this check trustworthy on itself.
ALLOW_PATHS: dict[str, str] = {}

#: Never opened as text. A binary that happens to contain these bytes is not
#: the failure this guards against.
#:
#: This affects CONTENTS only -- the PATH of a binary is still checked, and has
#: to be, because a mesh or matrix named after its site is the single most
#: likely way one of these identifiers gets committed. ``.ply`` is in this set
#: and a site-named ``.ply`` is exactly the case.
#:
#: ``.tf`` and ``.tpc`` are deliberately ABSENT although the equivalent table in
#: ``src/python/tests/test_tracked_tree_hygiene.py`` has them: SPICE text
#: kernels are text (only ``.bsp``/``.bpc`` are binary), and a frame kernel
#: defining a site-fixed frame is a plausible place for one of these to land.
BINARY_SUFFIXES = {
    ".png", ".jpg", ".jpeg", ".gif", ".pdf", ".bsp", ".bpc",
    ".bin", ".fits", ".ply", ".vtk", ".obj", ".h5", ".hdf5", ".npy", ".npz",
    ".tif", ".tiff", ".zip", ".gz", ".ico", ".icns", ".woff", ".woff2",
}

#: End-of-line escape hatch, for the case this cannot anticipate: a genuine
#: coincidence in text that must stay. It has to be visible on the line, so a
#: reviewer sees the claim being made. There is deliberately no equivalent for
#: a PATH: a file whose NAME carries one of these shapes has to be renamed.
INLINE_ALLOW = "site-id-shape: allow"

#: ``re.IGNORECASE`` because the deny list this stands in for matches
#: case-insensitively (``test_tracked_tree_hygiene.py``), and the two halves
#: disagreeing on case would leave the gap open on exactly the axis where this
#: half is the only one that runs: a lowercased filename, a URL slug or a
#: ``str.lower()``-ed fixture key would ship past CI and be caught only on a
#: machine that has the local list.
COMPILED = [(re.compile(p, re.IGNORECASE), why) for p, why in PATTERNS]


class GitUnavailable(RuntimeError):
    """The tracked file list could not be obtained. Never a pass."""


def tracked_files(root="."):
    """Every tracked path, via whichever git works here.

    The CommandLineTools path first because on this project's Mac Pro
    ``/usr/bin/git`` is the xcrun shim, which exits 1 -- the same code it uses
    for "not a repository" -- so a caller that trusted the first git on PATH
    would read a broken git as an empty tree. On a CI runner that path does not
    exist and plain ``git`` is correct.

    Raises rather than returning an empty list: "no files" and "could not ask"
    are the same value and opposite meanings.
    """
    override = os.environ.get("PTM3D_GIT", "").strip()
    candidates = [override] if override else [
        "/Library/Developer/CommandLineTools/usr/bin/git", "git"]
    problems = []
    for git in candidates:
        try:
            out = subprocess.run([git, "-C", root, "ls-files", "-z"],
                                 capture_output=True, check=True)
        except (OSError, subprocess.CalledProcessError) as exc:
            problems.append(f"{git}: {exc}")
            continue
        names = [n for n in out.stdout.decode("utf-8", "replace").split("\0")
                 if n]
        if not names:
            problems.append(f"{git}: listed no tracked files")
            continue
        return names
    raise GitUnavailable(
        "could not list the tracked files, so this check has nothing to scan "
        "and will NOT report clean: " + "; ".join(problems))


#: ``check()`` uses line number 0 for a hit in the file's own PATH, which has
#: no line. Named rather than written as a bare 0 at three call sites.
PATH_LINE = 0


class Unreadable(RuntimeError):
    """A path was named but could not be read. Never a pass."""


def check(path, root=".", grandfathered=None):
    """``[(line number, why), ...]`` for one file. Line 0 means the PATH.

    For a Jupyter notebook the "line number" is a string naming where in the
    notebook the hit is (``notebook_text``), because nearly all of a notebook's
    bytes sit on a few very long lines of base64 image data.

    Reports the line and what SHAPE matched -- never the matching text. A gate
    against publishing an identifier must not print the identifier into a CI
    log, which is public to everyone who can read the build.

    **The path is checked first, and always.** A file NAMED after a site is the
    likeliest way one of these identifiers is committed -- a mesh, a matrix, a
    per-site config -- and the contents of such a file are usually innocent
    numbers. A first draft of this checker scanned contents only: a tracked
    mesh whose filename carried both a catalog identifier and a coordinate
    field reported ``clean``, doubly so because ``.ply`` skips the read. The
    path has no line to point at, so it is reported as line ``PATH_LINE`` (0).

    Raises ``Unreadable`` rather than returning ``[]`` when a named path cannot
    be opened. "I found nothing" and "I could not look" are the same value and
    opposite meanings, and this file exists because a check that could not run
    was mistaken for one that passed.
    """
    hits = []
    normalized = path.replace(os.sep, "/")
    for rx, why in COMPILED:
        if rx.search(normalized):
            hits.append((PATH_LINE, f"{why}, in the file's own NAME"))
            break

    full = os.path.join(root, path) if root != "." else path
    if os.path.splitext(path)[1].lower() in BINARY_SUFFIXES:
        return hits
    try:
        with open(full, "r", encoding="utf-8") as handle:
            text = handle.read()
    except UnicodeDecodeError:
        # Genuinely not text after all. The path was still checked above.
        return hits
    except OSError as exc:
        raise Unreadable(f"{path} could not be read ({exc})") from exc
    if os.path.splitext(path)[1].lower() == ".ipynb":
        pieces = notebook_text(text)
        if pieces is not None:
            for where, piece in pieces:
                for number, why in _check_lines(path, piece, grandfathered):
                    hits.append((f"{where} line {number}", why))
            return hits
        # Not valid JSON: scan it as plain text, line by line, as below.
    hits.extend(_check_lines(path, text, grandfathered))
    return hits


def _check_lines(path, text, grandfathered):
    """``[(line number, why), ...]`` for the lines of one piece of text."""
    hits = []
    for number, line in enumerate(text.splitlines(), start=1):
        if INLINE_ALLOW in line:
            continue
        for rx, why in COMPILED:
            if rx.search(line):
                key = line_key(path, line)
                if key in ALLOW_LINES:
                    # Exempt, but never silently: the caller counts these and
                    # the report says how many lines were let through.
                    if grandfathered is not None:
                        grandfathered.append(key)
                    break
                # One report per line. A site-shaped path on a line matches
                # BOTH patterns, and counting it twice made the summary say
                # "2 lines in 1 file" about a single line.
                hits.append((number, why))
                break
    return hits


# --- Jupyter notebooks ---------------------------------------------------------

#: Output and attachment entries whose value is a base64-encoded image. The
#: value is binary data, not text: scanned as text, hundreds of kilobytes of
#: base64 contain these short shapes by chance. Only the image's own text
#: metadata is scanned (see ``image_text``).
IMAGE_MIMES = {"image/png", "image/jpeg", "image/gif"}


def _join(value):
    """A notebook string field: one string, or a list of line strings."""
    if isinstance(value, list) and all(isinstance(v, str) for v in value):
        return "".join(value)
    return value if isinstance(value, str) else None


def notebook_text(text):
    """``[(where, text), ...]`` for every piece of text in a notebook.

    ``None`` if ``text`` is not a JSON notebook; the caller then scans the raw
    lines instead, so a broken notebook is never skipped.

    Every string is returned -- cell sources, text and HTML outputs, error
    tracebacks, metadata, attachment names -- and every dictionary KEY, except
    the VALUE of an ``IMAGE_MIMES`` entry. For those, the text metadata inside
    the image is returned instead. ``where`` is a path into the JSON, with
    zero-based indices (``cells[3].outputs[0].data["text/plain"]``), so a hit
    can be found without a line number in a file that is mostly one-line
    base64.
    """
    try:
        nb = json.loads(text)
    except ValueError:
        return None
    if not isinstance(nb, dict):
        return None
    pieces = []

    def walk(node, where):
        if isinstance(node, str):
            pieces.append((where, node))
        elif isinstance(node, list):
            for i, item in enumerate(node):
                walk(item, f"{where}[{i}]")
        elif isinstance(node, dict):
            for key, value in node.items():
                here = (f"{where}.{key}" if key.isidentifier()
                        else f"{where}[{json.dumps(key)}]")
                pieces.append((f"{here} (key)", key))
                payload = _join(value) if key in IMAGE_MIMES else None
                if payload is None:
                    walk(value, here)
                    continue
                found = image_text(key, payload)
                if found is None:
                    # Not decodable as the image it claims to be: scan it as
                    # text, exactly as before this parser existed.
                    pieces.append((here, payload))
                else:
                    for label, chunk in found:
                        pieces.append((f"{here} {label}", chunk))

    walk(nb, "notebook")
    return pieces


def image_text(mime, payload):
    """``[(label, text), ...]``: the human-readable metadata in one image.

    ``payload`` is base64. Returns ``None`` if it does not decode to the format
    ``mime`` names, so the caller falls back to scanning it as text.

    * PNG: the ``tEXt``, ``zTXt`` and ``iTXt`` chunks (keyword and text,
      decompressed where compressed);
    * JPEG: comment (``COM``) and application (``APPn``) segments, which carry
      EXIF and XMP text;
    * GIF: comment and application extension blocks.

    Pixel data is not text and is not returned.
    """
    try:
        data = base64.b64decode("".join(payload.split()), validate=True)
    except (binascii.Error, ValueError):
        return None
    try:
        if mime == "image/png":
            return _png_text(data)
        if mime == "image/jpeg":
            return _jpeg_text(data)
        return _gif_text(data)
    except (IndexError, ValueError, UnicodeDecodeError):
        return None


def _png_text(data):
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG")
    found, pos = [], 8
    while pos < len(data):
        if pos + 8 > len(data):
            raise ValueError("truncated PNG chunk header")
        length = int.from_bytes(data[pos:pos + 4], "big")
        kind = data[pos + 4:pos + 8]
        body = data[pos + 8:pos + 8 + length]
        if len(body) != length:
            raise ValueError("truncated PNG chunk")
        pos += 12 + length
        name = kind.decode("latin-1")
        if kind == b"tEXt":
            keyword, _, value = body.partition(b"\0")
            found.append((name, (keyword + b"\n" + value).decode("latin-1")))
        elif kind == b"zTXt":
            keyword, _, rest = body.partition(b"\0")
            value = zlib.decompress(rest[1:])
            found.append((name, (keyword + b"\n" + value).decode("latin-1")))
        elif kind == b"iTXt":
            keyword, _, rest = body.partition(b"\0")
            compressed = rest[0]
            language, _, rest = rest[2:].partition(b"\0")
            translated, _, value = rest.partition(b"\0")
            if compressed:
                value = zlib.decompress(value)
            found.append((name, "\n".join(
                part.decode("utf-8", "replace")
                for part in (keyword, language, translated, value))))
        elif kind == b"IEND":
            break
    return found


def _jpeg_text(data):
    if data[:2] != b"\xff\xd8":
        raise ValueError("not a JPEG")
    found, pos = [], 2
    while pos + 4 <= len(data):
        if data[pos] != 0xFF:
            raise ValueError("bad JPEG marker")
        marker = data[pos + 1]
        if marker == 0xDA:               # start of scan: pixel data follows
            break
        length = int.from_bytes(data[pos + 2:pos + 4], "big")
        if length < 2:
            raise ValueError("bad JPEG segment length")
        body = data[pos + 4:pos + 2 + length]
        if marker == 0xFE or 0xE0 <= marker <= 0xEF:
            label = "COM" if marker == 0xFE else f"APP{marker - 0xE0}"
            found.append((label, body.decode("latin-1")))
        pos += 2 + length
    return found


def _gif_text(data):
    if data[:6] not in (b"GIF87a", b"GIF89a"):
        raise ValueError("not a GIF")
    found = []
    flags = data[10]
    pos = 13 + (3 * 2 ** ((flags & 7) + 1) if flags & 0x80 else 0)

    def blocks(p):
        out = bytearray()
        while data[p]:
            out += data[p + 1:p + 1 + data[p]]
            p += 1 + data[p]
        return bytes(out), p + 1

    while pos < len(data):
        tag = data[pos]
        if tag == 0x3B:                  # trailer
            break
        if tag == 0x21:                  # extension
            label = data[pos + 1]
            body, pos = blocks(pos + 2)
            if label in (0xFE, 0xFF):    # comment, application
                found.append(("comment" if label == 0xFE else "application",
                              body.decode("latin-1")))
        elif tag == 0x2C:                # image descriptor
            flags = data[pos + 9]
            pos += 10 + (3 * 2 ** ((flags & 7) + 1) if flags & 0x80 else 0)
            _pixels, pos = blocks(pos + 1)   # skip the LZW code size byte
        else:
            raise ValueError("bad GIF block")
    return found


def redact(text):
    """``text`` with every matching span replaced by ``<site-id>``.

    For PRINTING a path: a hit in a file's NAME used to print the name whole,
    i.e. the identifier, into the CI log this check exists to keep it out of.
    The rest of the path is kept so the file can still be found.
    """
    for rx, _why in COMPILED:
        text = rx.sub("<site-id>", text)
    return text


def scan(paths, root=".", grandfathered=None):
    """``[(path, line number, why), ...]`` over the given tracked paths.

    ``grandfathered``, if a list, receives the ``ALLOW_LINES`` key of every
    line that matched and was exempted by that table.
    """
    hits = []
    for path in sorted(set(paths)):
        if path in ALLOW_PATHS:
            continue
        for number, why in check(path, root, grandfathered):
            hits.append((path, number, why))
    return hits


# --- commit messages (PTM-820) -----------------------------------------------

#: The hash git reports for "no commit": the ``before`` of a push that created
#: a ref. There is no range to scan from it.
NULL_SHA = "0" * 40


class RangeUnavailable(RuntimeError):
    """The commits to scan could not be listed. Never a pass."""


def _git_candidates():
    override = os.environ.get("PTM3D_GIT", "").strip()
    return [override] if override else [
        "/Library/Developer/CommandLineTools/usr/bin/git", "git"]


def _git(args, root="."):
    """stdout of the first git that runs ``args``; raises RangeUnavailable."""
    problems = []
    for git in _git_candidates():
        try:
            out = subprocess.run([git, "-C", root] + args,
                                 capture_output=True, check=True)
        except OSError as exc:
            problems.append(f"{git}: {exc}")
            continue
        except subprocess.CalledProcessError as exc:
            err = exc.stderr.decode("utf-8", "replace").strip()
            problems.append(f"{git}: exit {exc.returncode}: {err}")
            continue
        return out.stdout.decode("utf-8", "replace")
    raise RangeUnavailable("; ".join(problems))


def message_lines(text):
    """The lines git would actually store for a message being edited.

    Drops ``#`` comment lines and everything below a ``--verbose`` scissors
    line, as git's default cleanup does. Those lines name the branch and the
    changed files, so scanning them would make the hook refuse a commit whose
    file NAMES the tree check already refuses, and report the wrong thing.
    """
    kept = []
    for line in text.splitlines():
        if line.startswith("# ------------------------ >8"):
            break
        if line.startswith("#"):
            kept.append("")          # keep numbering aligned with the editor
            continue
        kept.append(line)
    return kept


def check_message(lines):
    """``[(line number, why), ...]`` for one commit message.

    There is no inline escape hatch for a message, unlike a file: a message is
    rewritten before it is merged, never annotated. Only shapes are reported.
    """
    hits = []
    for number, line in enumerate(lines, start=1):
        for rx, why in COMPILED:
            if rx.search(line):
                hits.append((number, why))
                break
    return hits


def scan_commits(rev_range, root="."):
    """``(n commits, [(sha, line number, why), ...])`` for a revision range.

    ``rev_range`` is anything ``git rev-list`` accepts as one argument
    (``A..B``, ``A^!`` for one commit) or a list of such arguments
    (``[TIP, "^BASE"]``). A bare sha is NOT one commit: ``git rev-list SHA``
    lists everything reachable from it. Merge commits are included: GitHub's
    "Merge pull request #N from owner/branch" message carries the branch NAME.
    """
    revs = list(rev_range) if isinstance(rev_range, (list, tuple)) \
        else [rev_range]
    try:
        shas = _git(["rev-list"] + revs + ["--"], root).split()
    except RangeUnavailable as exc:
        raise RangeUnavailable(
            f"could not list the commits in {rev_range!r}, so this check has "
            f"nothing to scan and will NOT report clean: {exc}") from exc
    hits = []
    for sha in shas:
        body = _git(["log", "-1", "--format=%B", sha], root)
        for number, why in check_message(body.splitlines()):
            hits.append((sha, number, why))
    return len(shas), hits


def _peel_to_commit(obj, root="."):
    """The commit sha an object (commit, annotated tag) points at."""
    try:
        return _git(["rev-parse", "--verify", "--quiet", f"{obj}^{{commit}}"],
                    root).strip()
    except RangeUnavailable as exc:
        raise RangeUnavailable(f"the pushed object {obj} does not resolve to "
                               f"a commit in this clone: {exc}") from exc


def _comparison_base(event, root="."):
    """``(base commit sha, ref name)`` that a ref-creating push is measured
    against, or ``(None, reason)`` when there is none to trust.

    The base is the repository's default branch as this clone's remote sees
    it (``refs/remotes/origin/<default>``, the default taken from the event
    payload, else ``main``). Only the default branch is used, not every other
    ref in the clone: a branch or tag pushed in the same ``git push`` is
    already in the clone when this runs, and excluding it would let that
    push's new commits go unscanned by both events.
    """
    default = ((event.get("repository") or {}).get("default_branch")
               or "main")
    ref = f"refs/heads/{default}"
    if event.get("ref") == ref:
        return None, (f"the pushed ref IS the default branch {default!r}, so "
                      f"it cannot be its own comparison base")
    try:
        shallow = _git(["rev-parse", "--is-shallow-repository"], root).strip()
    except RangeUnavailable as exc:
        return None, f"could not ask whether the clone is shallow ({exc})"
    if shallow == "true":
        return None, "the clone is shallow, so reachability is unknown"
    remote = f"refs/remotes/origin/{default}"
    try:
        base = _git(["rev-parse", "--verify", "--quiet",
                     f"{remote}^{{commit}}"], root).strip()
    except RangeUnavailable:
        return None, f"{remote} is not in this clone"
    return base, remote


def _new_ref_range(event, after, reason, root="."):
    """The range for a push with no usable ``before``: what is new to the
    default branch, or the peeled tip alone (warned) when that is unknown."""
    tip = _peel_to_commit(after, root)
    base, where = _comparison_base(event, root)
    if base is None:
        return f"{tip}^!", (f"WARNING: {reason} and {where}; only the pushed "
                            f"tip commit is scanned")
    return [tip, f"^{base}"], (f"the commits this push adds that are not "
                               f"already on {where} ({reason})")


def range_from_github_event(event_name, event, root="."):
    """``(revision range, explanation)`` for one GitHub Actions event.

    Only commits this event ADDS are scanned, so history already on the
    default branch (which carries the short form in a few older messages that
    stay as they are) never fails, and a genuinely new commit always is:

    * ``pull_request``: ``base.sha..head.sha``, the PR's own commits;
    * ``push``: ``before..after``;
    * a push that CREATED the ref (``before`` is all zeros: a new tag or
      branch), or whose ``before`` is not in the clone (a force push): the
      pushed object is peeled to its commit, and the commits reachable from
      it but not from the default branch are scanned. A tag on a commit that
      is already on the default branch therefore scans nothing, and says so.
      When the default branch cannot serve as the base (the clone is shallow,
      the remote ref is missing, or the pushed ref is the default branch
      itself), only the peeled tip commit is scanned, with a warning -- never
      the whole history;
    * anything else (``workflow_dispatch``): the checked-out commit alone,
      likewise stated.

    The message of an annotated tag is published text too; ``tag_messages``
    covers it, and ``--github-event`` scans it alongside the commits.
    """
    if event_name == "pull_request":
        pr = event.get("pull_request") or {}
        base = (pr.get("base") or {}).get("sha")
        head = (pr.get("head") or {}).get("sha")
        if not base or not head:
            raise RangeUnavailable("pull_request event without base/head sha")
        return f"{base}..{head}", "the commits of this pull request"
    if event_name == "push":
        before, after = event.get("before"), event.get("after")
        if not after or after == NULL_SHA:
            raise RangeUnavailable("push event without an 'after' sha "
                                   "(a deleted ref has nothing to scan)")
        if not before or before == NULL_SHA:
            return _new_ref_range(
                event, after, "this push created the ref, so there is no "
                "before..after range", root)
        try:
            _git(["cat-file", "-e", f"{before}^{{commit}}"], root)
        except RangeUnavailable:
            return _new_ref_range(
                event, after, "the push's 'before' commit is not in this "
                "clone (force push?)", root)
        return f"{before}..{after}", "the commits this push added"
    return "HEAD^!", (f"WARNING: event {event_name!r} has no commit range; "
                      f"only the checked-out commit is scanned")


def tag_messages(obj, root="."):
    """``[(tag sha, message lines), ...]`` for every annotated tag object
    met while peeling ``obj`` to a commit; empty for a commit or a
    lightweight tag. A signature block below the message is not scanned (it
    is base64, where short shapes occur by chance)."""
    found, seen = [], set()
    while obj not in seen:
        seen.add(obj)
        kind = _git(["cat-file", "-t", obj], root).strip()
        if kind != "tag":
            break
        raw = _git(["cat-file", "tag", obj], root)
        header, _, body = raw.partition("\n\n")
        lines = []
        for line in body.splitlines():
            if line.startswith("-----BEGIN "):
                break
            lines.append(line)
        found.append((obj, lines))
        target = [h for h in header.splitlines() if h.startswith("object ")]
        if not target:
            break
        obj = target[0].split()[1]
    return found


def _describe(rev_range):
    return " ".join(rev_range) if isinstance(rev_range, (list, tuple)) \
        else rev_range


def _report_message_hits(hits, label):
    for where, number, why in hits:
        print(f"{label(where)} message line {number}: matches {why}")
    print(f"\n{len(hits)} commit-message line(s) carry the SHAPE of a "
          f"landing-site identifier.\n"
          f"This project does not publish those, in a file or in a commit "
          f"message. Reword the message (git commit --amend, or an "
          f"interactive rebase for an older commit) before it is merged. "
          f"There is no escape hatch for a message.")


def _parser():
    p = argparse.ArgumentParser(
        prog="check_site_id_patterns.py",
        description="Refuse landing-site identifier SHAPES in the tracked tree.")
    p.add_argument("paths", nargs="*", metavar="PATH",
                   help="specific tracked paths; default is the whole tracked tree")
    p.add_argument("--root", default=".",
                   help="repository root to scan (default: the current directory)")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--commits", metavar="RANGE",
                      help="scan the commit MESSAGES in a git revision range "
                           "(e.g. origin/main..HEAD) instead of the tree")
    mode.add_argument("--github-event", action="store_true",
                      help="scan the commit messages of the GitHub Actions "
                           "event in GITHUB_EVENT_NAME/GITHUB_EVENT_PATH "
                           "(the PR range, or the push before..after)")
    mode.add_argument("--message-file", metavar="FILE",
                      help="scan one commit message file (the commit-msg hook)")
    return p


def _main_messages(args):
    if args.message_file:
        try:
            with open(args.message_file, "r", encoding="utf-8",
                      errors="replace") as handle:
                lines = message_lines(handle.read())
        except OSError as exc:
            print(f"check_site_id_patterns: {args.message_file} could not be "
                  f"read ({exc}); refusing to report clean", file=sys.stderr)
            return 2
        hits = [("", n, w) for n, w in check_message(lines)]
        if hits:
            _report_message_hits(hits, lambda _w: "commit message:")
            return 1
        return 0

    try:
        if args.github_event:
            name = os.environ.get("GITHUB_EVENT_NAME", "")
            path = os.environ.get("GITHUB_EVENT_PATH", "")
            if not name or not path:
                raise RangeUnavailable("GITHUB_EVENT_NAME / GITHUB_EVENT_PATH "
                                       "are not set (not a GitHub Actions run?)")
            try:
                with open(path, "r", encoding="utf-8") as handle:
                    event = json.load(handle)
            except (OSError, ValueError) as exc:
                raise RangeUnavailable(f"could not read the event payload "
                                       f"{path}: {exc}") from exc
            rev_range, why = range_from_github_event(name, event, args.root)
            print(f"check_site_id_patterns: event {name!r}: scanning {why} "
                  f"({_describe(rev_range)})")
            if why.startswith("WARNING"):
                print(f"::warning::{why}")
        else:
            rev_range = args.commits
        count, commit_hits = scan_commits(rev_range, args.root)
        hits = [(f"commit {sha[:12]}", n, w) for sha, n, w in commit_hits]
        tags = []
        if args.github_event and name == "push":
            tags = tag_messages(event["after"], args.root)
            for sha, lines in tags:
                hits.extend((f"tag {sha[:12]}", n, w)
                            for n, w in check_message(lines))
    except RangeUnavailable as exc:
        print(f"check_site_id_patterns: {exc}", file=sys.stderr)
        return 2
    if hits:
        _report_message_hits(hits, lambda where: where)
        return 1
    if count == 0 and isinstance(rev_range, list):
        print("check_site_id_patterns: this push adds no commits that are "
              "not already on the default branch; nothing new to scan")
    tag_note = (f" and {len(tags)} annotated tag message(s)" if tags else "")
    print(f"check_site_id_patterns: clean ({count} commit message(s) in "
          f"{_describe(rev_range)}{tag_note})")
    return 0


def main(argv):
    # STRICT, for the same reason check_american_english.py is: a hand-rolled
    # loop that files an unrecognized flag under "a path" reports clean for
    # arguments it never understood.
    args = _parser().parse_intermixed_args(argv[1:])
    if args.commits or args.github_event or args.message_file:
        if args.paths:
            print("check_site_id_patterns: PATH arguments cannot be combined "
                  "with a commit-message mode", file=sys.stderr)
            return 2
        return _main_messages(args)
    try:
        paths = args.paths or tracked_files(args.root)
    except GitUnavailable as exc:
        print(f"check_site_id_patterns: {exc}", file=sys.stderr)
        return 2
    grandfathered = []
    try:
        hits = scan(paths, args.root, grandfathered)
    except Unreadable as exc:
        print(f"check_site_id_patterns: {exc}; refusing to report clean on a "
              f"file this check could not open", file=sys.stderr)
        return 2
    for path, number, why in hits:
        shown = redact(path)
        where = f"{shown}:{number}" if number != PATH_LINE else f"{shown} (name)"
        print(f"{where}: matches {why}")
    if hits:
        files = len({p for p, _n, _w in hits})
        print(f"\n{len(hits)} line(s) in {files} file(s) carry the SHAPE of a "
              f"landing-site identifier.\n"
              f"This project does not publish those. Remove them, or -- for a "
              f"genuine coincidence on a LINE -- put '{INLINE_ALLOW}' on it. A "
              f"hit in a file's NAME has no escape hatch: rename the file.\n"
              f"For test fixtures, invent a tag that does not match these "
              f"shapes; for real local paths, register them in "
              f"~/.ptm3d/test_fixtures.json and look them up through "
              f"src/python/tests/local_fixtures.py.\n"
              f"This is the shape-only half of the gate. On a machine with the "
              f"local identifier list, also run "
              f"src/python/tests/test_tracked_tree_hygiene.py.")
        return 1
    print(f"check_site_id_patterns: clean ({len(paths)} tracked files)")
    _report_grandfathered(grandfathered, whole_tree=not args.paths)
    return 0


def _report_grandfathered(seen, whole_tree):
    """Say how many lines ALLOW_LINES let through, and which entries are stale.

    A notice, never a failure: the entries exist ahead of the lines they cover
    (see ``ALLOW_LINES``), so "not present" is expected on a tree the stack
    has not reached. Named so the table cannot outlive its lines unnoticed.
    """
    if seen:
        print(f"check_site_id_patterns: {len(seen)} matching line(s) exempted "
              f"by ALLOW_LINES (PTM-820 grandfathered; see the table)")
    if whole_tree:
        missing = sorted(set(ALLOW_LINES) - set(seen))
        if missing:
            print(f"check_site_id_patterns: notice: {len(missing)} ALLOW_LINES "
                  f"entr{'y is' if len(missing) == 1 else 'ies are'} not "
                  f"present in this tree; remove an entry once its line is "
                  f"gone for good:")
            for path, digest in missing:
                print(f"  {path} sha256:{digest[:12]}")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
