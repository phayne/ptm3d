# Contributing to PTM3D

## Scope of the repository

This repository contains the PTM3D distribution, which is the model, its tools, the GUI, the
Python package, the documentation, tutorials and examples, and the tests that check them. A
file belongs here if a user needs it to build, run, test, learn or extend PTM3D for their own
problem, or if a test in this repository runs it. A file does not belong here if it serves one
study, campaign, site, mission or paper. Examples are the scripts that run a particular
campaign, configurations and data for a particular site or study, scripts for the figures of a
paper, and diagnostics written for a single use. These are kept in a separate research
repository that depends on this one.

An example configuration must run with data that are in the repository or that
`scripts/fetch_data.py` or `scripts/fetch_kernels.py` provides. A tool that was written for a
study may be added once its study-specific defaults are removed and it has a test. If you are
unsure whether a file belongs here, ask before adding it.

Landing-site identifiers and coordinates must not appear in this repository, in any file,
name, image or notebook output, or in a commit message. Continuous integration checks the
tracked files and the message of every new commit for them
(`scripts/check_site_id_patterns.py`). To make the same check when you commit, install the
optional hook once per clone with `git config core.hooksPath .githooks`.

This repository was started on 2026-10-29 from the development repository, which is archived.
Commit ids and references to "the archived repository" in comments and documentation refer to
that archive.

The rest of this page describes three conventions.

## American English

PTM3D is written in American English. This applies to documentation, tutorials, code
comments, docstrings, and anything a user can see, including GUI labels and command-line
help.
<!-- american-english: allow-begin (this section has to quote what it forbids) -->
`color`, not `colour`; `behavior`, not `behaviour`; `normalize`, not
`normalise`; `meter`, not `metre`; `gray`, not `grey`.
<!-- american-english: allow-end -->

```bash
python3 scripts/check_american_english.py
python3 scripts/test_check_american_english.py   # after editing the checker
```

Continuous integration runs both on every push. Notebooks are parsed as JSON, so the spelling
in cells and in printed output is checked, and base64 image data is not.

Three kinds of British spelling are correct and are kept.

- Quoted text you do not own: a paper title in a citation, third-party
  license text, vendored code.
<!-- american-english: allow-begin -->
- External format literals: the OGC WKT unit name is spelled `metre`, and
  GDAL parses that string.
<!-- american-english: allow-end -->
- Code identifiers: a function, attribute, config key or saved-data field
  whose British name cannot change without breaking callers or existing runs.

Record such a case in the `ALLOW` table at the top of
`scripts/check_american_english.py`, which pairs a path with the substring that
makes the line acceptable, or put `american-english: allow` at the end of the
line. Several consecutive lines, such as prose that quotes the spellings it forbids, go
between `american-english: allow-begin` and `american-english: allow-end` in a comment, as in
this section. Do not make a pattern of the checker more general in order to pass one line.

Two words are already correct and should not be changed. `cancellation` keeps both L's in
American usage (only `canceled` and `canceling` drop one), and `analyses` is the plural of
`analysis`.

## Changing an engine configuration key

The configuration keys and their defaults are defined by the YAML parser of the model
(`src/cpp/tools/config_parser.cpp`, with the struct defaults in `config_parser.h`). The Python
package keeps a copy of that table in `src/python/ptm/engine_defaults.py`, and the two must
agree.

If you add a key, remove a key, or change a default, regenerate the table and commit the
result.

```bash
python3 -m ptm.engine_config --regenerate
git add src/python/ptm/engine_defaults.py
```

A change to an existing default needs particular care. No YAML file changes, no key is added
or removed, and every configuration that was accepted before is still accepted. A
configuration that omits the key then runs a different simulation, while the table in the
Python package still holds the previous value until it is regenerated.

The table is regenerated from the C++ sources, so no compiler, CSPICE or CMake step is needed,
and the model does not have to be built. Only the sources of your own checkout are read. Make
the C++ change first, then regenerate, because a table regenerated from sources without the
change will not contain it.

`src/python/tests/test_config_parity.py` derives the table from the C++ sources each time it
runs and fails if the result differs from the committed table. A table that was not
regenerated therefore fails in continuous integration. The test also fails when it can no
longer parse the structure of the parser, and when an example configuration in the repository
contains a key that the model does not read.

The test is run directly, without pytest.

```bash
python3 src/python/tests/test_config_parity.py
```

One of its five checks reads the configuration writer of the GUI and needs PyQt6. Without
PyQt6 it reports `SKIP`, and the summary counts skipped checks separately from failed ones.
Continuous integration installs PyQt6 and runs all five.

## The Python layer

`pip install -e .` from the repository root installs the `ptm` package and the `ptm` command.
The package contains no physics. It starts the model and reads its output.

Run the following tests before opening a pull request that changes `src/python/ptm/` or
`gui/`.

```bash
python3 src/python/tests/test_config_parity.py
python3 src/python/tests/test_results_api.py
python3 src/python/tests/test_engine_process.py
```

These are the three tests of the `python-layer` job in continuous integration. They need no
compiled model, SPICE kernels or GPU.

## Continuous integration

A failed Linux job blocks a merge. The macOS jobs temporarily do not (since 2026-11-04),
because the capacity of GitHub's macOS arm64 runners is limited. They run when a runner is
available, but a failed or missing macOS result does not fail the workflow. If your change
could behave differently on macOS (paths, compilers, Homebrew dependencies, PyQt6), open the
macOS job and check its result before merging. The workflow run carries a warning while this
arrangement is in effect, and the macOS jobs will block a merge again when capacity
recovers.
