#!/usr/bin/env python3
"""Which TeX macros GitHub's math renderer is known to handle (tools/github_math_macros.txt).

    from github_math_macros import classify
    blocked, unverified = classify(r'a \\thickspace = b')   # ([('\\thickspace', reason)], [])

    python3 tools/github_math_macros.py --derive SAVED.html   # list macros seen in a saved page

The list was measured from the MathML GitHub serves for a Markdown file view, not from
stock MathJax (which accepts macros GitHub's renderer does not). See the data file header.
"""
import html
import re
import sys
from pathlib import Path

DATA = Path(__file__).resolve().with_name('github_math_macros.txt')
# A control word, \begin{env} / \end{env} (classified by environment), or a control symbol.
MACRO = re.compile(r'\\(?:(?:begin|end)\{[A-Za-z*]+\}|[A-Za-z]+|.)', re.S)


def _key(tok):
    return '\\<space>' if tok in ('\\ ', '\\\n', '\\\t') else tok


def load(path=DATA):
    allow, block = set(), {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        kind, rest = line.split(None, 1)
        parts = rest.split(None, 1)
        if kind == 'allow':
            allow.add(parts[0])
        elif kind == 'block':
            block[parts[0]] = parts[1] if len(parts) > 1 else ''
        else:
            raise ValueError(f'{path}: unknown entry {line!r}')
    return allow, block


_ALLOW, _BLOCK = load()


def macros_in(tex):
    out = []
    for tok in MACRO.findall(tex):
        if tok.startswith('\\end{'):
            continue
        out.append(_key(tok))
    return out


def classify(tex):
    """(blocked [(macro, reason)], unverified [macro]) for one math expression, in order."""
    blocked, unverified = [], []
    for m in macros_in(tex):
        if m in _BLOCK:
            if m not in (b for b, _ in blocked):
                blocked.append((m, _BLOCK[m]))
        elif m not in _ALLOW and m not in unverified:
            unverified.append(m)
    return blocked, unverified


def derive(saved_html):
    """Macros in a saved GitHub file view: (rendered, red)."""
    h = Path(saved_html).read_text(encoding='utf-8', errors='replace')
    rendered, red = set(), set()
    for m in re.findall(r'<math-renderer[^>]*>(.*?)</math-renderer>', h, re.S):
        if '<mlabeledtr' in m:
            # No red text, but browsers do not lay out labeled rows (see the data file).
            red.add('\\tag')
        for attrs in re.findall(r'<\w+([^>]*data-latex="[^"]*"[^>]*)>', m):
            lat = html.unescape(re.search(r'data-latex="([^"]*)"', attrs).group(1))
            first = MACRO.match(lat)
            if 'mathcolor="red"' in attrs:
                red.add(_key(lat.strip()))
            elif first:
                rendered.add(_key(first.group(0)))
        if 'mathcolor="red"' not in m and '<merror' not in m:
            top = re.search(r'<math[^>]*data-latex="([^"]*)"', m)
            if top:
                rendered.update(macros_in(html.unescape(top.group(1))))
    return rendered - red - {'\\tag'}, red


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--derive':
        ok, red = derive(sys.argv[2])
        for m in sorted(ok):
            print(f'allow {m}' + ('   # new' if m not in _ALLOW else ''))
        for m in sorted(red):
            print(f'red   {m}')
        sys.exit(0)
    print(__doc__)
    sys.exit(2)
