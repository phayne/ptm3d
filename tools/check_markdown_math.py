#!/usr/bin/env python3
"""Flag math in Markdown that GitHub's renderer garbles.

    python3 tools/check_markdown_math.py                 # docs/, README.md, examples/, notebooks
    python3 tools/check_markdown_math.py FILE [FILE...]

GitHub parses Markdown before it hands math to MathJax, so some characters inside $...$ and
$$...$$ are consumed by the Markdown parser first. The rules below were measured against the
HTML GitHub serves for Markdown files (2026-11):

  star        '*' inside $...$ or $$...$$ pairs with another '*' as emphasis, and the math
              falls apart. Write \\ast (D^{\\ast} renders the same as D^*).
  escape      A backslash before ASCII punctuation (\\{ \\} \\, \\; \\! \\| \\_ \\% \\# ...) is a Markdown
              escape: GitHub drops the backslash. A '\\\\' line break survives only at the end
              of a line of a $$ block. Use \\lbrace \\rbrace \\thinspace \\negthinspace \\Vert
              (not \\thickspace: see macro), or write the expression as $`...`$ or in a ```math block.
  underscore  An '_' after a closing brace or parenthesis ('}_', ')_') can open emphasis, and a
              later '_' in the same paragraph can close it; everything between them is then
              italic text instead of math. Write that expression as $`...`$.
  dash        A '$' right after a dash does not open math: in $a$–$b$ the second expression is
              left as text. Write $a\\text{–}b$ or put spaces around the dash.
  table-pipe  '|' inside math in a table row ends the cell. Use \\lvert \\rvert or \\vert.
  display     A $$ block must stand alone, with blank lines before and after it and no
              indentation. Inside a list item or right after a line of text it is not math.
  macro       A macro GitHub's math renderer does not handle (tools/github_math_macros.txt,
              measured from GitHub's MathML): \\thickspace is printed as red text, and \\tag
              collapses the equation into a column. Write the operator alone (a = b), \\thinspace
              between a number and a unit, and \\qquad \\text{(B1)} for an equation number.
              Macros not yet seen working on GitHub are listed as warnings ("unverified"); they
              do not change the exit status.
  notebook    Jupyter does not understand $`...`$ or ```math; notebooks use $...$ and $$...$$.

Exit status 1 if anything other than an unverified macro is flagged.

Earlier versions passed pages that were broken on GitHub because they only modeled the
Markdown stage (is the math recognized?), not the renderer's macro set.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from github_math_macros import classify  # noqa: E402

FENCE = re.compile(r'^(\s*)(```+|~~~+)\s*(\S*)')
MATH_FENCE = re.compile(r'^[ \t]*```+math[ \t]*\n(.*?)^[ \t]*```', re.M | re.S)
ESCAPABLE = r'!"#%&\'()*+,\-./:;<=>?@\[\]^_`{|}~'
ESCAPE = re.compile(r'\\([' + ESCAPABLE + r'])')


def paragraphs(text):
    """Yield (first_line_number, lines) for every block outside code fences."""
    block, start, fence = [], None, None
    for n, line in enumerate(text.splitlines(), 1):
        m = FENCE.match(line)
        if fence is not None:
            if m and m.group(2)[0] == fence[0] and len(m.group(2)) >= len(fence):
                fence = None
            continue
        if m and m.group(3) != 'math':
            if block:
                yield start, block
            block, start, fence = [], None, m.group(2)
            continue
        if m and m.group(3) == 'math':
            if block:
                yield start, block
            block, start = [], None
            if m.group(1):
                yield n, ['\0indented-math-fence']
            fence = m.group(2)
            continue
        if not line.strip():
            if block:
                yield start, block
            block, start = [], None
            continue
        if start is None:
            start = n
        block.append(line)
    if block:
        yield start, block


def strip_code_spans(s):
    """Replace inline code (but not $`...`$ math) with spaces, keeping offsets."""
    out = list(s)
    i = 0
    while i < len(s):
        if s[i] == '\\':
            i += 2
            continue
        if s[i] == '`':
            n = len(s[i:]) - len(s[i:].lstrip('`'))
            close = s.find('`' * n, i + n)
            is_math = i > 0 and s[i - 1] == '$' and close >= 0 and s[close + n:close + n + 1] == '$'
            if close < 0:
                i += n
                continue
            if not is_math:
                for k in range(i, close + n):
                    out[k] = ' '
            i = close + n
            continue
        i += 1
    return ''.join(out)


def math_spans(par):
    """(kind, start, end) of math in one paragraph string; end is exclusive of the delimiter."""
    spans = []
    i = 0
    s = par
    while i < len(s):
        c = s[i]
        if c == '\\':
            i += 2
            continue
        if s.startswith('$`', i):
            close = s.find('`$', i + 2)
            if close >= 0:
                spans.append(('protected', i + 2, close))
                i = close + 2
                continue
        if s.startswith('$$', i):
            close = s.find('$$', i + 2)
            if close >= 0:
                spans.append(('display', i + 2, close))
                i = close + 2
                continue
            i += 2
            continue
        if c == '$' and i + 1 < len(s) and not s[i + 1].isspace():
            j = i + 1
            while j < len(s) and s[j] != '\n':
                if s[j] == '\\':
                    j += 2
                    continue
                if s[j] == '$' and not s[j - 1].isspace() and not s[j + 1:j + 2].isdigit():
                    break
                j += 1
            if j < len(s) and s[j] == '$':
                spans.append(('inline', i + 1, j))
                i = j + 1
                continue
        i += 1
    return spans


def _flank(s, k):
    """(left_flanking, right_flanking) for the single delimiter character at s[k]."""
    prev = s[k - 1] if k > 0 else ' '
    nxt = s[k + 1] if k + 1 < len(s) else ' '
    punct = lambda ch: not ch.isalnum() and not ch.isspace()
    left = not nxt.isspace() and (not punct(nxt) or prev.isspace() or punct(prev))
    right = not prev.isspace() and (not punct(prev) or nxt.isspace() or punct(nxt))
    return left, right


def macro_issues(line_no, tex):
    blocked, unverified = classify(tex)
    issues = [(line_no, 'macro', f'{m}: {why}') for m, why in blocked]
    issues += [(line_no, 'unverified', m) for m in unverified]
    return issues


def check_paragraph(start, lines, in_table):
    issues = []
    if lines == ['\0indented-math-fence']:
        return [(start, 'display', 'indented ```math block (inside a list) is shown as code')]
    par = '\n'.join(lines)
    par = strip_code_spans(par)
    spans = math_spans(par)
    line_of = lambda off: start + par.count('\n', 0, off)
    for kind, a, b in spans:
        issues += macro_issues(line_of(a), par[a:b])
        if kind == 'inline' and a >= 2 and par[a - 2] in '-\u2013\u2014':
            issues.append((line_of(a), 'dash', par[a:b]))
        if kind == 'protected':
            if in_table and '|' in par[a:b] and par.lstrip().startswith('|'):
                issues.append((line_of(a), 'table-pipe', par[a:b]))
            continue
        tex = par[a:b]
        if '*' in tex:
            issues.append((line_of(a), 'star', tex))
        for m in ESCAPE.finditer(tex):
            if m.group(1) == '\\':
                continue
            issues.append((line_of(a + m.start()), 'escape', tex))
            break
        if '\\\\' in tex:
            ok = kind == 'display' and all(
                seg.rstrip().endswith('\\\\') or '\\\\' not in seg for seg in tex.split('\n'))
            if not ok:
                issues.append((line_of(a), 'escape', tex))
        if in_table and '|' in tex and par.lstrip().startswith('|'):
            issues.append((line_of(a), 'table-pipe', tex))
        if kind == 'display':
            whole = par.strip()
            if not (whole.startswith('$$') and whole.endswith('$$') and whole.count('$$') == 2):
                if '\n' in tex or len(whole) > b - a + 4:
                    issues.append((line_of(a), 'display', tex))
    # underscore emphasis: an opener inside plain $ math followed by any closer in the paragraph
    plain = list(par)
    for kind, a, b in spans:
        if kind == 'protected':
            for k in range(a - 2, b + 2):
                plain[k] = ' '
    plain = ''.join(plain)
    unders = [k for k, ch in enumerate(plain) if ch == '_' and not (k > 0 and plain[k - 1] == '\\')]
    # Only inline spans: inside a $$ block GitHub writes underscore emphasis back as the same
    # characters, so it reaches MathJax unchanged.
    in_math = lambda k: any(kind == 'inline' and a <= k < b for kind, a, b in spans)
    for k in unders:
        if not in_math(k):
            continue
        left, right = _flank(plain, k)
        prev = plain[k - 1] if k > 0 else ' '
        can_open = left and (not right or (not prev.isalnum() and not prev.isspace()))
        if not can_open:
            continue
        for k2 in unders:
            if k2 <= k + 1:
                continue
            l2, r2 = _flank(plain, k2)
            nxt = plain[k2 + 1] if k2 + 1 < len(plain) else ' '
            if r2 and (not l2 or (not nxt.isalnum() and not nxt.isspace())):
                span = next(par[a:b] for kind, a, b in spans if a <= k < b)
                issues.append((line_of(k), 'underscore', span))
                break
        else:
            continue
        break
    return issues


def check_markdown(path):
    text = Path(path).read_text()
    issues = []
    lines = text.splitlines()
    for m in MATH_FENCE.finditer(text):
        issues += macro_issues(text.count('\n', 0, m.start(1)) + 1, m.group(1))
    for start, block in paragraphs(text):
        in_table = any(l.lstrip().startswith('|') for l in block)
        issues += check_paragraph(start, block, in_table)
        # display block placement: $$ opening a line inside a block that is not only that block
        for off, line in enumerate(block):
            # A one-line $$...$$ inside a list item is still math; a multi-line one is not.
            one_line = line.strip().startswith('$$') and line.strip().endswith('$$') and line.strip().count('$$') == 2
            if line.lstrip().startswith('$$') and line != line.lstrip() and not one_line \
                    and block[0] != '\0indented-math-fence':
                issues.append((start + off, 'display', line.strip()))
        if len(block) > 1 and any(l.strip() == '$$' for l in block):
            first, last = block[0].strip(), block[-1].strip()
            if not (first.startswith('$$') and last.endswith('$$')):
                issues.append((start, 'display', block[0].strip()))
    return issues


def check_notebook(path):
    nb = json.loads(Path(path).read_text())
    issues = []
    for i, cell in enumerate(nb.get('cells', [])):
        if cell.get('cell_type') != 'markdown':
            continue
        src = ''.join(cell.get('source', []))
        if re.search(r'\$`', src) or re.search(r'^\s*```math', src, re.M):
            issues.append((f'cell {i}', 'notebook', 'GitHub-only math syntax in a notebook'))
        for start, block in paragraphs(src):
            for line_no, rule, tex in check_paragraph(start, block, False):
                if rule in ('star', 'display', 'macro', 'unverified'):
                    issues.append((f'cell {i} line {line_no}', rule, tex))
    return issues


def default_files(root):
    files = [root / 'README.md']
    for sub in ('docs', 'examples'):
        files += sorted((root / sub).rglob('*.md'))
        files += sorted(p for p in (root / sub).rglob('*.ipynb') if '.ipynb_checkpoints' not in p.parts)
    return [f for f in files if f.exists()]


def main(argv):
    root = Path(__file__).resolve().parent.parent
    files = [Path(a) for a in argv] if argv else default_files(root)
    total = warned = 0
    for f in files:
        issues = check_notebook(f) if f.suffix == '.ipynb' else check_markdown(f)
        for where, rule, tex in issues:
            shown = re.sub(r'\s+', ' ', tex).strip()
            if rule == 'unverified':
                print(f'{f}:{where}: warning: [unverified] {shown} (not yet seen working on GitHub)')
                warned += 1
                continue
            print(f'{f}:{where}: [{rule}] {shown[:160]}')
            total += 1
    print(f'{len(files)} files checked, {total} problem(s), {warned} unverified macro(s)')
    return 1 if total else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
