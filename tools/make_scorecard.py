#!/usr/bin/env python
"""make_scorecard.py -- build scorecard.svg and SCORECARD.md: every sealed test of this project, failures included.

The nine seals and their outcomes are those of the lab ledger (paper Section 9). The four HOLD tests are checked against
the ledger export results/ledger-lightball.json (their settle rows); the script stops if an outcome disagrees. The
FAILED and VOID attempts are listed from the paper's Section 9 table; their seal and result files stay in the lab and are
not shipped (they certify nothing).
Run from anywhere:  python tools/make_scorecard.py
"""
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, 'results')

# (section, [(test, seal, what, outcome, detail, source or None)])
TESTS = [
    ('EXISTENCE (INTERVAL ARITHMETIC)', [
        ('E', '9ed31dfab0', 'Existence, attempt 1', 'FAILED', 'k = 1/2 warm-up did not certify', None),
        ('E2', '95f6de5cf2', 'Existence, attempt 2', 'FAILED', 'warm-up ladder did not contract', None),
        ('E3', '8aec628fd2', 'Existence, attempt 3', 'VOID',
         'settled HOLD, then corrected to VOID: a library series cap truncated the sonic tail input', None),
        ('E4', '7928715efe', 'Existence, attempt 4', 'FAILED', 'k = 1/2 warm-up did not contract', None),
        ('E5', 'd0af75803f', 'Existence, attempt 5 (mean-value head bound R8)', 'HOLD',
         'every item HOLD; controls and three content nulls pass', 'e5/BAR-LB-E5.md'),
        ('E5 dup', '9592b015f0', 'The same E5 file sealed twice', 'VOID', 'duplicate, never run or used', None),
        ('B5C', '71904980c4', 'Blind second existence proof (closure checks)', 'HOLD',
         'blind seat, Chebyshev + Poincare-Miranda; agrees with E5 to 19 digits (POST HOC)', 'b5c/RESULT.md'),
    ]),
    ('THE COUNT AND THE EXPONENT', [
        ('M', 'b1b0ee1cde', 'Exactly one growing spherical mode; kappa*', 'HOLD',
         'count N = 1 on the half-disc of radius 40 + Far Wall to infinity', 'm/BAR-LB-M.md'),
        ('G2', '8a1c27c7fb', 'Blind second kappa* (on the B5 solution)', 'HOLD',
         'inside M\'s enclosure; instruments pinned after the run (disclosed)', 'g2/BAR-LB-G2.md'),
    ]),
]


def export_outcomes():
    rows = json.load(open(os.path.join(R, 'ledger-lightball.json'), encoding='utf-8'))
    return {r['seal_id']: r['outcome'] for r in rows if r.get('kind') == 'settle' and 'outcome' in r}


FONT = 'Segoe UI, Helvetica, Arial, sans-serif'
INK, MUTED, ACCENT, BAND, LINE, BORDER = '#e8edf2', '#9aa4ae', '#39ff14', '#262b31', '#22272e', '#2f353d'
CHIP = {'HOLD': ('#12301c', '#39ff14'), 'FAILED': ('#3a1a1a', '#ff6b6b'), 'VOID': ('#33290f', '#ffc94d')}


def tw(s, per=7.4):
    return round(len(s) * per)


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def build():
    exp = export_outcomes()
    counts = {'HOLD': 0, 'FAILED': 0, 'VOID': 0}
    for _, tests in TESTS:
        for tid, seal, _, outcome, _, _ in tests:
            counts[outcome] += 1
            if outcome == 'HOLD' and exp.get(seal) != 'HOLD':
                sys.exit('ledger export disagrees for %s %s: %r' % (tid, seal, exp.get(seal)))
    if set(exp) - {s for _, ts in TESTS for _, s, *_ in ts}:
        sys.exit('ledger export has a settled seal not in the table')
    n = sum(counts.values())
    title = 'Scorecard: %d sealed tests (%d HOLD, %d FAILED, %d VOID)' % (n, counts['HOLD'], counts['FAILED'], counts['VOID'])
    sub = 'Sealed = prediction and pass bar hashed before the run. Failures are shown as failures. Details: SCORECARD.md.'
    W, y, body = 900, 66, []
    for head, rows in TESTS:
        body.append(f'<rect x="10" y="{y}" width="{W - 20}" height="22" fill="{BAND}"/>')
        body.append(f'<text x="20" y="{y + 16}" font-size="11.5" font-weight="600" fill="{ACCENT}" letter-spacing="0.5">{esc(head)}</text>')
        y += 22
        for tid, seal, what, outcome, _, _ in rows:
            body.append(f'<text x="20" y="{y + 22}" font-size="12.5" font-family="Consolas, Menlo, monospace" fill="{INK}">{esc(tid + " " + seal)}</text>')
            body.append(f'<text x="190" y="{y + 22}" font-size="12.5" fill="{INK}">{esc(what)}</text>')
            fill, fg = CHIP[outcome]
            t = tw(outcome)
            body.append(f'<rect x="620" y="{y + 8}" width="{t + 20}" height="20" rx="10" fill="{fill}"/>')
            body.append(f'<text x="630" y="{y + 22}" font-size="11.5" font-weight="600" fill="{fg}" '
                        f'textLength="{t}" lengthAdjust="spacingAndGlyphs">{outcome}</text>')
            y += 30
            body.append(f'<line x1="10" y1="{y + 1}" x2="{W - 10}" y2="{y + 1}" stroke="{LINE}"/>')
        y += 8
    H = y + 6
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="{FONT}" '
           f'role="img" aria-label="{esc(title)}">\n<title>{esc(title)}</title>\n'
           f'<rect width="{W}" height="{H}" rx="10" fill="#1c2026" stroke="{BORDER}"/>\n'
           f'<text x="20" y="34" font-size="18" font-weight="600" fill="{INK}">{esc(title)}</text>\n'
           f'<text x="20" y="54" font-size="12" fill="{MUTED}">{esc(sub)}</text>\n' + '\n'.join(body) + '\n</svg>\n')
    with open(os.path.join(ROOT, 'scorecard.svg'), 'w', encoding='utf-8', newline='\n') as f:
        f.write(svg)
    chip = {k: '![%s](brand/chips/%s.svg)' % (k, k.lower()) for k in CHIP}
    md = ['# Scorecard', '', title + '.', '',
          'Generated by `tools/make_scorecard.py`. A sealed test is one whose predictions and pass bars were hashed before',
          'its run. The four HOLD outcomes are checked against the ledger export `results/ledger-lightball.json`; the',
          'FAILED and VOID attempts are listed as in paper Section 9 (their files certify nothing and are not shipped).', '',
          '| Section | Test | Seal | What was predicted | Verdict | As scored | Source |',
          '|---|---|---|---|---|---|---|']
    for head, rows in TESTS:
        sec = head.split(' (')[0].capitalize()
        for tid, seal, what, outcome, detail, src in rows:
            link = f'[results/{src}](results/{src})' if src else 'paper.md Section 9'
            md.append(f'| {sec} | {tid} | `{seal}` | {what} | {chip[outcome]} | {detail} | {link} |')
    md += ['', 'Notes:',
           '- Amendments sealed under the HOLD tests: bee365c702 (B5C), b92e091b0c and 4971bef2bf (M), b62c6aecc9 (G2).',
           '- E3 is the instructive failure: it first scored HOLD, and the lab corrected it to VOID when a library default',
           '  (a series cap of 10 terms) was found to have truncated an input. The printed tail 3.7e-59 was really about',
           '  3.7e-35. E5 re-did the step with the cap raised and a run-time guard.',
           '- The clean-room replay of all four HOLD tests is in results/CLEANROOM-REPORT.md.']
    with open(os.path.join(ROOT, 'SCORECARD.md'), 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(md) + '\n')
    print(title)


if __name__ == '__main__':
    build()
