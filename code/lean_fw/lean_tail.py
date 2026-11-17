"""Lean extract of the Far Wall Collatz rows for the UNIFORM TAIL matrix (all K >= K_T = 2000, M = 3.5).
NEW file; engine/cert/fw_tail2.py is imported UNCHANGED and its bulk_tail(K_T) is called exactly as its
`run` / `bulk` modes do (runs6/fw_tail2_KT2000.log: weights [3.125, 0.6, 250.0, 19/20, 11/20], Collatz 0.54071).

A_up = decimal upper bounds strictly above the arb uppers (asserted, as in lean_all.py), y = rounded float
Perron vector of A_up, Lean theorems  A_up y < y  (3 rows) and  y > 0.  Output one .lean file (+ .json summary).
usage (from <LAB>, under memcap 8):  py -3.11 -u lean_fw/lean_tail.py lean_fw/fw_collatz_tail.lean
"""
import sys, decimal, time, json, hashlib
from fractions import Fraction
OUT = sys.argv[1]
KT = 2000.0
sys.argv = ['fw_tail2.py']
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import fw_tail2 as FT
import numpy as np
arb = FT.arb
T0 = time.time()
decimal.getcontext().prec = 60


def dec_up(a, digits=6):
    u = a.upper()
    d = decimal.Decimal(u.mid().str(40, radius=False).replace('e', 'E'))
    q = decimal.Decimal(1).scaleb(d.adjusted() - digits + 1)
    return str(d.quantize(q, rounding=decimal.ROUND_CEILING) + q)


def q(s):
    f = Fraction(s)
    return '(%d/%d : ℚ)' % (f.numerator, f.denominator)


FT.provenance()
b = FT.bulk_tail(KT, verbose=True)
assert b['ok']
print('[%6.0fs] bulk_tail weights %s, Collatz %s' % (time.time() - T0, b.get('w'), b['collatz'].str(6)), flush=True)
A = [[dec_up(b['A'][i][j]) for j in range(3)] for i in range(3)]
for i in range(3):
    for j in range(3):
        assert arb(A[i][j]) > b['A'][i][j].upper()      # strictly above the arb upper
        assert Fraction(A[i][j]) >= 0
Af = np.array([[float(x) for x in row] for row in A])
w, V = np.linalg.eig(Af)
yv = np.abs(V[:, np.argmax(abs(w))].real); yv = yv / yv.max()
for dig in (4, 6, 8, 10):
    ys = [('%.' + str(dig) + 'f') % max(t, 10.0 ** -dig) for t in yv]
    lhs = [sum(Fraction(A[i][j]) * Fraction(ys[j]) for j in range(3)) for i in range(3)]
    if all(lhs[i] < Fraction(ys[i]) for i in range(3)) and all(Fraction(s) > 0 for s in ys):
        break
else:
    raise SystemExit('no rational Collatz vector for the tail matrix')
ratio = max(float(lhs[i] / Fraction(ys[i])) for i in range(3))
lines = ['import Mathlib', '',
         '/-! Far Wall tail (fw_tail2, K_T = 2000, M = 3.5): Collatz rows A_up * y < y and y > 0 for the UNIFORM',
         '    bulk matrix valid for every K >= 2000 (PROOF-B4-v4-AMEND2-TAIL sec. 2). Entries are decimal upper bounds',
         '    strictly above the ball-arithmetic (arb) uppers of fw_tail2.bulk_tail(2000) (asserted in lean_fw/lean_tail.py).',
         '    Only these rational inequalities are kernel-checked here. -/', '']
for i in range(3):
    terms = ' + '.join('%s * %s' % (q(A[i][j]), q(ys[j])) for j in range(3))
    lines.append('theorem fw_tail_row%d : %s < %s := by norm_num' % (i, terms, q(ys[i])))
lines.append('theorem fw_tail_ypos : (0 : ℚ) < %s ∧ (0 : ℚ) < %s ∧ (0 : ℚ) < %s := by norm_num' % (q(ys[0]), q(ys[1]), q(ys[2])))
open(OUT, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
summary = dict(KT=KT, weights=[str(x) for x in b.get('w')], collatz_ia=b['collatz'].str(6), A_arb=[[a.str(8) for a in row] for row in b['A']],
               A_up=A, y=ys, rational_ratio=ratio,
               sha256_fw_tail2=hashlib.sha256(open('<LAB>/engine/cert/fw_tail2.py', 'rb').read()).hexdigest())
json.dump(summary, open(OUT.replace('.lean', '.json'), 'w'), indent=1)
print('A_up =', A)
print('y =', ys)
print('wrote %s: rational Collatz ratio %.6f (arb %s)' % (OUT, ratio, b['collatz'].str(6)))
