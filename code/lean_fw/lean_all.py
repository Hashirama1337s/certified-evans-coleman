"""Lean extract of the Far Wall Collatz rows for EVERY bulk K-cell of both fw_cert3 bands (NEW file;
fw_cert3 / fw_cert2 imported unchanged; same setup() and bulk_cells_K() as the sealed runs).

For each K-cell: A_up = decimal upper bounds (strictly above the arb uppers, asserted), y = rounded float Perron
vector (more digits if needed), and Lean theorems  A_up y < y  (3 rows)  and  y > 0.  Output one .lean file.
usage (from <LAB>, under memcap):  py -3.11 -u lean_fw/lean_all.py lean_fw/fw_collatz_all.lean
"""
import sys, decimal, time, json
from fractions import Fraction
OUT = sys.argv[1]
sys.argv = ['fw_cert3.py']
sys.path.insert(0, '<LAB>/engine/cert')
import fw_cert3 as F3
import fw_cert2 as FW2
import numpy as np
arb = F3.arb
T0 = time.time()
decimal.getcontext().prec = 60
BANDS = [(2.5, 40.0, 55.0), (3.5, 55.0, 2000.0)]   # (M, Ka, Kmax) exactly as runs6/fw_cert3_M2.5_Ka40_K55.log, _M3.5_Ka55_K2000.log


def dec_up(a, digits=6):
    u = a.upper()
    d = decimal.Decimal(u.mid().str(40, radius=False).replace('e', 'E'))
    q = decimal.Decimal(1).scaleb(d.adjusted() - digits + 1)
    return str(d.quantize(q, rounding=decimal.ROUND_CEILING) + q)


def q(s):
    f = Fraction(s)
    return '(%d/%d : ℚ)' % (f.numerator, f.denominator)


lines = ['import Mathlib', '',
         '/-! Far Wall (fw_cert3): Collatz rows A_up * y < y and y > 0 for EVERY bulk K-cell of the two bands',
         '    [40, 55] (M = 2.5) and [55, 2000] (M = 3.5). Entries are decimal upper bounds strictly above the',
         '    interval-arithmetic uppers (asserted in lean_fw/lean_all.py). Only these rational inequalities are',
         '    kernel-checked here; the interval entries themselves come from ball arithmetic (arb). -/', '']
summary = []
n = 0
for (M, Ka, Kmax) in BANDS:
    cells, son, cdat, Kmax2 = F3.setup(M, Ka, Kmax)
    bulks = FW2.bulk_cells_K(cells, son, M, Ka, Kmax2, 1.03)
    assert bulks is not None
    for (a, b), r in bulks:
        A = [[dec_up(r['A'][i][j]) for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(3):
                assert arb(A[i][j]) >= r['A'][i][j].upper()
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
            raise SystemExit('no rational Collatz vector for K-cell [%s, %s]' % (a, b))
        ratio = max(float(lhs[i] / Fraction(ys[i])) for i in range(3))
        tag = 'c%03d' % n
        for i in range(3):
            terms = ' + '.join('%s * %s' % (q(A[i][j]), q(ys[j])) for j in range(3))
            lines.append('theorem fw_%s_row%d : %s < %s := by norm_num' % (tag, i, terms, q(ys[i])))
        lines.append('theorem fw_%s_ypos : (0 : ℚ) < %s ∧ (0 : ℚ) < %s ∧ (0 : ℚ) < %s := by norm_num' % (tag, q(ys[0]), q(ys[1]), q(ys[2])))
        summary.append(dict(cell=tag, M=M, K=[a, b], collatz_ia=r['collatz'].str(6), rational_ratio=ratio))
        n += 1
    print('[%6.0fs] band M=%s [%s, %s]: %d K-cells' % (time.time() - T0, M, Ka, Kmax, len(bulks)), flush=True)
open(OUT, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
json.dump(summary, open(OUT.replace('.lean', '.json'), 'w'), indent=0)
print('wrote %s: %d K-cells, worst rational ratio %.6f' % (OUT, n, max(s['rational_ratio'] for s in summary)))
