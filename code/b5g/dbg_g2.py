import sys
sys.path.insert(0, '<LAB>/b5g')
from flint import arb, ctx, arb_mat
import pipeline_g2 as PG
from tm import TM
ctx.prec = 200
sig0 = arb(arb('[8.9013232777556435191172160567756856 +/- 5.31e-35]').mid())
T = PG.ball(0, arb('0.05'))
Psi = PG.psi_centre(sig0)
for n in range(3):
    D = [TM([arb(0)] * (n + 2), T), TM([arb(0)] * (n + 2), T)]
    out = Psi(D)
    print(n, [o.c[n] for o in out])
    D = [TM([arb(0)] * n + [arb(1), arb(0)], T), TM([arb(0)] * (n + 2), T)]
    out = Psi(D)
    print(n, 'e0', [o.c[n] for o in out])
import pipeline_g2 as PG2
orig_solve = arb_mat.solve
for N in (26, 30, 36):
    try:
        yC, uC = PG.nl_solve(PG.psi_centre(sig0), [sig0, arb(-1) / 2], T, N, 2, 'affine')
        print('N', N, 'ok', yC[1].ev(arb('0.05')))
    except Exception as ex:
        print('N', N, 'fail', ex)
d = 2
N = 26
u = [[None] * (N + 1) for _ in range(d)]
def Dlow(n, extra=None):
    Ds = []
    for i in range(d):
        cs = [u[i][j] / (j + 1) for j in range(n)]
        cs.append(arb(0) if extra is None else extra[i])
        cs.append(arb(0))
        Ds.append(TM(cs, T))
    return Ds
for n in range(N + 1):
    base = Psi(Dlow(n))
    c = [base[i].c[n] for i in range(d)]
    L = arb_mat(d, d)
    for jj in range(d):
        e = [arb(0)] * d
        e[jj] = arb(1)
        col = Psi(Dlow(n, e))
        for i in range(d):
            L[i, jj] = col[i].c[n] - c[i]
    print(n, [b.n for b in base], ['%.2e %.2e' % (float(abs(x).mid()), float(x.rad())) for x in c])
    Mm = arb_mat(d, d)
    for i in range(d):
        for jj in range(d):
            Mm[i, jj] = (1 if i == jj else 0) - L[i, jj] / (n + 1)
    sol = Mm.solve(arb_mat([[ci] for ci in c]))
    for i in range(d):
        u[i][n] = sol[i, 0]
