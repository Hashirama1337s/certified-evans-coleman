"""B4 speedup seat (NEW file; fw_cert.py is imported, never edited): Far Wall with a lower start K_a.

fw_cert.certify (PROOF-B4-v4 sec. 6) is reused verbatim for every mathematical object: bulk_ia, centre_theta,
centre_test_ia, the background cells and the centre data.  This file only changes DISCRETISATION knobs and
adds diagnostics:
  * ns   : number of log-cells of the centre Volterra majorant (fw_cert default 150);
  * NB1  : cells of the signed Born quadrature B1 (default 1500);
  * nth  : theta-cells on [0, pi/2] (default 160) and an adaptive bisection of failing theta-cells (each
           half is an ordinary fw_cert theta-cell, so the union still covers [0, pi/2]);
  * kratio: K-cell ratio (default 1.03), with adaptive bisection of failing K-cells.
Every accepted (theta-cell, K-cell) pair passes exactly fw_cert.centre_test_ia; the certificate statement is
the fw_cert one on the union of the cells.  Diagnostics: per theta-cell the worst ratio lower/rhs.

usage:
  fw_cert2.py probe M Ka Kmax nth step ns nb1        (sampled theta-cells, ratios; NON-COUNTING)
  fw_cert2.py run   M Ka Kmax nth ns nb1 [depth]      (full certificate with adaptive bisection depth)
"""
import sys, os, math, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import fw_cert as FW
from flint import arb, acb
from core import sym

T0 = time.time()
def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)


def setup(M, Ka, Kmax):
    T1 = time.time()
    cells = FW.eval_cells(FW.tube_cells())
    tube_n = len(cells)
    xbar_K = float((arb(M)/cells[-1]['xlo'].exp()).lower())
    if Kmax > xbar_K:
        th_min = float((FW.PS*arb(M)/arb(Kmax)).lower())*0.97
        cells = cells + FW.centre_bulk_cells(th_min)
        log('[setup] bulk extended below xbar: %d centre-manifold cells (theta >= %.3g)' % (len(cells) - tube_n, th_min))
    son = FW.sonic_piece(); cells, son = FW.prep_bulk(cells, son)
    xmin = cells[-1]['xlo']
    Kcap = float((arb(M)/xmin.exp()).lower())
    if Kmax > Kcap:
        log('Kmax %.1f lowered to %.2f (deepest bulk cell)' % (Kmax, Kcap)); Kmax = Kcap
    cm = FW.centre_manifold_data(); FW.cdat_C01, FW.cdat_C12 = cm['C01'], cm['C12']
    vneed = FW.PS*arb(M)/arb(Ka)
    cdat = cm['data'] + [FW.bt_centre_data(cm['C01'], cm['C12'])] + FW.tube_centre_data(cells, cm['C01'], cm['C12'], vneed)
    log('[setup] cells %d, sonic g=%s, A_max=%s, centre data %d cells (%.0fs)' % (len(cells), son['g'].str(4), son['Amax'].str(4), len(cdat), time.time() - T1))
    return cells, son, cdat, Kmax


def bulk_cells_K(cells, son, M, Ka, Kmax, kratio):
    Ks = [Ka]
    while Ks[-1] < Kmax: Ks.append(min(Kmax, Ks[-1]*kratio))
    out = []
    for a, b in zip(Ks[:-1], Ks[1:]):
        r = FW.bulk_ia(cells, son, M, arb(a), arb(b))
        if not r['ok']:
            log('BULK FAILS on K in [%.3f, %.3f]: collatz %s' % (a, b, r['collatz'].str(4))); return None
        out.append(((a, b), r))
    log('[bulk] %d K-cells in [%.2f, %.2f] contractive; worst collatz %s, max e_b %s, max nbar %s' % (
        len(out), Ka, Kmax, max(r['collatz'] for _, r in out).str(4), max(r['eb'] for _, r in out).str(4), max(r['nbar'] for _, r in out).str(4)))
    return out


def theta_cell(cells, son, M, Ka, cdat, bulks, th_lo, th_hi, ns, depth, maxdepth, stats):
    """certify one theta-cell against all K-cells; bisect theta on failure (depth <= maxdepth).
    returns (ok, worst ratio lower/rhs over the K-cells, worst margin)"""
    cen = FW.centre_theta(th_lo, th_hi, Ka, M, cdat, ns=ns)
    stats['centre'] += 1
    worst = None; wm = None; failK = None
    if cen['w'] is not None and cen['N'] is not None:
        for (a, b), bk in bulks:
            ok, info = FW.centre_test_ia(bk, cen, a, b, th_lo, th_hi)
            if isinstance(info, str):
                ratio = 0.0
            else:
                ratio = float((info['lower']/info['rhs']).lower()) if info['rhs'] > 0 else float('inf')
                m = float((info['lower'] - info['rhs']).lower())
                wm = m if wm is None else min(wm, m)
            if worst is None or ratio < worst[0]: worst = (ratio, a, b)
            if not ok:
                failK = (a, b); break
    else:
        failK = 'super'
    if failK is None:
        stats['ns_used'][ns] = stats['ns_used'].get(ns, 0) + 1
        return True, worst, wm
    # ns ladder (env FW2_NS_LADDER, e.g. "150,300,600"): retry the SAME cell with a finer centre majorant
    ladder = [int(x) for x in os.environ.get('FW2_NS_LADDER', '').split(',') if x]
    finer = [x for x in ladder if x > ns]
    if finer:
        return theta_cell(cells, son, M, Ka, cdat, bulks, th_lo, th_hi, finer[0], depth, maxdepth, stats)
    if depth >= maxdepth:
        return False, worst, wm
    mid = (th_lo + th_hi)/2
    stats['split'] += 1
    ns0 = ladder[0] if ladder else ns
    o1 = theta_cell(cells, son, M, Ka, cdat, bulks, th_lo, mid, ns0, depth + 1, maxdepth, stats)
    o2 = theta_cell(cells, son, M, Ka, cdat, bulks, mid, th_hi, ns0, depth + 1, maxdepth, stats)
    w = min([x for x in (o1[1], o2[1]) if x is not None], default=None, key=lambda t: t[0])
    ms = [x for x in (o1[2], o2[2]) if x is not None]
    return o1[0] and o2[0], w, (min(ms) if ms else None)


def run(M, Ka, Kmax, nth, ns, nb1, maxdepth, sub=None, kratio=float(os.environ.get("FW2_KRATIO", "1.03"))):
    FW.NB1 = nb1
    log('fw_cert2: M=%s Ka=%s Kmax=%s theta cells=%d ns=%d NB1=%d bisection depth %d kratio %s prec %d (fw_cert objects unchanged)' % (
        M, Ka, Kmax, nth, ns, nb1, maxdepth, kratio, FW.PREC))
    cells, son, cdat, Kmax = setup(M, Ka, Kmax)
    bulks = bulk_cells_K(cells, son, M, Ka, Kmax, kratio)
    if bulks is None:
        log('VERDICT: bulk not certified'); return None
    pi2 = arb.pi()/2
    ths = [pi2*k/nth for k in range(nth + 1)]
    idx = range(nth) if sub is None else sub
    stats = dict(centre=0, split=0, ns_used={})
    fails = []; wmin = None; rmin = None
    for k in idx:
        t1 = time.time()
        ok, worst, wm = theta_cell(cells, son, M, Ka, cdat, bulks, ths[k], ths[k + 1], ns, 0, maxdepth, stats)
        if wm is not None: wmin = wm if wmin is None else min(wmin, wm)
        if worst is not None: rmin = worst if rmin is None else min(rmin, worst, key=lambda t: t[0])
        if not ok: fails.append(k)
        if sub is not None or not ok:
            log('  theta cell %d [%.4f, %.4f]: %s  worst lower/rhs %.4f at K in [%.2f, %.2f]  (%.0fs)' % (
                k, float(ths[k].mid()), float(ths[k + 1].mid()), 'ok' if ok else 'FAIL', worst[0] if worst else float('nan'),
                worst[1] if worst else 0, worst[2] if worst else 0, time.time() - t1))
    log('ns ladder %s, accepted cells by ns %s' % (os.environ.get('FW2_NS_LADDER', '-'), stats['ns_used']))
    log('theta cells %d (+%d bisections, %d centre evaluations), K-cells %d, failures %s, worst lower/rhs %s, worst margin %s' % (
        len(idx), stats['split'], stats['centre'], len(bulks), fails, rmin, wmin))
    if sub is None and not fails:
        log('VERDICT (CERTIFIED, interval arithmetic on the E5 enclosure, fw_cert objects): no physical eigenvalue kappa with Re kappa >= 0 and %.4g <= |kappa| <= %.4g' % (Ka, Kmax))
    elif sub is None:
        log('VERDICT: NOT certified at Ka=%s (failing theta cells %s)' % (Ka, fails))
    return fails


if __name__ == '__main__':
    mode = sys.argv[1]
    M = float(sys.argv[2]); Ka = float(sys.argv[3]); Kmax = float(sys.argv[4]); nth = int(sys.argv[5])
    if mode == 'probe':
        st = int(sys.argv[6]); ns = int(sys.argv[7]); nb1 = int(sys.argv[8])
        md = int(sys.argv[9]) if len(sys.argv) > 9 else 0
        sub = sorted(set(list(range(0, nth, st)) + [nth - 1]))
        if os.environ.get('FW2_SUB'):
            a_, b_ = os.environ['FW2_SUB'].split('-'); sub = list(range(int(a_), int(b_) + 1, st))
        run(M, Ka, Kmax, nth, ns, nb1, md, sub=sub)
    elif mode == 'run':
        ns = int(sys.argv[6]); nb1 = int(sys.argv[7]); md = int(sys.argv[8]) if len(sys.argv) > 8 else 0
        run(M, Ka, Kmax, nth, ns, nb1, md)
