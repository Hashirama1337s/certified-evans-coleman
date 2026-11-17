"""Draw the README figures of certified-evans-coleman from the data the repository ships.

Usage:  python figures/make_figures.py <repo root> [--png]
Needs Python 3.11 and matplotlib 3.11 (nothing else). Writes figures/*.svg next to this script
(or into ./figures when the script does not sit in a folder named figures). --png also writes PNG previews.

Every plotted number is read from shipped files:
  results/m/runs7/count_R40.jsonl       the 99 certified contour pieces of the sealed count (endpoints, arg increment)
  results/m/runs7/verdict_M.json        N = 1, Delta/pi, kappa*, gamma, kappa_g (sealed test M)
  results/m/runs6/fw_cert3_*.log,
  results/m/runs6/fw_tail2_KT2000.log   the three Far Wall verdict lines (band limits 40, 55, 2000, infinity)
  results/e5/runs6/target.log           E5 enclosures of V0 and W0 (PE1 lines)
  code/b5/results_summary.json          B5 enclosures of V_s and sigma0 (blind second proof)
  results/pickles/lin_fw4_bg3.pkl       certified background cover (E5 instrument): start boxes of every integration
                                        step from the sonic point inward; checked against MANIFEST.sha256 and loaded
                                        with an unpickler that refuses every class (plain lists, dicts and strings only)
"""
import hashlib
import math
import io
import json
import os
import pickle
import re
import sys
from decimal import Decimal, getcontext

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle
from matplotlib.lines import Line2D

getcontext().prec = 80

# ---------------------------------------------------------------- brand style (figures/STYLE.md)
BG, BORDER, GRID, AXC, TXT, ACC = "#1c2026", "#2f353d", "#22272e", "#9aa4ae", "#e8edf2", "#39ff14"
CAT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]
W_IN = 12.5  # 900 pt viewBox width

plt.rcParams.update({
    "svg.fonttype": "none", "svg.hashsalt": "moki-julio",
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 16, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15,
    "axes.facecolor": BG, "figure.facecolor": "none", "savefig.facecolor": "none",
    "axes.edgecolor": AXC, "axes.labelcolor": TXT, "xtick.color": AXC, "ytick.color": AXC,
    "text.color": TXT, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.fontsize": 14, "lines.linewidth": 2.0, "axes.axisbelow": True,
})


def new_fig(h_in):
    fig = plt.figure(figsize=(W_IN, h_in))
    fig.patches.append(FancyBboxPatch((0.004, 0.006), 0.992, 0.988, boxstyle="round,pad=0,rounding_size=0.015",
                                      transform=fig.transFigure, facecolor=BG, edgecolor=BORDER, linewidth=1.5,
                                      zorder=-10))
    return fig


def title(fig, main, sub):
    h = fig.get_figheight()
    fig.text(0.035, 1 - 0.42 / h, main, fontsize=21, fontweight="bold", color=TXT, va="top")
    fig.text(0.035, 1 - 0.88 / h, sub, fontsize=15, color=AXC, va="top")


PNG = "--png" in sys.argv


def save(fig, out, name):
    path = os.path.join(out, name)
    fig.savefig(path, format="svg", metadata={"Date": None, "Creator": None})
    with open(path, encoding="utf-8") as fh:  # mathtext pieces name one font only; give them the full stack
        svg = fh.read().replace("font-family: 'Segoe UI';", "font-family: 'Segoe UI', 'Helvetica', 'Arial', sans-serif;")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(svg)
    if PNG:
        prev = os.path.join(out, os.pardir, "_preview")
        os.makedirs(prev, exist_ok=True)
        fig.savefig(os.path.join(prev, name.replace(".svg", ".png")), format="png", dpi=144)
    plt.close(fig)
    print("wrote figures/" + name)


def rd(root, rel):
    with open(os.path.join(root, rel), encoding="utf-8-sig", errors="replace") as fh:
        return fh.read()


BALL = re.compile(r"\[\s*([-+0-9.eE]*)\s*\+/-\s*([0-9.eE+-]+)\s*\]")


def ball(s):
    """'[mid +/- rad]' (arb print) -> (Decimal mid, Decimal rad); '[+/- r]' has mid 0."""
    m = BALL.fullmatch(s.strip())
    if not m:
        return Decimal(s.strip()), Decimal(0)
    return Decimal(m.group(1) or "0"), Decimal(m.group(2))


SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci_up(x, nd=1):
    """like sci, but the mantissa is rounded up (radii, widths, upper bounds are never understated)."""
    x = Decimal(str(x))
    e = x.adjusted()
    m = (x.scaleb(-e) * 10 ** nd).to_integral_value(rounding="ROUND_CEILING") / Decimal(10 ** nd)
    if m >= 10:
        m, e = m / 10, e + 1
    return "%.*f × 10%s" % (nd, m, str(e).translate(SUP))


def rounded_common(ends, nmax=40):
    """largest n <= nmax such that every endpoint rounds (half-even) to the same n-decimal string; that string."""
    from decimal import ROUND_HALF_EVEN
    for n in range(nmax, 0, -1):
        q = Decimal(1).scaleb(-n)
        vals = {str(Decimal(x).quantize(q, rounding=ROUND_HALF_EVEN)) for x in ends}
        if len(vals) == 1:
            return vals.pop(), n
    return None, 0


def sci(x, nd=1):
    """2.42e-5 -> '2.4 × 10⁻⁵' (plain text for labels)."""
    m, e = ("%.*e" % (nd, float(x))).split("e")
    return "%s × 10%s" % (m, str(int(e)).translate(SUP))


# ---------------------------------------------------------------- figure 1: the count (hero)
def load_count(root):
    lines = rd(root, "results/m/runs7/count_R40.jsonl").splitlines()
    head = json.loads(lines[0])
    assert head.get("header") and head["R1"] == 40.0
    pieces = [json.loads(l) for l in lines[1:] if l.strip()]
    assert len(pieces) == head["ndisc_plan"] == 99
    verdict = json.loads(rd(root, "results/m/runs7/verdict_M.json"))
    return head, pieces, verdict


def far_wall(root):
    """[(lo, hi)] from the three certified VERDICT lines (hi = None means infinity)."""
    out = []
    for rel in ("results/m/runs6/fw_cert3_M2.5_Ka40_K55.log", "results/m/runs6/fw_cert3_M3.5_Ka55_K2000.log",
                "results/m/runs6/fw_tail2_KT2000.log"):
        v = [l for l in rd(root, rel).splitlines() if "VERDICT (CERTIFIED" in l]
        assert len(v) == 1, rel
        m = re.search(r"Re kappa >= 0 and (\d+) <= \|kappa\| <= (\d+)", v[0])
        if m:
            out.append((int(m.group(1)), int(m.group(2))))
        else:
            m = re.search(r"Re kappa >= 0 and \|kappa\| >= (\d+)", v[0])
            out.append((int(m.group(1)), None))
    return out


def fig_count(root, out):
    head, pieces, verdict = load_count(root)
    bands = far_wall(root)
    assert bands[0][0] == 40 and bands[-1][1] is None
    pts = [complex(pieces[0]["a"])] + [complex(p["b"]) for p in pieces]
    # cumulative arg increment (ball sums: mids add, radii add)
    mids, rads, cum_m, cum_r = [], [], Decimal(0), Decimal(0)
    for p in pieces:
        m, r = ball(p["arg"])
        cum_m += m
        cum_r += r
        mids.append(cum_m)
        rads.append(cum_r)
    pi = Decimal("3.14159265358979323846264338327950288419716939937510")
    total, total_r = cum_m / pi, cum_r / pi
    assert abs(total - 1) + total_r < Decimal("0.25") and verdict["PM1"]["N_R40"] == 1
    kstar = ball(verdict["PM2"]["kappa_star"])[0]
    kg = ball(verdict["controls"]["C1c_kappa_g"])[0]
    gamma = verdict["PM2"]["gamma"]

    fig = new_fig(8.0)
    title(fig, "Exactly one way to fall",
          "Counting the growing disturbances of the collapsing ball of light, with a guaranteed error bar")

    # left: the plane of growth rates (upper half; the lower half is its mirror image)
    ax = fig.add_axes([0.075, 0.115, 0.40, 0.70])
    R0 = bands[0][0]
    lim = 64
    th = [math.pi / 2 - math.pi / 2 * k / 90 for k in range(91)]
    ring = [(R0 * math.cos(t), R0 * math.sin(t)) for t in th]
    ring += [(lim, 0.0), (lim, lim), (0.0, lim)]  # closed inside the plotted window: no reliance on clip paths
    ax.add_patch(Polygon(ring, closed=True, facecolor=CAT[2], alpha=0.25, edgecolor=CAT[2], lw=1, zorder=0))
    xs, ys = [z.real for z in pts], [z.imag for z in pts]
    ax.plot(xs, ys, color=CAT[0], lw=2.6, zorder=3)
    for i0, i1 in ((10, 12), (30, 32)):
        ax.annotate("", xy=(pts[i1].real, pts[i1].imag), xytext=(pts[i0].real, pts[i0].imag),
                    arrowprops=dict(arrowstyle="-|>", color=CAT[0], lw=2, mutation_scale=24), zorder=4)
    ax.annotate("", xy=(pts[-20].real, pts[-20].imag), xytext=(pts[-24].real, pts[-24].imag),
                arrowprops=dict(arrowstyle="-|>", color=CAT[0], lw=2, mutation_scale=24), zorder=4)
    ax.plot([float(kstar)], [0], "o", ms=11, color=CAT[1], zorder=5)
    ax.set_xlim(-6, lim)
    ax.set_ylim(-6, lim)
    ax.set_aspect("equal")
    ax.axhline(0, color=AXC, lw=0.8, zorder=1)
    ax.axvline(0, color=AXC, lw=0.8, zorder=1)
    ax.set_xlabel("growth rate κ: real part (> 0 means growing)")
    ax.set_ylabel("imaginary part (ringing)")
    ax.text(52, 52, "Far Wall:\nproved empty\nfor |κ| ≥ %d,\nall the way out" % R0, color=TXT, fontsize=14,
            ha="center", va="center")
    ax.text(5, 30, "counting\npath,\nradius %d" % R0, color=CAT[0], fontsize=14, ha="left", va="center")
    ax.text(5, 9, "inside: exactly\none zero", color=TXT, fontsize=14, ha="left", va="center")
    ax.annotate("κ*", (float(kstar), 0), xytext=(9, -4), color=CAT[1], fontsize=14, va="center",
                arrowprops=dict(arrowstyle="-", color=CAT[1], lw=1.2))

    # right top: zoom on the origin (the one mode and the removed clock-shift mode)
    ix = fig.add_axes([0.58, 0.685, 0.39, 0.105])
    ix.plot([0, 0], [0, 1.2], color=CAT[0], lw=2.6)
    ix.axhline(0, color=AXC, lw=0.8)
    ix.plot([float(kstar)], [0], "o", ms=11, color=CAT[1], zorder=5)
    ix.plot([float(kg)], [0], "o", ms=10, mfc=BG, mec=AXC, mew=2, zorder=5)
    ix.set_xlim(-0.2, 4.6)
    ix.set_ylim(-1.75, 1.1)
    ix.set_yticks([])
    ix.spines["left"].set_visible(False)
    ix.grid(False)
    ix.tick_params(labelsize=13)
    ix.text(float(kstar) + 0.12, 0.08, "κ* = %.5f\nthe one growing mode" % float(kstar), color=CAT[1],
            fontsize=13, ha="left", va="bottom")
    ix.text(float(kg) + 0.10, -0.36, "κ_g = %.5f: clock-shift mode,\nnot physical, removed exactly"
            % float(kg), color=AXC, fontsize=13, ha="left", va="top")
    ix.set_title("zoom on the real axis near 0", fontsize=13, color=AXC, loc="left", pad=4)

    # right bottom: the running total of the turning of D along the upper path
    bx = fig.add_axes([0.58, 0.115, 0.39, 0.43])
    s = [0.0]
    for a, b in zip(pts[:-1], pts[1:]):
        s.append(s[-1] + abs(b - a))
    n_arc = sum(1 for p in pieces if abs(complex(p["a"]).real) > 1e-9)  # pieces before the imaginary axis
    y = [0.0] + [float(m / pi) for m in mids]
    bx.plot(s, y, color=ACC, lw=2.4, zorder=3)
    bx.axhline(1, color=AXC, lw=1, ls=(0, (2, 3)), zorder=1)
    bx.axvline(s[n_arc], color=AXC, lw=0.8, zorder=1)
    top = max(y)
    bx.text(s[n_arc] / 2, top * 1.10, "along the arc", color=AXC, fontsize=13, ha="center")
    bx.text((s[n_arc] + s[-1]) / 2, top * 1.10, "down the imaginary axis", color=AXC, fontsize=13, ha="center")
    bx.set_ylim(-0.05, top * 1.24)
    bx.set_xlim(0, s[-1])
    bx.set_xlabel("distance travelled along the counting path")
    bx.set_ylabel("net turning of detA·D\n(half-turns)")
    vm, vr = ball(verdict["PM1"]["delta_over_pi"])
    assert abs(vm - total) <= vr + total_r  # the ledger sum agrees with the sealed verdict
    bx.annotate("ends at %s\n± %s\n= one zero inside" % (vm, sci_up(vr, 2)),
                (s[-1], float(total)), xytext=(s[n_arc] * 0.6, 2.3), color=ACC, fontsize=13.5, ha="center",
                arrowprops=dict(arrowstyle="-", color=ACC, lw=1.2))
    bx.text(s[-1] * 0.995, top * 0.99, "it swings round\nmany times on the\nway; only where it\nends counts",
            color=AXC, fontsize=13, va="top", ha="right")
    save(fig, out, "one-way-to-fall.svg")
    return dict(n_pieces=len(pieces), total=str(total)[:8], total_r=float(total_r), verdict=str(vm) + " +/- " + str(vr),
                top=top, kstar=str(kstar),
                kg=str(kg), gamma=gamma, bands=bands, n_arc=n_arc, vertex=pts[0].real)


# ---------------------------------------------------------------- figure 2: two proofs, one answer
def e5_values(root):
    t = rd(root, "results/e5/runs6/target.log")
    out = {}
    for name in ("V0", "W0"):
        m = re.findall(r"\b%s: K = (\[[^\]]+\])" % name, t)
        assert len(m) == 1, name
        out[name] = ball(m[0])
    return out


def b5_values(root):
    d = json.loads(rd(root, "code/b5/results_summary.json"))["EC"]
    out = {}
    for name in ("V_s", "sig0"):
        lo, hi = ball(d[name]["lo"]), ball(d[name]["hi"])
        out[name] = (lo[0] - lo[1], hi[0] + hi[1])  # outward
    return out


def shared_digits(a, b):
    sa, sb = str(a), str(b)
    n = 0
    for ca, cb in zip(sa, sb):
        if ca != cb:
            break
        n += 1
    frac = sa[:n].split(".", 1)
    return sa[:n], (len(frac[1]) if len(frac) > 1 else 0)


def fig_agree(root, out):
    e5, b5 = e5_values(root), b5_values(root)
    rows = [("sonic velocity V₀", e5["V0"], b5["V_s"], Decimal("1e-20"), "10⁻²⁰"),
            ("central density W₀", e5["W0"], b5["sig0"], Decimal("1e-18"), "10⁻¹⁸")]
    fig = new_fig(7.2)
    title(fig, "Two proofs, one answer",
          "Two programs, the second written without seeing the first, each trap the solution in a box")
    info = []
    for i, (lab, (em, er), (blo, bhi), unit, ulab) in enumerate(rows):
        ax = fig.add_axes([0.05, 0.50 - i * 0.37, 0.90, 0.15])
        bmid = (blo + bhi) / 2
        lo, hi = float((em - er - bmid) / unit), float((em + er - bmid) / unit)
        ax.add_patch(Rectangle((lo, -0.5), hi - lo, 1.0, facecolor=CAT[0], alpha=0.30, edgecolor=CAT[0], lw=1))
        ax.axvline(0, ymin=0.05, ymax=0.95, color=CAT[1], lw=2.4)
        span = hi - lo
        ax.set_xlim(lo - 0.06 * span, hi + 0.62 * span)
        ax.set_ylim(-1.0, 1.0)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.grid(False)
        ax.set_xlabel("distance from the second proof's value, in units of %s" % ulab, fontsize=15)
        pre, nd = rounded_common([em - er, em + er, blo, bhi])
        ax.text(lo - 0.06 * span, 1.35, lab + ":  both boxes round to  %s  (%d decimals)" % (pre, nd), fontsize=16,
                color=TXT, va="bottom")
        ax.text(lo + 0.02 * span, 0, "first proof (E5):\nbox of width %s" % sci_up(2 * er), color=TXT,
                fontsize=13, va="center", ha="left")
        ax.text(hi + 0.03 * span, 0, "second, blind proof (B5):\nbox of width %s,\nfar thinner than the orange line"
                % sci_up(bhi - blo), color=CAT[1], fontsize=13, va="center", ha="left")
        info.append(dict(name=lab, e5_mid=str(em), e5_rad=str(er), b5_lo=str(blo)[:45], b5_hi=str(bhi)[:45],
                         prefix=pre, decimals=nd, b5_width=float(bhi - blo)))
    save(fig, out, "two-proofs-agree.svg")
    return info


# ---------------------------------------------------------------- figure 3: the shape of the solution
class _NoClasses(pickle.Unpickler):
    def find_class(self, module, name):  # the cover holds only dicts, lists, tuples, strings and floats
        raise pickle.UnpicklingError("refusing class %s.%s" % (module, name))


def load_cover(root):
    rel = "results/pickles/lin_fw4_bg3.pkl"
    data = open(os.path.join(root, rel), "rb").read()
    pin = [l.split()[0] for l in rd(root, "MANIFEST.sha256").splitlines() if l.strip().endswith(rel)]
    assert pin and hashlib.sha256(data).hexdigest() == pin[0], "pickle does not match MANIFEST.sha256"
    return _NoClasses(io.BytesIO(data)).load()


def fig_profile(root, out):
    d = load_cover(root)
    e5 = e5_values(root)
    pts = []  # (x, (N, om, V) as (mid, rad))
    for x_mid, _h, _bx, st in d["tau_rec"]:
        pts.append((float(st[3][0]), [(float(m), float(r)) for m, r in st[:3]]))
    for x0, _h, _bx, st in d["rec"]:
        pts.append((x0, [(float(m), float(r)) for m, r in st[:3]]))
    pts.sort(key=lambda p: -p[0])
    xs = [p[0] for p in pts]
    V = [p[1][2] for p in pts]
    om = [p[1][1] for p in pts]
    maxrad = max(max(v[1] for v in V), max(o[1] for o in om))
    # V changes sign between consecutive certified boxes (both boxes strictly signed)
    cross = [(xs[i + 1], xs[i]) for i in range(len(xs) - 1)
             if V[i][0] - V[i][1] > 0 and V[i + 1][0] + V[i + 1][1] < 0]
    assert len(cross) == 1
    v0 = float(e5["V0"][0])

    fig = new_fig(8.2)
    title(fig, "The shape of the critical collapse",
          "Fluid velocity and matter density of the certified solution, from the sonic point inward")
    ax = fig.add_axes([0.10, 0.50, 0.86, 0.33])
    bx = fig.add_axes([0.10, 0.11, 0.86, 0.33], sharex=ax)
    ax.plot(xs, [v[0] for v in V], color=CAT[0], lw=2.4)
    ax.plot([0], [v0], "o", ms=10, color=ACC, zorder=5)
    ax.axhline(0, color=AXC, lw=0.8)
    ax.set_ylabel("fluid velocity V")
    ax.set_ylim(min(v[0] for v in V) * 1.25, v0 * 1.45)
    ax.annotate("sonic point (x = 0): the fluid moves at the\nspeed of sound relative to the shrinking pattern",
                (0, v0), xytext=(-0.12, v0 * 1.22), color=TXT, fontsize=14, va="center", ha="right",
                arrowprops=dict(arrowstyle="-", color=AXC, lw=1.2, shrinkB=7))
    xc = (cross[0][0] + cross[0][1]) / 2
    ax.annotate("V changes sign exactly once", (xc, 0), xytext=(-2.2, 0.045), color=TXT, fontsize=14,
                va="center", arrowprops=dict(arrowstyle="-", color=AXC, lw=1.2, shrinkB=3))
    plt.setp(ax.get_xticklabels(), visible=False)
    bx.plot(xs, [o[0] for o in om], color=CAT[1], lw=2.4)
    bx.set_ylabel("scaled matter density ω")
    bx.set_xlabel("similarity coordinate x = ln(−r/t)   (centre at x → −∞)")
    for a in (ax, bx):
        a.axvline(0, color=AXC, lw=0.8, ls=(0, (2, 3)))
    bx.set_xlim(min(xs) - 0.05, 0.08)
    bx.text(min(xs), max(o[0] for o in om) * 0.93, "← toward the centre", color=AXC, fontsize=14, va="center")
    save(fig, out, "shape-of-the-collapse.svg")
    return dict(n_points=len(pts), x_min=min(xs), x_max=max(xs), maxrad=maxrad, cross=cross[0],
                V_min=min(v[0] for v in V), om_max=max(o[0] for o in om), v0=v0)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        sys.exit(__doc__)
    root = args[0]
    here = os.path.dirname(os.path.abspath(__file__))
    out = here if os.path.basename(here) == "figures" else os.path.join(here, "figures")
    os.makedirs(out, exist_ok=True)
    info = dict(count=fig_count(root, out), agree=fig_agree(root, out), profile=fig_profile(root, out))
    if "--info" in sys.argv:
        print(json.dumps(info, indent=1, default=str))


if __name__ == "__main__":
    main()
