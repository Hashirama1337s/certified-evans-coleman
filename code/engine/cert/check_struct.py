"""exact check: structured sonic field == expanded polynomial field (fields.pkl)."""
import sys
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import system as S0
from pipeline import load_expr, FIELDS
Ms, bs = S0.Msym, S0.bsym
det = Ms[0, 0]*Ms[1, 1] - Ms[0, 1]*Ms[1, 0]
gst = [det*S0.gN, Ms[1, 1]*bs[0] - Ms[0, 1]*bs[1], Ms[0, 0]*bs[1] - Ms[1, 0]*bs[0], det]
ps = [load_expr(s) for s in FIELDS['pol_s']]
print('structured - expanded == 0:', [sp.simplify(sp.cancel(sp.together(a - b))) == 0 for a, b in zip(gst, ps)])
