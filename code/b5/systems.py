"""Symbolic systems for the validated integrator (built from system2.json, derived in derive.py/system.py/system2.py)."""
import json
import sympy as sp
from rigode import System

d = json.load(open('<LAB>/b5/system2.json'))
k, z, Om, V, sig, U, N = sp.symbols('k z Om V sig U N')
loc = {'k': k, 'z': z, 'Om': Om, 'V': V, 'sig': sig, 'U': U, 'N': N}
E = {nm: sp.sympify(d[nm], locals=loc) for nm in d}
S, P, nO, nV, MU = E['S'], E['P'], E['nO'], E['nV'], E['mu']
FOm = Om * nO / (S * z * P)
FV = (V**2 - 1) * nV / ((k + 1) * S * z * P)
U0 = -2 / (3 * (1 + k))


def centre_system(with_var=True):
    """unknowns (sig, U[, sig_s, U_s]) on [0, z1]; Om = z^2 sig, V = z U."""
    dsig = E['dsig']
    g = sp.cancel(E['dU'] * z)        # z dU/dz = g
    ss, Us = sp.symbols('sig_s U_s')
    G = [dsig, g + 3 * (U - U0)]
    ops = ['int', 'K3']
    us = [sig, U]
    if with_var:
        G += [sp.diff(dsig, sig) * ss + sp.diff(dsig, U) * Us,
              sp.diff(g, sig) * ss + (sp.diff(g, U) + 3) * Us]
        ops += ['int', 'K3']
        us += [ss, Us]
    return System(z, us, k, G, ops)


def regular_system(with_var=True):
    """unknowns (Om, V[, Om_s, V_s]) in independent variable z"""
    Os, Vs = sp.symbols('Om_s V_s')
    G = [FOm, FV]
    us = [Om, V]
    if with_var:
        G += [sp.diff(FOm, Om) * Os + sp.diff(FOm, V) * Vs, sp.diff(FV, Om) * Os + sp.diff(FV, V) * Vs]
        us += [Os, Vs]
    return System(z, us, k, G, ['int'] * len(G))


def far_system():
    """unknowns (Om, V) in independent variable N = 1/z (regular through N = 0, i.e. t = 0)"""
    G = [E['dOm_dN'], E['dV_dN']]
    return System(N, [Om, V], k, G, ['int', 'int'])


def tau_field():
    """desingularised polynomial field (z, Om, V)' with dx/dtau = S*Q"""
    return [z * P * S, Om * nO, (V**2 - 1) * nV / (k + 1)]
