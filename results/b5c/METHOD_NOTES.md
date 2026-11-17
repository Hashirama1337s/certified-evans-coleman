# B5 method notes (working file; RESULT.md is the deliverable)

## Equations (derived, not copied)
- derive.py: G_{mu nu} = 8 pi T_{mu nu} and div T = 0 for ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2
  (Evans-Coleman gr-qc/9402041 eq. (1); KHA gr-qc/9503007 eq. (16)), perfect fluid p = k rho,
  u^t = W/alpha, u^r = W V/a (KHA sec. 3).  Self-similar ansatz x = ln r - ln(-t), rho = Om/(4 pi r^2).
- Variables: mu = 2m/r = 1 - 1/a^2, Om = 4 pi r^2 rho, V, z = a xi/alpha (xi = r/(-t)).
  z is the speed of the homothetic lines relative to the Eulerian (polar-slicing) observers;
  z = Maison's B, Om = Maison's rho-tilde, A_Maison^2 = 1 - mu, N_KHA = 1/z, A_KHA = 1/(1-mu), omega_KHA = Om/(1-mu).
- system.py checks: Hamiltonian + momentum constraint give mu = 2 Om (1 + kV^2 + (1+k)V/z)/(1-V^2) == Maison eq. (8);
  z_x/z == Maison eq. (4) == KHA (18b); KHA (18d) satisfied identically; constraint propagation residual 0.
  KHA (18c) as transcribed from the PDF image is NOT satisfied (unresolved transcription/misprint question; my system is
  confirmed independently by the exact flat-FRW solution for k = 1/3, 1/2, 3/5, frw_check.py).
- system2.py: Om_x = Om nO/(S Q), V_x = (V^2-1) nV/((k+1) S Q), z_x = z P/Q, with
  S = k(1+Vz)^2 - (V+z)^2 (sonic factor: relative speed fluid vs homothetic line = sqrt k), Q = -z(1-V^2)(1-mu).
  Desingularised polynomial field (dx/dtau = S Q):  z' = z P S, Om' = Om nO, V' = (V^2-1) nV/(k+1).
  Centre form Om = z^2 sig, V = z U; U(0) = -2/(3(1+k)) forced, sig(0) = sig0 free (Fuchs exponents 0 and -3).
  Far form N = 1/z: regular at N = 0 (t = 0).
- sonic_sym.py: sonic line e(zs) = (zs, Om_s(zs), V_s(zs)), V_s = (sqrt k - zs)/(1 - sqrt k zs); identities
  S(e)=0, nV(e)=0, J t = 0 (t = de/dzs) proved by cancelling to a zero numerator: a line of equilibria, J has eigenvalue 0.

## Rigorous pieces
- rigode.py: Chebyshev models (cheb.py) + sup-norm Picard contraction, singular K3 operator at the centre
  (u = U0 + x^-3 int_0^x s^2 G), rational right-hand sides via exact polynomial residual / denominator lower bound.
- sonic.py: parameterization of the analytic invariant curve through e(zs) tangent to the eigenvector of the
  positive eigenvalue; exact (closed-form) eigenvalues/eigenvectors; weighted (|th|^2) Volterra contraction.
- prove.py: Poincare-Miranda on H = C (Gamma(sig0) - c(zs)) in a box.
