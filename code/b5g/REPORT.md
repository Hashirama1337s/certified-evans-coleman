# B5G REPORT: blind second implementation of the Evans-Coleman (k = 1/3) growing mode (sealed test LIGHTBALL-G2)

Saved by the lead from the blind seat's final message; the harness refused the seat's own write. The lead has edited for clarity and length only; the numbers are unchanged, and every number was checked against result_g2.json and the refine logs.

## In plain words
I wrote my own equations for small spherical perturbations of the k = 1/3 Evans-Coleman self-similar collapse solution and solved them with interval arithmetic. I did not use the first seat's equations or code.

There is exactly one growing rate in the window I searched, apart from a coordinate artefact at rate 1:
- **kappa* = 2.810552548827(2)**, so the critical exponent is **gamma = 1/kappa* = 0.355801922443(1)**. Two runs with different settings agree on every printed digit.
- The required certified box has width 3.8e-7: kappa* is in [2.8105522155761715, 2.8105525970458984].

The zero at rate 1 is a gauge artefact, the mode you get by sliding the clock origin, t -> t + const. Its computed eigenfunction matches the closed form component by component. On [0.40, 2.70] there are no other zeros. After dividing out the rate-1 factor, the function is certified negative on the whole range.

Known-answer control: the same pipeline, run on the flat Friedmann (radiation) universe, brackets that background's hand-derived rate 1. That rate is also the time-shift mode, so this control is weaker than a physical rate would be (see Surprises).

## Brief hygiene and independence
- **Number seen before computing.** The brief gave the literature exponent (about 0.3558), so I expected kappa near 2.81. I did not use it as a fit target or starting guess.
  - I located the root by a coarse scan of F on [0, 4] in steps of 0.2 (cert_A.json key coarse_scan_0_4).
  - The scan shows exactly two sign changes: at 1 (gauge) and in (2.8, 3.0).
  - Plain bisection on [2.8, 3.0] then narrowed the root.
  - Caveat: my first development scan window, [0.4, 3.0], was chosen knowing roughly where the root should be.
- **Forbidden paths: none opened.**
  - In b5/ I read only ARGUMENTS.md, METHOD_NOTES.md, cfg_ec.json, cfg_frw.json, result_ec.json, systems.py and frw_check.py. Not RESULT.md. No web searches.
  - The harness put the lab memory index into my context automatically; it contains no number relevant to this task.
  - One attempted read of a saved tool-output copy was refused by the guard, so nothing was read.
- **Formulation differs from b5.**
  - b5 uses a desingularised polynomial field, Chebyshev models and Poincare-Miranda, and defines no perturbations.
  - Mine: polar-areal slicing with the central-lapse gauge alpha(t, r=0) = 1, and unknowns (delta mu, delta Omega, delta V, eta = delta N/N).
  - Algorithm: interval Taylor models with a causal tail bootstrap, Lohner QR, and analytic Frobenius solutions at the centre and at the sonic point.
- **Reused from b5:** only the certified sig0 = [8.9013232777556435191172160567756856 +/- 5.31e-35] and zs = [0.49718669611129068678126112853988237484 +/- 2.61e-39]. I re-integrated the background myself.

## Method
- **Equations, derived from scratch** (derive_g2.py: sympy Einstein tensor and div T).
  - Metric: ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2. Fluid: p = rho/3.
  - Variables: tau = -ln(-t), x = ln(r/(-t)). Unknowns: mu = 2m/r, Omega = 4 pi r^2 rho, V, N = alpha(-t)/(a r).
  - Hamiltonian: mu_x = 2 Omega (1+kV^2)/(1-V^2) - mu.
  - Slicing: (ln N)_x = (1 - 2mu + (1-k) Omega)/(mu - 1).
  - Euler: linear in the derivatives.
  - Sonic factor: proportional to k(1+Vz)^2 - (V+z)^2, with z = 1/N.
  - The exact flat FRW solution satisfies the system (sympy).
- **Perturbations.**
  - Each field is written as f(x) + e^{kappa tau} delta f(x), so kappa > 0 means growth. Then Y_x = (P0 + kappa P1) Y, integrated in z = 1/N.
  - The momentum constraint is used only as a control.
- **Centre.**
  - In scaled unknowns, z Y' = M(z) Y, where M(0) has eigenvalues {0, -2, -3, -3}, independent of kappa, so there is no resonance.
  - The regular solution is y0 = (2/3, 1, -kappa/(4 sig0), 0).
  - eta(0) = 0 is the gauge condition.
- **Sonic point.**
  - s Y' = M(s) Y with exponents {0, 0, 0, lambda(kappa)}, where lambda = tr M(0).
  - The analytic solutions form the 3-dimensional kernel of M(0).
  - For EC, lambda < 0 is certified on every interval used, so there is no sonic resonance for kappa >= 0. For FRW, lambda = 1 - 2 kappa.
- **Characteristic function.**
  - F(kappa) = det[Y_centre(zm), Y_sonic^(1..3)] / (c_p w_p^2), by Leibniz expansion.
  - The normalisation makes F independent of which pivot builds the kernel basis. During development a pivot switch had flipped the sign and created a fake zero near 1.7.
  - F' comes from the kappa-derivative system.
- **Gauge bookkeeping.**
  - Polar slicing with alpha(t, 0) = 1 leaves only t -> t + const. That gives a single gauge mode at kappa = 1, with eigenfunction e^tau (mu_x, Omega_x, V_x, (ln N)_x + 1), so F(1) = 0 exactly.
  - Reparametrisations with kappa != 1 are excluded by eta(0) = 0.
- **Rigour.**
  - python-flint arb, with flint.ctx.cap = 80 set at module top (series lengths at most 37).
  - Taylor models with exact low coefficients; the remainder is folded into the top coefficient; causal tail bootstrap.
  - Constant terms are dropped only where an exact identity holds (proved in sympy).
  - Lohner QR in the regular region. The sonic slope is checked by Krawczyk.
  - On an interval K, F(K) is enclosed in F(c) + F'(K)(K - c). If that evaluation succeeds, F is continuous on K.

## Controls
1. **Background overlap with b5:** the centre-side and sonic-side enclosures overlap at zm in both runs (for example, run B: [0.179724664526690409566492422684 +/- 6.91e-32]).
2. **FRW background:** matches the closed form V = -(1 - sqrt(1-z^2))/z, and overlaps it at zm.
3. **Gauge eigenfunction:** parallel to the closed form at kappa = 1, for EC and FRW, in both runs.
4. **Momentum constraint residual contains 0:**
   - EC: [+/- 1.98e-18] at kappa = 0.5 and [+/- 2.48e-16] at kappa = 2.0;
   - FRW: [+/- 1.66e-12] at kappa = 0.3.
5. **Generated code:** the arb code matches exact sympy rational evaluation.

## Results (result_g2.json)
- **Item 1.**
  - kappa* in [2.8105522155761715, 2.8105525970458984], width 3.81e-7. F' on the interval is [4047.2 +/- 0.0824], so the zero is simple and unique there. lambda = [-2.647411 +/- 3.49e-7].
  - gamma = [0.3558019 +/- 6.47e-8].
  - Refined, in both runs: **kappa* = [2.8105525488272 +/- 7.32e-14]** (bracket [2.8105525488271268040, 2.8105525488272177534]) and **gamma = [0.35580192244308 +/- 8.03e-15]**.
- **Item 2 (FRW known answer).** The rate 1 (t -> t + eps), derived by hand, is enclosed in [0.99999, 1.00001] in both runs, with certified sign change and unique zero. The eigenfunction check passes.
- **Item 3.** The only zero of the undeflated F on [0.40, 2.70] is kappa = 1 exactly, with enclosure [0.9999, 1.0001] and F' = [-2.6 +/- 0.0777]. It is gauge (time translation).
- **Item 4.** G = F/(kappa - 1) is certified negative on all of [0.40, 2.70], using 124 intervals and 263 evaluations.

**Convergence.** Runs A and B differ in every setting:
- Taylor order: 30/28 vs 36/34.
- Step fractions, h0 and z1: all changed.
- zm: 0.4947 vs 0.4935.
- Precision: 200 vs 256 bits.

Both give the identical bracket at width 9.1e-14. Trusted digits: 13 for kappa* and 12 for gamma. No identification of either number is attempted. Run times: cert 279 / 163 s, refine 41 / 81 s.

## Surprises and deviations
1. **The FRW sonic point is degenerate** (double root of the slope equation; resonant recursion). For the FRW control only, the sonic-side background is the closed form, checked by overlap. The perturbation pipeline is unchanged.
2. **The physical FRW modes sit exactly on sonic resonances.** With lambda = 1 - 2 kappa, the modes kappa = -1/2, -1, -2, ... give lambda = 2, 3, 5, ..., so F cannot detect them. F has a pole at kappa = 0. Rate 1 is the only clean closed-form FRW rate F can see.
3. **UNTESTED curiosity: the EC sonic exponent is exactly affine.** It fits lambda(kappa) = r - (1+r) kappa, with r = lambda_-/lambda_+ (b5's background eigenvalue ratio).
   - lambda(0) = [-0.090105970507942775019 +/- 3.41e-22], equal to b5's ratio.
   - lambda(1) = -1 to 24 digits.
   - FRW fits the same rule with r = 1.
   - Unproven (curio_g2.py); it may be the known kink-exponent relation.
4. **Pitfalls fixed:**
   - finite-difference D_n coefficients blew up (now computed once);
   - arb_mat.det gives a crude bound when an LU pivot ball contains 0 (now Leibniz);
   - direct Taylor stepping wraps near the sonic point (now Lohner QR);
   - a pivot sign flip (now normalised).
5. **The interval F'(K) is overestimated** near kappa = 1, which needed geometric refinement.
6. **Checkpointing:** none needed (every run under 5 min). No item is float-only.
