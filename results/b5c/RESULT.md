# LIGHTBALL B5: blind second implementation (existence of the k = 1/3 CSS solution)

## Disclosure (read first)
- **Targets seen.** The brief gave Evans-Coleman's m/r = 0.0596. I also opened Evans-Coleman gr-qc/9402041 (p. 4), which
  quotes n ≈ 1.1485, a ≈ 1.07, m/r → 0.0596 and Ω → 9.56e-3. ROADMAP-RAID.md quotes κ = 2.81055255 and the gauge mode
  0.35699; KHA and Maison quote β ≈ 0.3558. None of these numbers was used as an input, a fit target or a starting
  guess. The starting guess came from my own float scan (σ0 ≈ 8.9, z_s ≈ 0.497). No functional form was fitted.
- **Forbidden paths.** I did not open anything under engine/, runs*/, theory/, DESIGN.md, PROPOSAL-V2/V3.md, LOG.md,
  BAR-LB-E*, AMEND-LB-E-*, seats/THEORY-B2.md or lab memory. An `ls` of <LAB>/ showed those file *names*
  (only names, no contents). One `ls` of my own harness task folder printed other sessions' output file names (none
  opened, none from LIGHTBALL).
- I did not use identify.py: no constant was identified. Nothing was published.

## Object and variables (my formulation)
Spherical Einstein-Euler, p = kρ, k = 1/3, metric ds² = -α²dt² + a²dr² + r²dΩ² (EC eq. (1), KHA eq. (16)). Self-similar
in x = ln r - ln(-t). My unknowns are Maison's variables, not KHA's:
- z = a ξ/α (ξ = r/(-t)), the speed of the homothetic lines relative to the Eulerian (polar-slicing) observers. This equals
  Maison's B and 1/N_KHA.
- Ω = 4πr²ρ (Maison's ρ̃).
- V = the Eulerian 3-velocity. This is the same V as KHA's.
- μ = 2m/r = 1 - 1/a², eliminated by the constraint.

z is used as the **independent variable** (z_x > 0 is certified along the solution), not x. That removes the gauge freedom
(t → ct shifts x) entirely: all reported numbers are gauge invariants.

**Derivation** (derive.py, system.py, system2.py; sympy, from G = 8πT and ∇T = 0 with no copied equations). Checks against
the papers:
- The Hamiltonian and momentum constraints combine to μ = 2Ω(1 + kV² + (1+k)V/z)/(1-V²). This is exactly Maison
  gr-qc/9504008 eq. (8) (A² = 1-μ, B = z, ρ̃ = Ω).
- z_x/z is exactly Maison eq. (4) and KHA eq. (18) second line.
- KHA eq. (18) fourth line holds identically on my system.
- **Unresolved.** KHA eq. (18) third line, as I read it from the PDF image, is NOT satisfied by my system, and I could not
  find a misreading that fixes it. My system is independently confirmed in two ways:
  1. The exact flat-FRW solution satisfies it identically for k = 1/3, 1/2, 3/5 (frw_check.py).
  2. The constraint propagates (residual 0).
- The reduced form is Ω_x = Ω nO/(SQ) and V_x = (V²-1) nV/((k+1)SQ), with:
  - z_x = zP/Q;
  - sonic factor S = k(1+Vz)² - (V+z)², i.e. the fluid speed relative to the homothetic lines equals √k;
  - Q = -z(1-V²)(1-μ).
- Desingularised polynomial field (dx/dτ = SQ): z' = zPS, Ω' = Ω nO, V' = (V²-1) nV/(k+1).
- Sonic line e(z_s) with V_s = (√k - z_s)/(1 - √k z_s). sympy proves exactly that it is a line of equilibria: J has the
  exact eigenvalue 0 along the tangent (sonic_sym.py).
- Centre: Ω = z²σ, V = zU. U(0) = -2/(3(1+k)) is forced, and σ(0) = σ0 is the free parameter (Fuchs exponents 0 and -3).
- t = 0 corresponds to N = 1/z → 0, which is a **regular** point of the system.

## Method (different from Taylor-model/Lohner shooting)
1. **Centre to matching point** (rigode.py, cheb.py): Chebyshev models (polynomial in the Chebyshev basis plus a sup-norm
   remainder, all arithmetic in arb balls). Each step is proved by a **sup-norm Picard contraction** (a radii-polynomial /
   Newton-Kantorovich-type fixed point).
   - The defect is computed exactly: rational right-hand sides are handled by an exact polynomial residual divided by a
     certified lower bound of the denominator.
   - The Lipschitz constants come from interval Jacobians.
   - The first step [0, 0.05] uses the singular Volterra operator U = U0 + z⁻³∫₀^z s²(...) ds (norm 1/3), so the regular
     centre is built in with no series or asymptotics.
   - A variational run over the whole σ0-box gives the mean-value Jacobian.
2. **Sonic point** (sonic.py): the analytic curve through e(z_s) is computed by the **parameterization method in the
   desingularised field**, F(P(θ)) = λθP'(θ).
   - Exact closed-form eigenvalues and eigenvectors are used, plus a 50-term head.
   - The remainder is bounded by a **weighted (|θ|²) Volterra contraction in eigen-coordinates**.
   - z_s is handled by 128 validated sub-boxes. No ODE integration is done on the sonic side.
3. **Matching at z_m = 0.4945** (a different matching point from shooting to the sonic point at x = 0). Existence comes from
   the **Poincaré-Miranda theorem** applied to H = C(Γ(σ0) - c(z_s)) on a box, with C an approximate inverse Jacobian.
4. **Outer part**: the other half of the same analytic curve is followed to z = 3, then to N = 1/z = 0 (t = 0) with the
   same validated integrator.

## Controls
- **Positive control (exact flat FRW, k = 1/2 > 1/3), same pipeline** (result_frw.json). Miranda HOLDS. The certified boxes
  contain the exact values:
  - σ0 = 8/27, box [0.2962962962962962962962962963 ± 1.38e-29];
  - z_s = 9√2/14;
  - V_s = -2√2/5.
  - Also Ω_s = 0.48 and 2m/r = 0.32 to 33 digits. The FRW sonic point is a certified node (λ = 0.2041870402855547551767,
    0.1633496322284438041414).
  - V has 0 zeros, as it should.
- **Positive twin P3**: the k = 1/2 centre-regular solution with σ0 = 8/27 contains the FRW closed form at z = 0.6, with
  residual [±9.26e-79], [±9.99e-78].
- **Nulls** (result_null.json; each must give a certified nonzero matching residual):
  - N1, off-solution (σ0 × (1+1e-10)): residual (-7.8607443e-12 ± 4e-20, 1.1772641e-11 ± 4e-19). Nonzero, a strict catch.
  - N2, k = 17/50 equations at the k = 1/3 solution: (0.0057762592, -0.016812277). Nonzero.
  - N3, FRW closed form against the k = 49/100 solution: (0.00032247715, -0.0021857193). Nonzero (P3 is its twin).
- **Checkpoint**: the centre enclosures are cached as arb.str(radius=True) strings. A write-and-reload in the EC run showed
  the rebuilt balls contain the originals (the strings are not byte-identical after re-rounding). `"resume_tested": true`.
- **50-digit cross-check** (independent code path: mpmath Taylor stepping in z, finite-difference eigenvectors, Frobenius
  centre series): see the section at the end.

## Certified results, k = 1/3 (result_ec.json, results_summary.json; enclosures [lo, hi])
Miranda HOLDS on σ0 ∈ σc ± 1e-35, z_s ∈ zc ± 1e-39, so a solution exists in that box. Uniqueness inside the box is not
certified.

**(1) Sonic-point data.**

| quantity | lo | hi | width |
|---|---|---|---|
| z_s (Maison's B) | 0.497186696111290686781261128539882374837 | 0.497186696111290686781261128539882374843 | 5.2e-39 |
| **KHA's V at the sonic point (Eulerian 3-velocity)** | **0.11243940138809833622376302760313511665** | **0.11243940138809833622376302760313511667** | 1.2e-38 |
| Ω_s = 4πr²ρ | 0.17496734385182408831760021052316466859 | 0.17496734385182408831760021052316466861 | 1.3e-38 |
| 2m/r at the sonic point | 0.46277767055006198093880722579836030116 | 0.46277767055006198093880722579836030124 | 6e-38 |
| KHA's N_s = 1/z_s | 2.0113168912632351728022976840069364054 | 2.0113168912632351728022976840069364056 | 4.5e-38 |
| KHA's A_s = a² | 1.8614267225710816421012892674058759908 | 1.8614267225710816421012892674058759912 | 3.2e-37 |
| KHA's ω_s | 0.32568888942306840447866967391091572242 | 0.32568888942306840447866967391091572258 | 1.6e-37 |
| σ0 = lim Ω/z² at the centre | 8.9013232777556435191172160567756855 | 8.9013232777556435191172160567756857 | 1.1e-34 |

**(2) t = 0 (far side, N = 1/z = 0).**

| quantity | lo | hi | width |
|---|---|---|---|
| m/r | 0.05945140493553022004816764060 | 0.05945140493553022004816764140 | 7.8e-28 |
| Ω | 0.009173972967152708160094838143 | 0.009173972967152708160094838257 | 1.1e-28 |
| V | 0.8968377749150561991476567348 | 0.8968377749150561991476567352 | 4.0e-28 |
| a | 1.065339652212510926461788352 | 1.065339652212510926461788354 | 9.9e-28 |

Comparison with Evans-Coleman (p. 4):
- m/r: EC quote 0.0596; ours minus EC = -1.48595e-4 (0.25%).
- a: EC quote ≈ 1.07; ours rounds to 1.07, consistent.
- Ω: EC quote 9.56e-3; ours is 9.174e-3, about 4% lower.

**(3) Sonic-point type: SADDLE** of the desingularised field (certified c2 = λ₊λ₋ < 0; the third eigenvalue is exactly 0
along the sonic line).
- λ₊ ∈ [0.16429349223752482469695388538395, 0.16429349223752482469695388538405]
- λ₋ ∈ [-0.01480382456620133709303399554384, -0.01480382456620133709303399554376]
- ratio λ₋/λ₊ ∈ [-0.090105970507942775018659844809, -0.090105970507942775018659844807]

The solution enters along the λ₊ (unstable-in-τ) curve. In x it is attracting from the inner side; the other transverse
direction repels weakly. The classification does not depend on the positive desingularising factor.

**(4) Zeros of V: exactly 1 on the inner part** (centre to sonic point). The proof has three pieces:
- V = zU < 0 on (0, 0.05] (U certified negative);
- exactly one simple zero inside a single validated step, where V is monotone and the end signs differ;
- V > 0 on the rest, including the sonic piece, where V ∈ [0.11 ± 4.9e-3].

On the outer part V > 0 all the way to t = 0. The physical checks also pass (1 - μ > 0, P < 0 so z is monotone, V² < 1),
inner and outer.

Convergence: every number above is a rigorous enclosure, so all printed digits are certified. The independent 50-digit
Taylor computation agrees: see the cross-check section.

## Surprises and doubts about definitions
1. **Saddle vs node.** Maison (p. 3) writes that the sonic point "changes from an attractive node to an attractive focus"
   at k_max. The Living Review (2507.07636, l.1286-1287) says it "becomes a nodal point for k > 0.89, rather than a focal
   point". These two sources contradict each other. My certified solution passes through a **saddle**, close to the
   saddle-node point of the sonic line (z ≈ 0.51, where λ₋ crosses 0). The sonic line of the k = 1/3 system has saddle
   (z_s < ≈0.51), node (≈0.51-0.58), focus (≈0.6-0.85) and node segments (float scan, not certified).
2. **Is this "the" Evans-Coleman solution?** m/r differs from EC's 0.0596 by 0.25%, and Ω differs from EC's 9.56e-3 by 4%.
   A centre-regular family also reaches the **node** segment in a C¹ (generally non-analytic) way. A float scan cannot
   tell analytic node solutions from the rest. A 256-bit Newton search for a weak-branch analytic solution near the scan's
   candidates (z_s ≈ 0.546, σ0 ≈ 9.7) did not converge in 8 iterations (residual stalled at ~1e-6; hp13_weak.log). That is
   evidence, not proof, that no analytic node solution sits there. I have **not** excluded other analytic solutions with
   one zero of V. If one exists, "analytic + regular + one zero" does not pick a unique object.
3. KHA eq. (18) third line does not match my derived system (see above).
4. Failed attempts are kept (run_*_attempt*.out):
   - Miranda at 1e-36/1e-38 boxes failed because of preconditioned-face wrapping. Fixed by doing the scalar product before
     the interval product.
   - The EC box was off-centre by 6e-35 in σ0. Fixed by Newton re-centring (refine.py).
   - The first FRW null at k = 51/100 hit a focus, and the second hung. It was replaced by the closed-form null.

## Files (all in <LAB>/b5/)
- derive.py, system.py, system2.py, sonic_sym.py: derivation and identities.
- cheb.py, rigode.py, systems.py, sonic.py, prove.py: the rigorous pipeline.
- cfg_ec.json, cfg_frw.json; result_ec.json, result_frw.json, result_null.json, results_summary.json; run_*.out logs.
- numerics.py, series.py, crosscheck50.py, crosscheck50_ec.json: the independent 50-digit check.
- locate_float.py, branch_scan*.py, solve_hp.py, refine.py: exploration and box centring (not part of the proof).

## 50-digit cross-check (crosscheck50.py, crosscheck50_ec.json; non-rigorous, independent code path)
This is mpmath at 60 digits: Taylor stepping in z (60 terms, step at 0.2 of the radius), a power-series sonic branch with
finite-difference eigenvectors, a 60-term Frobenius centre series, and mpmath odefun on the far part.
- Integrating the sonic branch at the certified z_s inward to z = 0.02 and comparing with the centre series at the
  certified σ0 gives a residual of (Ω: -3.8e-35, V: -1.1e-31). The inward integration amplifies errors along the z⁻³
  centre mode, and the stored σ0/z_s strings carry about 35 and 39 digits.
- t = 0 values: m/r = 0.059451404935530220048167640937…, Ω = 0.0091739729671527081600948381807…,
  V = 0.89683777491505619914765673504…. All three lie **inside** the certified enclosures, at 6e-29, 2e-29 and 5e-29 from
  the midpoints (xc_check.py).
- So the independent numerics confirm about 28 significant digits of m/r. That is the limit set by the certified width,
  7.8e-28.

## Proof-chain notes (what is certified, and what is argued)
- **Certified by computation (arb ball arithmetic):**
  - the Picard contractions on every step (radii r, contraction factor ≤ 0.5);
  - the sonic-branch remainder contraction;
  - θ* by interval Newton;
  - the Poincaré-Miranda sign conditions;
  - the eigenvalue signs;
  - the zero count of V;
  - the physical inequalities.
  - box_cover.py checks that the σ-ball used for the variational run, and the mean-value term, cover the Miranda faces.
- **Argued, not computed:**
  - (a) The contraction fixed point on the centre step is bounded at z = 0 with U(0) = -2/(3(1+k)). The other Fuchs
    exponent is -3, so the bounded solution is the regular, analytic one.
  - (b) The curve P(θ) consists of backward-τ orbits converging to e(z_s) at rate e^{λ₊τ}, tangent to the λ₊ eigenvector.
    Because λ₊ > 0 ≥ the other eigenvalues and nλ₊ never equals another eigenvalue, it is the (unique, analytic) strong
    unstable manifold of e.
  - (c) dx/dθ = S Q/(λθ) is analytic and nonzero (S/θ ∈ [-0.25 ± 7e-3] on the sonic piece). So the solution is analytic in
    x through the sonic point.
  - (d) Matching at z_m plus uniqueness of the regular initial-value problem glues the two pieces.
- **Not done:** uniqueness of the solution inside the box; exclusion of other analytic solutions with one zero of V
  (see Surprises 2); perturbation modes (not part of B5).
- **Compute**: each proof run took about 55 minutes of one core at 320/600-bit precision. Two runs went in parallel under
  memcap 8 GB, with a peak of 0.08 GB.
