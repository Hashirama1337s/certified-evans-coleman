# A certified Evans-Coleman critical solution: existence, one growing spherical mode, and its exponent

**Moki & Julio**

Version 1.0, 2026. Every number below is traced to a lab result file in `CLAIMS.md`. Balls are
printed as `[mid ± rad]` exactly as the certifying program stored them; interval ends are printed as stored.

---

## Abstract

We study the self-similar collapse of a self-gravitating perfect fluid with the radiation equation of state p = ρ/3 in
spherical symmetry, at the threshold of black-hole formation: the Evans-Coleman critical solution. With interval (ball)
arithmetic we give a computer-assisted certificate that a solution with its defining properties exists: regular at the
centre, real-analytic through a sonic point of saddle type, with exactly one zero of the fluid velocity inside the sonic
point, no second sonic point outside it, and a continuation to t = 0. In the variables of Koike, Hara and Adachi (KHA) its
sonic-point velocity is V0 ∈ [0.1124394013880983362 ± 3.31e-20]. A second, blind implementation with a different
formulation (Chebyshev models and the Poincaré-Miranda theorem in place of Taylor models and a Krawczyk test) proves the
existence of a solution whose sonic velocity agrees with the first to all 19 certified digits; local uniqueness is not
claimed by either. On the first enclosure we study linear spherical perturbations that are analytic through the sonic
point and regular at the centre. An exact gauge deflation removes the coordinate mode at κ_g ∈ [0.3556992037 ± 7.41e-11].
An argument-principle count on a half-disc of radius 40, together with an a-priori bound ("the Far Wall") that excludes
growing modes for every |κ| ≥ 40, shows that there is exactly one growing non-gauge mode in Re κ ≥ 0. It is real and
simple, κ* ∈ [2.810552549 ± 3.83e-10], so the exponent of the KHA/Maison scaling argument is
γ = 1/κ* ∈ [0.3558019224 ± 6.97e-11]. The first engine already agrees with the 9-digit value 2.81055255 printed in the text of Hara, Koike and Adachi (HKA,
1996): its enclosure [2.8105525486170, 2.8105525493830] rounds to it. A blind second computation, in a different gauge and
with different code, on the B5 solution encloses a real growth rate in [2.8105525488272 ± 7.32e-14] (γ ∈
[0.35580192244308 ± 8.03e-15]), inside the first enclosure; this 13-digit refinement is certified, beyond the sealed bar.
If the two solutions coincide (not claimed; see Section 11), HKA's table value 2.8105525488 is correct to its last digit
and two further digits are certified. We do not prove the critical-collapse mass law, nonlinear stability, or anything
about non-spherical perturbations. Limitations are listed in Section 11.

---

> **In plain words**
>
> When a cloud of light (radiation) collapses under its own gravity, a small push decides between two fates: it forms a
> black hole, or it disperses. Exactly at the dividing line, computer simulations since 1994 show a special collapse
> that keeps the same shape while it shrinks, a "ball of light" balanced on a knife-edge. Physicists have computed it
> numerically for thirty years. Nobody had a checkable proof that it exists, or that it has exactly one way to fall off
> the knife-edge.
>
> We give both, with a computer that carries a guaranteed error bar on every number. Two separate programs, built in
> different ways, the second written without seeing the first, each prove that the solution exists, and they agree to 19
> digits. Then we count
> how small disturbances grow. Exactly one kind grows, at a rate whose inverse is 0.3558019224, the famous critical
> exponent found numerically in 1995. To make the count finite we prove a "far wall": no disturbance that oscillates or
> grows very fast can survive, all the way to infinity.
>
> The idea in one picture: a ball balanced on a ridge. It has exactly one direction in which it rolls off. Where the
> picture breaks: our ball is not sitting still. It shrinks forever while keeping its shape, so "rolling off" means
> "growing relative to the shrinking," and we only test round (spherical) pushes.

---

## 1. Results

Notation (Section 2): KHA variables (N, A, ω, V) as functions of x = ln(−r/t); the sonic point sits at x = 0; linear
perturbations behave as e^{κ s} with s = −ln(−t). "Certified" means proved by a computer with outward-rounded ball
arithmetic (python-flint arb/acb), inside a sealed test (Section 9). "Proved" means a written argument checked by at least
one internal hostile referee seat.

**Theorem 1 (existence, first implementation; CERTIFIED, sealed test E5).** For k = 1/3 the self-similar Einstein-Euler
system in KHA variables has a solution that is regular at the centre (x → −∞), real-analytic through the sonic point
x = 0, and continues to x → +∞ (t = 0, r > 0), with:
- sonic data V0 ∈ [0.1124394013880983362 ± 3.31e-20], centre density parameter W0 ∈ [8.90132327775564352 ± 7.92e-18]
  and centre-chart shift xc ∈ [0.212365646758945385 ± 7.95e-19];
- a sonic point of saddle type, with eigenvalues λ ∈ [4.29174946397 ± 8.77e-12] and λ_b ∈ [−0.3867122506 ± 4.23e-11];
- exactly one sign change of V between the centre and the sonic point;
- no second sonic point for x > 0 (the sonic discriminant stays negative, closed by an invariant-region lemma beyond a
  finite x);
- m/r at t = 0 in [0.0595 ± 8.24e-5].

Local uniqueness is not claimed.

**Theorem 2 (existence, blind second implementation; CERTIFIED, sealed closure test B5C).** For k = 1/3 a solution
exists that is regular at the centre, crosses a saddle-type sonic point real-analytically, has exactly one zero of V
between the centre and the sonic point, and reaches t = 0. In Maison-type variables (Section 4) its central density
parameter is σ0 ∈ [8.9013232777556435191172160567756855, 8.9013232777556435191172160567756857], its KHA sonic velocity is
V_s ∈ [0.11243940138809833622376302760313511665, 0.11243940138809833622376302760313511667], and at t = 0,
m/r ∈ [0.05945140493553022004816764060, 0.05945140493553022004816764140]. Uniqueness is not claimed. Its agreement with
Theorem 1 (V0, and σ0 against W0) is a POST HOC comparison of two rigorous enclosures, which overlap.

**Theorem 3 (characteristic function and exact gauge deflation; PROVED, given two certified nonzero numbers and a certified sign).** Let
s_κ be the sonic-analytic, constraint-satisfying perturbation normalised by (Ā_p, N̄_p)(0) = (0, 1), and ρ1, ρ2 the two
centre-regular solutions. Then D(κ) = det[s_κ(x_m), ρ1(x_m), ρ2(x_m), e_V] is analytic on Re κ > κ_res, with
κ_res < 0 (d = b_−−(0) < 0 from the certified saddle; detA root [−0.0990290820550 ± 2.91e-14]), and, provided c_V(x_m) ≠ 0 and Ā_ss'(0) ≠ 0, D(κ) = 0 exactly when κ is a physical (non-gauge) eigenvalue.
HKA's normalisation gives D_HKA(κ) = −(κ − κ_g) D(κ)/Ā_ss'(0), so the gauge zero is removed exactly. Both
preconditions are certified: c_V(x_m) ∈ [0.19288155165634 ± 1.41e-15] at x_m = −4.6 and
Ā_ss'(0) ∈ [−0.198928186801772 ± 4.57e-16]. (Section 5; this is Theorem 1.6 of the proof file PROOF-B4.)

**Theorem 4 (the Far Wall; CERTIFIED).** On the Theorem 1 background there is no physical eigenvalue with Re κ ≥ 0 and
|κ| ≥ 40. (Section 7: three bands, [40, 55], [55, 2000] and [2000, ∞).)

**Theorem 5 (the count; CERTIFIED, sealed test M).** In the closed half-plane Re κ ≥ 0, the Theorem 1 solution has
exactly one growing non-gauge linear mode within the class of spherical perturbations that are analytic through the
sonic point and regular at the centre, counted modulo the gauge mode at κ_g ∈ [0.3556992037 ± 7.41e-11]. The zero of D is
simple, it is real, and there is no mode on Re κ = 0. (Section 6.)

**Corollary 6 (κ* and γ; CERTIFIED).** The growing mode has κ* ∈ [2.810552549 ± 3.83e-10]. Hence
γ = 1/κ* ∈ [0.3558019224 ± 6.97e-11]. Here γ is the exponent that the KHA/Maison scaling argument predicts from a
single growing mode; the mass-scaling law itself is not proved here. (Section 8.)

**Remark 7 (kink modes; CITED + CERTIFIED input).** The count does not cover perturbations that are non-analytic at the
sonic point ("kink" modes). Harada's criterion (gr-qc/0109042) says that, in collapse, solutions crossing a saddle sonic
point are stable against the kink mode; Theorem 1 certifies that the sonic point is a saddle. So there is no kink-mode
instability by Harada's criterion. Harada calls this a necessary condition for stability, so we do not write "stable".

| Item | Status | Seal | Section |
|---|---|---|---|
| Theorem 1 | certified, refereed | d0af75803f (HOLD) | 3 |
| Theorem 2 | certified, refereed, blind | 71904980c4 (HOLD) | 4 |
| Theorem 3 | proved, refereed (preconditions certified in every count run) | b1b0ee1cde (preconditions) | 5 |
| Theorem 4 | certified, refereed (three referee rounds) | files pinned under b1b0ee1cde | 7 |
| Theorem 5 | certified, refereed | b1b0ee1cde (HOLD) | 6 |
| Corollary 6 | certified; blind second enclosure HOLD, [2.8105525488272 ± 7.32e-14] | b1b0ee1cde; 8a1c27c7fb | 8 |
| Remark 7 | cited criterion applied to a certified fact | - | 1 |

## 2. Background and conventions

**Metric and variables (KHA eq. (16) and Sec. 3).** ds² = −α²dt² + a²dr² + r²dΩ², perfect fluid with p = kρ, k = 1/3
(KHA write p = (γ−1)ρ with γ = 4/3). Similarity coordinates s = −ln(−t), x = ln(−r/t). KHA's variables are
N = α a⁻¹ e⁻ˣ, A = a², ω = 4πr²Aρ, and V, the fluid 3-velocity seen by Eulerian observers (u^t = Γ/α, u^r = ΓV/a).

**Self-similar system.** For functions of x alone:
- Hamiltonian constraint (KHA (18a)): A_x/A = 1 − A + 2ω(1 + kV²)/(1 − V²);
- slicing condition (KHA (18b)): N_x/N = A − 2 − (1 − k)ω;
- the momentum constraint gives A algebraically (Maison eq. (8)): A = 1 + 2ω(1 + (1+k)NV + kV²)/(1 − V²);
- two fluid equations M(y)(ω_x, V_x)ᵀ = b(y), y = (N, ω, V).

The system is autonomous in x and three-dimensional. It was derived symbolically from the Einstein tensor and ∇·T = 0
(not transcribed from the papers); it reproduces KHA (18a), (18b) and Maison (8) exactly. det M =
−2(1+k)Nω(1−V²)²[k(N+V)² − (1+NV)²], so the **sonic point** is where (1 + NV)² = k(N + V)²: the fluid speed seen from a
line of constant x equals √k. We use the branch N = (1 − V√k)/(√k − V) and fix the residual gauge (t → λt, a shift in
x) by putting the sonic point at x = 0.

**Boundary conditions.** The centre x → −∞ is regular (A → 1, V → 0; KHA condition (ii)). In the blown-up chart
P = 1/N, W = ω/P², U = V/P, u = eˣ the regular centre is Z0 = (0, W0, −2/(3(1+k))), a one-parameter family in the
central density W0; the linearisation there has eigenvalues {1, 0, −3}, and −3 is the singular mode that the matching
must kill. At the sonic point we require real-analyticity (KHA condition (i)). No condition is imposed at t = 0
(x → +∞); we certify instead that the solution has no second sonic point there.

**Perturbations (KHA Sec. 5).** h = h_ss + ε e^{κs} h_p(x), with logarithmic perturbation variables
Y = (Ā_p, N̄_p, ω̄_p, V_p) (Ā = ln A, and so on). The system is a pencil Q(x)Y' = (P0(x) + κP1(x))Y plus an algebraic
identity I(Y) = c(x, κ)·Y = 0 (the linearised momentum constraint). KHA fix the gauge by N̄_p(0) = 0 at the sonic point.
The normalisation of κ: a growing mode has Re κ > 0, and KHA's β and Maison's γ are both 1/κ for the single relevant
mode (seats/GATE-REFRESH, normalisation paragraph).

## 3. Existence, first implementation (Theorem 1)

**Unknowns and map.** X = (V0, W0, xc) ∈ ℝ³. For given X the program builds (i) the analytic manifold leaving the sonic
point, (ii) the regular-centre manifold, and (iii) the flow between them, and returns a matching residual R(X) ∈ ℝ³ in
the centre variables Z = (P, W, U).

**Ingredients.**
1. *Invariant manifolds.* Both the sonic and the centre manifolds are computed by the parameterisation method in
   ball arithmetic, with a rigorous tail bound: if λθg' − μg = f with f = O(θ^{K+1}), then
   sup|g| ≤ sup|f|/|(K+1)λ − μ| whenever K + 1 − μ/λ > 0, plus a Cauchy defect tail and a Neumann bound for the
   nonlinear terms. On the certified box the sonic manifold tail is eps = [1.85525e-32 ± 1.91e-38] and the centre tail is
   eps = [3.28636e-56 ± 2.50e-62] (REFEREE-E5CAP).
2. *Flow.* A Lohner-type integrator with weighted QR (Taylor models with a complex Picard a-priori box and a Cauchy
   remainder), started on a short τ-time leg next to the sonic point (where x is a bad clock) and switched to x-time with
   a certified crossing.
3. *Repair R8 (the step that made the warm-up close).* In the sonic tail bound, the head coefficients are bounded by the
   mean value form |c_n(Ṽ0)| + r·|c_n'(box)| (the derivative is propagated with dual numbers through the sonic point and
   the eigenpair), in place of box evaluation. Box evaluation had booked the width of the coefficients, which is first
   order in the box radius r, as tail error. A referee checked the propagated derivative against central differences at
   all 525 head coefficients (relative error 7.08e-23 at k = 1/3).
4. *Existence test.* A 3×3 Krawczyk test with a scalar error addend E: if K(X) ⊂ int X, a zero exists in X by Brouwer.
   The tail error has no derivative bound, so the test proves existence only, not uniqueness.

**Run (runs6/target.log).** On the sealed radius ladder 1e-15 ... 1e-9 the first rung contracts: max rad(K)/rad(X) =
0.00703, Ebox = [1.93e-24 ± 7.05e-28], PE2 margin [1.385e-27 ± 3.78e-31]. The enclosures in Theorem 1 are the K box.
PE2 holds on 194 τ-boxes and 560 x-boxes: one V sign change, no ambiguous step, the sign of det certified on the inward τ
leg, x monotone on the manifold and centre pieces, and consistent monotone runs. The outward leg (731 steps) has
D_sonic < 0 and V > 0 throughout, and the far-side invariant-region lemma holds; m/r(t = 0) ∈ [0.0595 ± 8.24e-5].

**Controls (all under the same seal).**
- *Positive control (warm-up).* Flat Friedmann at k = 1/2 in the same chart is an exact solution with
  V0 = −2√2/5 and W0 = 8/27 (sealed amendment AMEND-LB-E-1, checked symbolically). Started from float64 seeds, the same
  pipeline contracts at r = 1e-15 (ratio 0.0247, Ebox [1.76e-22 ± 6.63e-26]) and the box contains the exact data:
  V0 ∈ [−0.565685424949238020 ± 5.73e-19], W0 ∈ [0.2962962962962963 ± 2.32e-17].
- *Content nulls.* Boxes of radius 1e-15 centred 1e-13 off the solution in V0, in W0 and in xc must not contract. All
  three refuse through the Krawczyk test itself.
- *Instrument controls.* The tail tests (including an injected 3·eps error that must be flagged, a derivative-versus-
  central-difference check, and a frozen-eigenpair null that must fail), the Lohner thin and box tests, and a structure
  check all pass.

**History (disclosed; Section 9).** Four earlier sealed attempts did not certify: E and E2 FAILED (the k = 1/2 warm-up did
not contract), E3 was first settled HOLD and then corrected to VOID (a library series cap silently truncated the sonic
tail input, so the printed tail 3.7e-59 was really about 3.7e-35), and E4 FAILED (the warm-up error addend was linear in
the box radius at k = 1/2; diagnosed and fixed by R8). The E5 bar disclosed all earlier numbers; the claim bands came
from published values and exact algebra.

## 4. Existence, blind second implementation (Theorem 2)

A separate seat wrote this proof without opening the first implementation's code, logs or design files (it saw only
their file names), and with no published target used as an input (the brief quoted Evans and Coleman's m/r = 0.0596).

**Formulation.** Maison-type variables: z = aξ/α (ξ = r/(−t); z = 1/N in KHA terms), Ω = 4πr²ρ, the same V, and
μ = 2m/r eliminated by the constraint μ = 2Ω(1 + kV² + (1+k)V/z)/(1 − V²), which is exactly Maison (8). The variable z is
the independent variable (z_x > 0 is certified), which removes the x-shift gauge, so every reported number is a gauge
invariant. The desingularised polynomial field is z' = zPS, Ω' = Ω n_O, V' = (V² − 1)n_V/(k+1) with sonic factor
S = k(1 + Vz)² − (V + z)², and the sonic set is a line of equilibria (proved symbolically).

**Method.**
1. Centre to the matching point z_m = 0.4945 by Chebyshev models with sup-norm Picard contractions; the first step
   [0, 0.05] uses the singular Volterra operator U = U0 + z⁻³∫₀^z s²(...)ds of norm 1/3, which builds in the regular
   centre (Appendix D.1).
2. The analytic curve through the sonic equilibrium by the parameterisation method, F(P(θ)) = λ₊θP'(θ), with a weighted
   Volterra remainder contraction; it is the analytic strong-unstable manifold (Appendix D.2).
3. Existence from the Poincaré-Miranda theorem for H = C(Γ(σ0) − c(z_s)) on the box σ0 ∈ σc ± 1e-35,
   z_s ∈ zc ± 1e-39, with an approximate inverse C.
4. The other half of the same analytic curve to z = 3, then to N = 1/z = 0 (t = 0), which is a regular point.

**Results.** Besides Theorem 2: z_s width 5.2e-39; sonic eigenvalues λ₊ ∈ [0.16429349223752482469695388538395,
0.16429349223752482469695388538405], λ₋ ∈ [−0.01480382456620133709303399554384, −0.01480382456620133709303399554376],
ratio λ₋/λ₊ ∈ [−0.090105970507942775018659844809, −0.090105970507942775018659844807] (a saddle); at t = 0,
Ω ∈ [0.009173972967152708160094838143, 0.009173972967152708160094838257] and
a ∈ [1.065339652212510926461788352, 1.065339652212510926461788354]. Each proof run took about 55 minutes on one core.

**Controls.** Flat Friedmann at k = 1/2 through the same pipeline: Miranda holds and the boxes contain σ0 = 8/27
(box [0.2962962962962962962962962963 ± 1.38e-29]), z_s = 9√2/14 and V_s = −2√2/5. Three nulls give certified nonzero
matching residuals (σ0 × (1 + 1e-10); the k = 17/50 equations at the k = 1/3 solution; the Friedmann closed form against
the k = 49/100 solution). A 60-digit floating-point cross-check by an independent code path lands inside every
certified t = 0 enclosure (not a control, a consistency check).

**Closure test B5C (after the referee, sealed, HOLD).** C1: on the outer sonic stretch θ ∈ [0, θ_o],
θ_o = [0.031099573504505126343 ± 4.02e-22], S/θ ≠ 0, P < 0, μ < 1 and V > 0 hold on 8 pieces with 0 ambiguous. C2: the
Miranda preconditioner has det C ∈ [−81.709996897608312500 ± 6.34e-20]. C3: P/z ∈ [−1 ± 0.0158] on [0, 0.05], over the
whole σ0 box. C4: all 20 cached objects fed to Miranda contain their reloaded copies. C5' (positive control, amended
before the run): the Friedmann k = 1/2 outer leg reaches the chart edge z = 1 with enclosures [0.96 ± 3e-30] and
[−0.8 ± 3e-30] containing the exact Ω = 24/25 and V = −4/5. Negative control: the same V > 0 test, run across the real V
zero of the certified solution, reports a certified sign change. Disclosed gap: the far t = 0 leg has no exact-solution
positive control (flat Friedmann never reaches t = 0 in this chart).

**The 19 digits.** Theorem 1 gives V0 ∈ [0.1124394013880983362 ± 3.31e-20]; Theorem 2 gives
V_s = 0.112439401388098336223763... with width 1.2e-38. The two enclosures overlap, and so do W0 and σ0. The comparison
is POST HOC (B5 had no sealed bar of its own). Since neither implementation proves local uniqueness, we do not claim the
two solutions are the same; we claim two independent proofs of existence of solutions with these properties, whose
enclosures coincide to 19 digits.

## 5. The linear problem, the characteristic function and gauge deflation (Theorem 3)

All statements in this section are from PROOF-B4 (sections 1.1-1.5), refereed with the count (REFEREE-B4FW2,
"what is counted: CORRECT WITH FIXES"; the fixes were the certified preconditions). Theorem 3 also had a line-by-line
stand-alone referee (seats/REFEREE-RELEASE.md, section 2: CORRECT, minor fixes, applied).

**Constraint propagation.** Along solutions of the pencil, dI/dx = c_I(x) I with
c_I = (8NVω + 2V²ω − 3V² + 6ω + 3)/(3(V² − 1)), independent of κ and regular on x ≤ 0, with c_I → −1 at the centre
(exact symbolic check). Hence I ≡ 0 as soon as I vanishes at one point, and every solution bounded at the centre with
Ā_p → 0 has I ≡ 0. The constrained solution space W_κ is three-dimensional.

**Sonic point (regular singular).** The Frobenius exponents at x = 0 are 0, 0, 0 and ρ(κ) = (d − κ)/λ', where
λ' = λ_−'(0) > 0 is the slope of the incoming characteristic speed and d < 0 (on the certified background
λ_−' ∈ [1.0986, 1.0995] and b_−− < 0 next to the sonic point, PROOF-B4-v4 sec. 3). Resonances ρ = n ≥ 0 sit at
κ = d − nλ', all in Re κ < 0. Imposing I = 0 leaves a two-dimensional space S_κ of sonic-analytic solutions,
parametrised by (Ā_p, N̄_p)(0) where a 2×2 normalising determinant is nonzero; its root is certified at
[−0.0990290820550 ± 2.91e-14] < 0 in every count run. Set κ_res := max(d, that root) < 0. Define s_κ ∈ S_κ by
(Ā_p, N̄_p)(0) = (0, 1) and s_κ^HKA by (1, 0); both are analytic in κ on Re κ > κ_res.

**Centre (regular singular in v = p eˣ).** The blown-up system v dY/dv = M_c(v, κ)Y has M_c(0) with eigenvalues
{−1, 0, 0, −2}, independent of κ. The bounded solutions form the two-dimensional regular space R_κ = span(ρ1, ρ2), with
convergent series whose coefficients are polynomials in κ; every other solution blows up through the −2 exponent
(V_p ~ v⁻²) or the −1 exponent (a point mass Ā_p ~ v⁻¹).

**Modes and gauge.** Following KHA's conditions, κ is a physical eigenvalue if S_κ ∩ R_κ contains a solution that is
not a multiple of the gauge solution G_κ = (Ā_ss', N̄_ss' + κ, ω̄_ss', V_ss') (the image of the residual x-shift),
which lies in S_κ ∩ R_κ for every κ. For x > 0 the coefficients are analytic (no second sonic point, Theorem 1), so the
problem on ℝ equals the problem on x ≤ 0.

**Lemma 5.1 (gauge).** Let κ_g := −N̄_ss'(0) = 2 − A0 + (2/3)ω0. Then
(a) G_κ = Ā_ss'(0) s_κ^HKA + (κ − κ_g) s_κ (both sides lie in S_κ with the same data at 0);
(b) if Ā_ss'(0) ≠ 0, G_κ and s_κ are independent for every κ;
(c) κ is a physical eigenvalue if and only if s_κ ∈ R_κ;
(d) for every linear functional Φ that vanishes on R_κ, Φ(s^HKA) = −(κ − κ_g)Φ(s)/Ā_ss'(0).

**Theorem 3, proof.** D(κ) := det[s_κ(x_m), ρ1(x_m), ρ2(x_m), e_V] with e_V = (0, 0, 0, 1). If
a s + b ρ1 + c ρ2 + d e_V = 0 at x_m, applying I (which vanishes on s, ρ1, ρ2) gives d·c_V(x_m) = 0, so d = 0; since ρ1,
ρ2 are independent at every regular point, a ≠ 0 and s ∈ span(ρ1, ρ2) at x_m, hence s ∈ R_κ by uniqueness, and Lemma
5.1(c) concludes. The converse is immediate. Analyticity follows from that of s_κ and the locally uniform convergence of
the centre series. Lemma 5.1(d) applied to D gives the deflation identity. ∎

The undeflated HKA function D_HKA has the gauge zero at κ_g; the deflated D does not. Discreteness: D is not identically
zero (it is nonzero on the whole count contour, Section 6), so its zeros in Re κ > κ_res are isolated and finitely many in
every compact set: the mode spectrum has no continuum. If the zero at κ* is simple, the mode is unique up to scale
(S_κ* = R_κ* contains G and one mode).

**What "modes" means here.** The count is about solutions of the mode problem in the class above (analytic through the
sonic point, regular at the centre). It is not a statement about the full spectrum of a linearised evolution operator on
a function space; absence of residual or continuous spectrum for such an operator is not proved (PROOF-B4 Remark 1.7).

## 6. The count (Theorem 5)

**Function counted.** The program encloses F(κ) = detA(κ)·D(κ), where detA is the affine normalising determinant of
Section 5 (root at −0.0990... < 0, so F and D have the same zeros in Re κ ≥ 0). F is real on the real axis, so
F(κ̄) = conj F(κ).

**Contour.** The half-disc {Re κ ≥ 0, |κ| ≤ 40}, upper half only by conjugate symmetry: the count is
N = Δ arg F(upper path)/π, which needs F real and nonzero at 0 and at the arc vertex on the real axis. The arc is a
polygon of 8 chords with vertices at radius R_p = 40/cos(Δφ/2)·(1 + 1e-12) = 40.1935428950926. In exact rational
arithmetic every chord stays at distance ≥ 40 from 0 (min distance² − 40² = 3.200e-09) and every vertex has |κ| ≤ 55, so
the region between the polygon and the circle |κ| = 40 lies in the Far Wall band [40, 55]. The plan has 99 discs and is
exactly continuous from the vertex to 0.

**Enclosure of F on a disc.** κ-Taylor models of the transition matrices along the background (degree m = 16 in the disc
variable, J = 4, nsub = 2, prec 128 bits), composed from the sonic start to the matching point. The remainder is taken
in the characteristic frame of the large-κ principal part A1 (Section 6.1). A point-value argument sharpening picks the
right branch of arg F (rigorous: the enclosure is narrower than 2π, so at most one translate meets it). Refinement was
fixed in advance (at most 3 halvings; any failure would have been INCONCLUSIVE and the test FAILED).

**Matching point.** Theorem 3 holds for any fixed x_m < 0 at which c_V(x_m) ≠ 0. The design notes first used
x_m = −2. The certified pipeline matches further in, at the end of the certified background tube, x_m = −4.6. That
keeps the centre series accurate at large |κ|, where the centre leg must reach the Bessel zone v_m ≲ 4/|κ|. The
precondition is certified for every x_m in the tube-end box −4.6 ± 3e-14, which contains the point actually used.

**Preconditions (certified at the start of the run).** c_V(x_m) ∈ [0.19288155165634 ± 1.41e-15] (x_m = −4.6),
Ā_ss'(0) ∈ [−0.198928186801772 ± 4.57e-16], detA root [−0.0990290820550 ± 2.91e-14], and the same normalisation index at
all 199 plan points.

**Result (runs7/verdict_M.json).** N(R1 = 40) = 1 with Δ/π ∈ [1.0000 ± 2.42e-5], parity (−1)^N = sign F(0)·sign F(R_p)
true, 99 discs, 0 refinements, 0 inconclusive. Runtime 3786 s.

**From the count to Theorem 5.** F has exactly one zero, counted with multiplicity, in the half-disc, and F ≠ 0 on the
contour, including the imaginary-axis segment. The Far Wall (Theorem 4) excludes zeros with |κ| ≥ 40 in Re κ ≥ 0.
Theorem 3 identifies zeros with physical eigenvalues. So there is exactly one, it has multiplicity one (simple), and it
is real (a non-real zero would come with its conjugate). No zero lies on Re κ = 0, where κ = 0 is the pure-gauge
x-translation.

**Controls (run first, same seal).** C1a: the undeflated D_HKA winds once around κ_g. C1b (null): the deflated D winds
zero times there. C1c: bisection encloses κ_g ∈ [0.3556992037 ± 7.41e-11], containing the closed form 2 − A0 + (2/3)ω0.
Positive control of the exact contour shape: the half-disc at R1 = 3 gives N = 1 with Δ/π ∈ [1.0000000 ± 3.42e-8] and
parity true.

### 6.1 The bridge: a change of frame turned a 21.5-hour count into a one-hour count

The first build enclosed the disc transition matrix as P_0 + Σ_j |P_j| r^j + tail, with absolute values taken entrywise
in the original variables (Ā, N̄, ω̄, V). The remainder then grows like exp(c r) in the disc radius r, which forced disc radii
r = 0.025-0.05 for |κ| > 40, about 1039 discs at R1 = 55 and a projected 21.5 hours on one process. The repair takes the
absolute values after a change of variables per step, T_s = [[I, 0], [−B⁻¹C, U]], built from A1 = [[0, 0], [C, B]] and the
eigenvectors U of B: columns 3-4 are the characteristic directions and columns 1-2 span the kernel of A1, which is the
slaving structure of the Far Wall (Section 7), frozen per step. Any invertible T is rigorous (T is an exact matrix and
T⁻¹ an enclosure; referee REFEREE-B4FW: "LinearTM3 is a valid enclosure"). At 20i with r = 0.1 the remainder ratio fell
from 1.65e-2 to 1.43e-7, imaginary-axis discs grew to r = 0.5, and, with the Far Wall lowered from 55 to 40, the plan fell
to 99 discs. The projection was 1.5-1.7 hours; the sealed count ran in 63 minutes.

## 7. The Far Wall (Theorem 4)

Large |κ| is where a finite computation cannot reach. The Far Wall turns "no eigenvalue for |κ| ≥ 40" into finitely many
interval inequalities that hold uniformly in |κ|. The proof is in PROOF-B4-v4 with two amendments (AMEND1: referee
repairs and the [40, 55] band; AMEND2-TAIL: the tail |κ| ≥ 2000). We write out the main estimates; details are in
Appendices A-C.

**Characteristic form (on x < 0).** With u = (ω̄_p, V_p), right eigenvectors r_± and left eigenvectors l_± of
Λ0 = M_s⁻¹M_x (speeds λ_± = 1 + N(V ± c)/(1 ± cV), c = 1/√3) and amplitudes w_± = l_± u:

  λ_± w_±' = (−κ + b_±±) w_± + b_±∓ w_∓ + g_±,   g_± = e_±A Ā + e_±N N̄,

with Ā fixed by the identity, Ā = −(c_N N̄ + α̃_+ w_+ + α̃_− w_−)/(κ − A), and N̄' = AĀ + H3(w_+ + w_−), N̄(0) = 0 (the
s^HKA normalisation, which excludes the gauge solution). Write ℓ_− = |λ_−| (λ_−(0) = 0 at the sonic point),
T = 1/λ_+ + 1/ℓ_−, β_+ = b_++/λ_+, β_− = b_−−/ℓ_−, σ̂ = b_−+/(ℓ_− T), τ̂ = b_+−/λ_+, K = |κ|.

**Mechanism in one line.** For large |κ| every wave is outgoing sound. The sonic point lets nothing in (the background
exponent n* < 0 is a certified saddle), so the incoming wave w_− is only an O(1/(Kξ)) echo of the outgoing w_+; a
regular centre needs an incoming wave as strong as the outgoing one. Both cannot hold.

**Step 1 (sonic representation).** For Re κ ≥ 0 > d every sonic-analytic solution satisfies
w_−(x) = ∫_x^0 E_−(s,x)[b_−+w_+ + g_−](s) ds/ℓ_−(s), E_−(s,x) = exp(−∫_x^s (κ − b_−−)/ℓ_−): the homogeneous solution
behaves like |x|^{−(κ−d)/λ'}, which is unbounded, so its coefficient is zero.

**Step 2 (second-order Levinson slaving; Lemma 2.1 of PROOF-B4-v4, re-derived symbolically by a referee).** Put
ψ = σ̂' + (β_+ + β_−)σ̂, χ = ψ/T, φ = χ' + (β_+ + β_−)χ, q = e_−A α̃_+/(ℓ_− T), δ = −q/(κ(κ − A)), and

  p = σ̂/κ + χ/κ² + δ,   z_− = w_− − p w_+.

Then exactly z_−' = a_− z_− − F with a_− = (κ − b_−−)/ℓ_− and
F = [φ/κ² + (β_+ + β_−)δ + δ'] w_+ + R_− + p h_+, where R_− = [−e_−A α̃_− w_−/(κ − A) + E_−N N̄]/ℓ_− and
h_+ = τ̂ w_− + g_+/λ_+. The choices of σ̂, χ and δ cancel the O(1), O(1/κ) and Ā-α̃_+ parts of the w_+ source, which is
left as O(1/(Kξ)²) with no total-variation term. p is analytic at s = 0, and z_−(t) = ∫_t^0 E_−(s,t)F(s) ds.

**Step 3 (bulk bootstrap).** Unknowns U = sup|w_+|/e, Z = sup|z_−|/(ρ_z e), N = sup K|N̄|/(ν e) on [x_c, 0], with
x_c = ln(M/K), e(s) = |E_+(s, x_c)| and free positive weights ρ_z, ν. Row 1 is the + equation integrated from x_c, row 2 is
Step 2, row 3 integrates N̄ by parts using the + equation (which gains 1/K in the lapse loop). This gives a 3×3
nonnegative matrix 𝔄 with v ≤ 𝔄v + e_1. **Collatz test:** if 𝔄y < y for some y > 0, then ρ(𝔄) < 1 and
v ≤ (I − 𝔄)⁻¹e_1, which yields at x_c: |z_−| ≤ e_b|w_+|, K|N̄| ≤ n̄|w_+|, and w_+(x_c) ≠ 0. Integrals are bounded by upper
Riemann sums with exact exponential weights on cells; next to the sonic point (|s| ≤ s0 = 0.0036260) we use
ℓ_−(s) = |s|λ_−'(ξ) with λ_−' ∈ [1.0986, 1.0995] and b_−− < 0, so Γ(s, −s0) ≤ e^{s0·max(0, β_+)}(|s|/s0)^g with
g = 0.0896 > 0; integrability does not need Re κ > 0.

**Step 4 (centre zone, uniform in K ≥ K_a).** On 0 < v ≤ v_a = pM/K_a write the centre system as a flat part plus a
perturbation Δ = Δ_L + Δ_R with |Δ| ≤ Δ̂(t)/K and |Δ − Δ_L| ≤ Δ̂_R(t)/K², using exact facts on the certified centre
manifold (C0(0) = C00, C1(0) = 0, C1'(0) = C11, P = θ + O(θ²); asserted in ball arithmetic, 128 checks). The flat
kernel has a closed form in entire functions (Appendix B). A Volterra supersolution at K = K_a scales as (K_a/K) for every
K ≥ K_a, so one evaluation per θ-cell covers all K ≥ K_a:
ρ_m(x_c) ∈ Y^f_m(Zε) + (1/K)B1^(1)_m(θ) ± E2_m/K², with Z = √3pM independent of K.

**Step 5 (contradiction).** Every centre-regular solution is α(ρ_ω + μρ_N) or ρ_N. The bulk bound |N̄| ≤ (n̄/K)|w_+|
confines μ to a disc of explicit radius r around μ0 = −a/b_N (and excludes μ = ∞); on that disc
|w_− − c_b w_+| ≥ |g0 + μ0 g1| − r|g1|. If this exceeds e_b(|A_w+ + μ0 B_w+| + r|B_w+|), no eigenvalue exists in that
(K-cell × θ-cell). Conjugate symmetry reduces θ = arg κ to [0, π/2].

**Certificates.**

| band (Re κ ≥ 0) | M | instrument | K-cells | worst Collatz | worst lower/rhs | log |
|---|---|---|---|---|---|---|
| 40 ≤ \|κ\| ≤ 55 | 2.5 | fw_cert3 | 11 | 0.7566 | 1.010825 (K in [40, 41.2]) | runs6/fw_cert3_M2.5_Ka40_K55.log |
| 55 ≤ \|κ\| ≤ 2000 | 3.5 | fw_cert3 | 122 | 0.5682 | 1.314537 (K in [55, 56.65]) | runs6/fw_cert3_M3.5_Ka55_K2000.log |
| \|κ\| ≥ 2000 | 3.5 | fw_tail2 | uniform | 0.54072 | 14.77 | runs6/fw_tail2_KT2000.log |

Each band uses 160 θ-cells. M is a free matching parameter and may differ per band. The thin margin is at the bottom of
the first band (lower/rhs 1.010825, absolute margin 7.48e-4, at θ next to π/2); the first band is the reason R1 could be
40 rather than 55.

**Tail (|κ| ≥ 2000; Appendix C).** With u = 1/K ∈ (0, 1/2000], the outer region [x_T, 0] (x_T = −6.3303) is monotone in u
for weights ρ_z ∝ u, so its bounds are evaluated once at K = 2000; the deep region [x_c(K), x_T] is handled by an exact
θ-Laurent blow-up of every coefficient on the certified centre manifold and closed-form integrals in y = Kξ. Uniform
result for all K ≥ 2000: Collatz bound 0.54072, e_b ≤ 0.06146, n̄ ≤ 0.48433, and the contradiction test passes on all 160
θ-cells with worst lower/rhs 14.77. Controls (all pass): envelope against cell enclosures on the 3 deep cells (147
comparisons, 0 non-overlaps; the same check on the 3 cells just above x_T gives 26 non-overlaps, so it can fail);
monomial bookkeeping against cell integrands within [0.9789, 1.0738]; measured outer monotonicity at K = 4000 and 1e9.

**Lean.** Lean 4 with Mathlib kernel-checks the Collatz step of the Far Wall, with standard axioms only:
- **Rows.** The rational Collatz rows 𝔄_up·y < y, y > 0 (decimal upper bounds at or above the interval entries), for
  every one of the 133 bulk K-cells of [40, 2000] (11 at M = 2.5, 122 at M = 3.5; 532 theorems; ledger 16cb6d873a).
  The same holds for the uniform tail matrix at K ≥ 2000 (Collatz ratio 0.540738 in rationals; ledger 6a2a9e1d2f).
- **Logic.** The first sentence of Corollary 3.1, for real 3×3 matrices: if 0 ≤ 𝔄 ≤ M ≤ 𝔄_up entrywise, these rows hold, v ≥ 0 and
  v ≤ 𝔄v + e₁, then I − M is invertible, (I − M)⁻¹ ≥ 0, every eigenvalue of M has modulus < 1, and
  v ≤ (I − M)⁻¹e₁ (ledger 750fcfda8b).
- **Teeth control.** The same lemma with A·y ≤ y in place of A·y < y fails to check, as it must (ledger 3078ef4ec1).

Not in Lean: the interval-arithmetic computation of the matrix entries, the analytic estimates that produce the matrix,
the arb evaluation of (I − 𝔄)⁻¹e₁ and the bounds e_b, n̄ computed from it, the three consequences at x_c (e_b, n̄,
w_+(x_c) ≠ 0), and the centre and μ-disc tests. The Far
Wall as a whole must not be described as Lean-checked.

## 8. κ* and γ (Corollary 6)

With N = 1 and κ* real, `count_b4v3.py kstar` encloses the zero by bisection on the real axis:
κ* ∈ [2.810552549 ± 3.83e-10], γ = 1/κ* ∈ [0.3558019224 ± 6.97e-11] (runs7/kstar.log, 67 s).

**Comparison with print (POST HOC).** Hara, Koike and Adachi (1996) printed κ = 2.81055255 (β ≃ 0.35580192)
in their text, which they "believe to be exact to the last digit" (gr-qc/9607010, l.2492-2496), and 2.8105525488 in
their table (l.2545); KHA (1995) printed β ≃ 0.3558019; Maison
(1995) printed γ = 0.3558; Evans and Coleman (1994) found ≈ 0.36 in simulations. Our κ* enclosure contains 2.8105525488
(difference about 2e-10) and our γ contains 0.3558019 and 0.3558. The first enclosure is about 10 significant digits wide.
A non-counting probe (kstar11) showed that its floor of about 1e-10 is the width of the stored background tube, not the
eigenvalue method. The blind second implementation (below) works on the 35-digit B5 background and so has no such floor.
Its enclosure [2.8105525488271268040, 2.8105525488272177534] rounds to 2.8105525488 at HKA's eleven digits. It was
computed on the B5 solution. If the B5 and E5 solutions coincide (not claimed; Section 11), HKA's table value is correct
to its last digit and two more digits are certified (13 significant digits in all). The first engine alone, on E5,
agrees with HKA's 9-digit text value 2.81055255 (enclosure [2.8105525486170, 2.8105525493830]).

**What γ means here.** γ = 1/κ* is the exponent that the renormalisation-group/scaling argument of KHA and Maison
predicts for the black-hole mass near threshold, M_BH ∝ (p − p*)^γ, from a single growing mode. That argument, and HKA's
conditional nonlinear Theorem 1, are not proved here.

**Blind second enclosure (sealed test G2, 8a1c27c7fb with amendment b62c6aecc9).** A blind seat, whose brief quoted the
literature exponent ≈ 0.3558 (contrary to the bar's no-target clause; disclosed) and which was given no path to our
engine, derives its own perturbation equations on the B5 background and encloses a real growth rate;
the bar requires overlap with κ*, a known-answer control on exact Friedmann at k = 1/3, and a certified constant sign on
[0.40, 2.70] after gauge removal.

**Result: HOLD.**
- **The blind formulation.** Polar-areal slicing with central lapse α(t, 0) = 1; unknowns (δμ, δΩ, δV, δN/N) in the
  independent variable z = 1/N; interval Taylor models with Lohner QR; Frobenius solutions at the centre and the sonic
  point. The certified sonic exponent λ(κ) < 0 on every interval used rules out sonic resonance.
- **PG1.** κ* ∈ [2.8105522155761715, 2.8105525970458984] (width 3.8e-7 ≤ 1e-6), overlapping our enclosure. Refined
  in two runs that differ in Taylor order, step sizes, matching point and precision: κ* ∈ [2.8105525488272 ± 7.32e-14]
  and γ ∈ [0.35580192244308 ± 8.03e-15], identical in both runs and inside our enclosure. This 13-digit refinement is
  certified, beyond the sealed bar (PG1 required width ≤ 1e-6), and it is computed on the B5 solution.
- **PG2' (known answer).** On exact flat Friedmann the pipeline encloses the hand-derived rate 1 in [0.99999, 1.00001].
  - This rate is the time-translation mode, so the control is weaker than a physical rate would be.
  - In this formulation, Friedmann's physical polynomial modes fall exactly on sonic resonances and cannot be detected.
- **PG3.** After dividing out the gauge zero, G = F/(κ − 1) is certified negative on all of [0.40, 2.70].
- **Gauge.** In polar slicing the only residual coordinate freedom is t → t + const, giving one gauge zero at κ = 1
  exactly; its eigenfunction matches the closed form. The first engine's gauge zero κ_g = 0.3556992037 belongs to its own
  slicing. The amendment written before launch anticipated exactly this.
- **Controls.** The background, re-integrated from the B5 parameters, overlaps from both sides at the matching point. The
  unused momentum constraint has residuals containing 0.
- **Disclosed deviations.** The seat's instruments were pinned after its runs, not before. The seat knew the published
  value ≈ 0.3558 and located the root by a coarse scan of [0, 4].

## 9. Sealed tests and verification

Every test below was sealed (prediction and pass bar hashed) before its run and scored as written. Amendments were
sealed before the runs they govern; AMEND-LB-E-2 and AMEND-LB-B5C-1 came after non-counting smokes and are labelled POST
HOC with respect to those smokes.

| test | seal | what | verdict |
|---|---|---|---|
| E | 9ed31dfab0 (+ 921909428f, 56439ea867) | existence, attempt 1 | FAILED (warm-up did not certify) |
| E2 | 95f6de5cf2 | existence, attempt 2 | FAILED (warm-up ladder did not contract) |
| E3 | 8aec628fd2 | existence, attempt 3 | VOID (settled HOLD, corrected: series-cap truncation in the tail input) |
| E4 | 7928715efe | existence, attempt 4 | FAILED (k = 1/2 warm-up did not contract) |
| B5C | 71904980c4 (+ bee365c702) | blind B5 closure checks | HOLD |
| E5 | d0af75803f | existence, attempt 5 (with R8) | HOLD |
| (E5 duplicate) | 9592b015f0 | the same E5 file sealed twice; this copy was never used or settled | VOID (duplicate, never run) |
| M | b1b0ee1cde (+ b92e091b0c, 4971bef2bf) | count on R1 = 40 + κ* | HOLD |
| G2 | 8a1c27c7fb (+ b62c6aecc9) | blind second κ* | HOLD (instruments pinned post-run, disclosed) |

Ledger: 9 seals (4 HOLD, 3 FAILED, 2 VOID); the duplicate 9592b015f0 is VOID (duplicate, never run). Four sealed tests
HOLD: E5, B5C, M and G2.

**Disclosures.**
- E5 disclosed the E3 and E4 numbers, the R8 smokes and the B5 proof. The bands came from published values and exact
  algebra.
- M disclosed non-counting smokes (including R1 = 3 dry runs with N = 1 and a smoke κ*); its prediction is the
  literature's claim of one growing mode. The wider half-plane wording (AMEND-LB-M-2) was sealed before any count data
  and made conditional on the tail referee; its condition was met (fixes rerun).
- B5 had no sealed bar; its agreement with E5 is POST HOC. The B5 brief quoted Evans and Coleman's m/r; the seat says no
  target was used as an input.
- Non-counting probes are not results: the kstar11 tightness probe (dropped, no seal) and all smoke logs.
- Early Far Wall logs (fw_cert, fw_cert2) were void because the series cap was never raised; they were superseded by the
  fw_cert3 runs. The first tail run's "positive control" was withdrawn as mislabelled and too loose, and replaced by the
  controls in Section 7.

**Internal referee seats.** "Seats" are separate internal working sessions of the lab with their own files. Referee and
blind seats are internal checks, not external peer review. A blind seat did not see our code or the values it was asked
to reproduce, except as disclosed. Seven hostile referee reports cover this work:

| seat | object | findings | what was done |
|---|---|---|---|
| REFEREE-E | E3 existence run | equations, tails, integrator CORRECT; Krawczyk existence-only, with four sub-1e-25 enclosure gaps; PE2 inward GAP; two controls could not fail; wording "first interval-arithmetic certificate", never "first existence proof" | v3.1 repairs, content nulls, controls at sealed parameters; E3 later VOID for a separate cap defect |
| REFEREE-R8 | repair R8 before E5 | R8 CORRECT; control file weak (injected test could not fail, no derivative check); a reported k = 1/2 "resonance" was an instrument artefact (90 digits) | control fixed and pinned before sealing E5 |
| REFEREE-B5 | blind b5/ | equations CORRECT; KHA (18) line 3 is a misprint; no float leaks; three unchecked pieces (outer sonic stretch, det C, z-monotonicity on the centre step); identity with Evans-Coleman's printed numbers only to 0.25% / 4% | closure test B5C (HOLD); analytic arguments written (Appendix D) |
| REFEREE-E5CAP | series-cap audit of E5 | E5 SAFE: every cap-sensitive object bit-identical at cap 10 and 1024; positive control fires | E5 HOLD stands |
| REFEREE-B4FW | Far Wall v4 and the frame change | mathematics CORRECT WITH FIXES; frame change CORRECT; certificate runs WRONG as run (series cap never raised, order-11/22 remainders dropped); switch-sliver GAP | fw_cert3 (cap 32, strict coefficients, radii rounded up, Lemma 4.1 asserted, sliver bootstrap), reruns, AMEND1 |
| REFEREE-B4FW2 | repairs + count wrapper | all repairs CORRECT; arc chords dipped 4.9e-15 below 40 in floating point; parity printed but not enforced; two preconditions not certified in code | vertex inflation with exact check, parity raises, preconditions certified at run start, fresh ledger with header |
| REFEREE-B4TAIL | tail \|κ\| ≥ 2000 | claim holds; positive control WRONG as a control; hygiene fixes | fw_tail2 with new controls, strict assertions, rerun; AMEND2 updated |

No report found an error in a certified value that survived into the results above. This release text had one internal
hostile referee round (seats/REFEREE-RELEASE.md, READY WITH FIXES, fixes applied).

**Clean-room replay.** From a fresh copy of the repository, using only the documented commands (`tools/reproduce.md`):
all 108 pinned files match their ledger sha256; the full tier regenerated all nine data files byte for byte and reproduced
every certified run of E5, B5C, M (including the radius-40 count and both Far Wall bands) and G2, with every output equal
to the lab's after removing run times only (`results/CLEANROOM-REPORT.md`; 3.5 hours on one core).

**Verification table (what checks what).**

| claim | first witness | second witness | literature (POST HOC) |
|---|---|---|---|
| existence | E5 (Taylor models + Krawczyk) | B5 (Chebyshev + Poincaré-Miranda), blind | HKA 1996 shooting assertion; numerics since 1994 |
| saddle sonic point | E5 eigenvalues | B5 eigenvalue ratio −0.0901059705... | Carr et al.; Harada (saddle for k ≲ 0.41) |
| m/r(t = 0) | E5 [0.0595 ± 8.24e-5] | B5 0.059451404935530220048... | Evans-Coleman 0.0596 (Section 10) |
| one growing mode | M count + Far Wall | G2: no other real zero on [0.40, 2.70] after gauge removal (real-axis part only) | HKA numerical uniqueness (0 ≤ Re κ ≤ 15, \|Im κ\| ≤ 14) |
| κ*, γ | M kstar | G2: [2.8105525488272 ± 7.32e-14], [0.35580192244308 ± 8.03e-15] | HKA 2.8105525488; KHA 0.3558019; Maison 0.3558 |
| κ_g | C1c enclosure | closed form 2 − A0 + (2/3)ω0 | HKA print 0.35699 (Section 10) |

## 10. Fairness notes and corrections

1. **KHA (1995) eq. (18), third line.** As printed (read from the page image) it is not satisfied by our derived system,
   by the blind seat's independently derived system (the two agree to 5e-38), or by the exact flat Friedmann solution,
   which satisfies lines 1, 2 and 4 exactly (REFEREE-B5). We read it as a misprint. Three single-term repairs tested by
   the referee all fail, so we cannot say which term is wrong.
2. **HKA's gauge mode 0.35699 (observation).** HKA define the gauge mode in their gauge N̄_p(0) = 0 as −dN̄_ss/dx(0) and
   print 0.35699 (gr-qc/9607010 l.2517; also gr-qc/9503007). With the same definition we certify
   κ_g ∈ [0.3556992037 ± 7.41e-11], which is the closed form 2 − A0 + (2/3)ω0. The printed value looks like 0.356992 with
   one digit dropped. We label this an observation; we did not contact the authors.
3. **Evans and Coleman (1994), m/r and Ω at t = 0 (observation, cause UNTESTED).** They quote m/r → 0.0596 and
   Ω → 9.56e-3. Both of our implementations give m/r = 0.0594514... (0.25% lower) and Ω = 0.0091740... (about 4% lower),
   with the same definitions (t = 0, Ω = 4πr²ρ; checked by the B5 referee). Their a ≈ 1.07 agrees with ours (1.0653...) at
   their precision. The most likely cause is 1994 numerical accuracy; we have not tested this.
4. **HKA (1996), "mathematically rigorously".** HKA wrote that the self-similar solution is seen "both numerically and
   mathematically rigorously" to exist (l.2231-2234), while calling their analysis not mathematically rigorous elsewhere
   (l.2749). Their shooting argument is an existence assertion without enclosures, not a checkable proof (REFEREE-E). So
   we claim a first interval-arithmetic (computer-assisted) certificate, not a first existence proof. HKA also proposed
   the route we followed: close analysis of the linearised spectrum "may lead to a computer-assisted proof" (l.2692-2697).
5. **Maison's "attractive node".** Maison's sentence about a node concerns k near k_max ≈ 0.888, not k = 1/3; at k = 1/3
   the sonic point is a saddle (Harada; Carr et al.; our certificates). No conflict.
6. **Our own corrections (none of this was published).**
   - Two lab notes said the Newtonian Hunter A solution had been proved to exist; Sandine (2025) says it is open.
   - A lab roadmap said "third critical solution ever proved"; the Bartnik-McKinnon (type I) solution was proved in 1991.
     The correct count is "third self-similar (type II)".
   - A lab menu said "first existence proof"; retired (item 4).
   - E3 was settled HOLD and later corrected to VOID (Section 9).
   - A relay said Maison's "Note Added" disputes the gauge mode; the seat that read it in full found no such statement.
   - Two Far Wall logs were void as certificates (series cap) and were rerun; the first tail control was withdrawn.

## 11. Limitations (what we do not claim)

- **Local uniqueness** of the solution is not claimed by either implementation, and we do not claim that the E5 and B5
  solutions are the same solution. Other analytic, centre-regular solutions with one zero of V near the node segment of the
  sonic line are not excluded (B5, disclosed GAP).
- **Linear, spherical only.** No nonlinear stability, no proof of the critical-collapse mass law, no statement about
  non-spherical perturbations (numerically studied by Gundlach, 1999). The count is in the class analytic through the sonic
  point and regular at the centre, modulo the gauge mode; non-analytic (kink) perturbations are addressed only by Harada's
  criterion (Remark 7).
- **Modes, not spectrum.** Absence of continuous or residual spectrum for a linearised evolution operator on a function
  space is not proved.
- **Trusted layer.** Ball arithmetic (python-flint arb/acb) with outward rounding is trusted. The Far Wall's interval
  certificates are the proof; Lean covers only the Collatz step stated in Section 7: the rows on every bulk cell and the tail, and the implication to the bootstrap bound. The matrix entries are ball-arithmetic enclosures and are not Lean-checked.
- **Classical facts used without opening a textbook:** Frobenius theory at regular singular points, the Briot-Bouquet
  theorem, the strong-unstable manifold theorem for analytic vector fields, the argument principle, Brouwer and
  Poincaré-Miranda.
- **One-sided controls.** The far t = 0 leg of B5 has no exact-solution positive control. The E5 outward leg uses the
  rigorous but looser v3.1 head. The count has no exact-solution positive control: the k = 1/2 Friedmann control C3
  failed as an instrument during the build and is not part of BAR-LB-M.
- **Agreement with the literature is POST HOC**, and so is the B5-E5 agreement.
- **Digits.** The first engine certifies κ* to about 10 significant digits. The 13-digit enclosure comes from the blind
  second implementation alone, run on the B5 solution (which we do not claim equals E5); the first engine confirms only
  its first 10 digits. We stopped there and attempted no
  identification of κ*.
- **Review status.** One internal hostile referee round on the release text (seats/REFEREE-RELEASE.md, READY WITH
  FIXES, fixes applied); prior-art refresh to 2026-11-16.

## 12. Related work and credit

- **Evans and Coleman (1994)** observed critical phenomena in the spherical collapse of a radiation fluid, found the
  self-similar solution and a critical exponent ≈ 0.36 in simulations; they credit the precise self-similar solution to
  Evans.
- **Koike, Hara and Adachi (1995)** gave the renormalisation-group picture β = (Re κ)⁻¹ and found the relevant mode
  numerically, β ≃ 0.3558019. **Hara, Koike and Adachi (1996)** gave the full numerical analysis, κ = 2.8105525488, wrote
  that their results "strongly suggest that the relevant mode is unique", and proposed a computer-assisted proof. The route
  we followed is theirs.
- **Maison (1995)** derived the linear perturbations independently and gave γ(k), with γ = 0.3558 at k = 1/3.
- **Harada (2001)** gave the analytic kink-mode criterion; the Evans-Coleman solution is stable for this mode for
  0 < k ≲ 0.89.
- **Bizoń and Wasserman (2002)** proved existence of the first (CSS) critical solution, by shooting.
  **Reiterer and Trubowitz (2012)** proved, with computer assistance, that Choptuik's DSS critical spacetime exists. With
  these, Theorems 1-2 make the Evans-Coleman solution the third self-similar (type II) critical solution with an existence
  proof in the list of the Living Review (§VII C), and the first for a perfect fluid. A type I critical solution
  (Bartnik-McKinnon) was proved to exist by Smoller, Wasserman, Yau and McLeod.
- **Gonçalves (2026)** proved with computer assistance that Choptuik's solution has exactly one unstable spherical mode,
  modulo gauge (λ ∈ [2.674040, 2.674115]), with a rank-two deflation of the gauge roots; this is the DSS analogue of our
  count. Ours is the first certified mode count for a continuously self-similar critical solution of gravitational collapse, and
  the first for a fluid.
- **Gundlach (1999)** studied non-spherical perturbations numerically: at k = 1/3 one spherical growing mode and no
  non-spherical growing modes.
- **Guo, Hadžić and Jang (2021)** proved existence of the relativistic Larson-Penston solution, a different (attractor,
  not critical) self-similar Einstein-Euler solution. **Guo, Hadžić, Jang and Schrecker (2025)** excluded growing modes of
  Newtonian Larson-Penston with an energy method at the frequency extremes and interval arithmetic in between, a frequency split close to ours
  (interval arithmetic in the middle; analytic bounds at the extremes, energy-based in theirs, a Levinson-slaving
  bootstrap in ours).
- **Goldstein and Jejjala (2026)** used Taylor models, a stable-manifold tail and Poincaré-Miranda for an Einstein-Weyl
  black hole. **Sandine (2023, 2025)** constructed Hunter-type implosion profiles for Euler-Poisson and states that the first
  Hunter solution remains open. **Carr, Coley, Goliath, Nilsson and Uggla (1999)** located the saddle region (k ≲ 0.41).
- The Living Review of Gundlach, Hilditch and Martín-García (2025) states that existence and stability of the relativistic
  Hunter solutions had not been proved.

## References

Labels: OPENED = the lab read the text (copy in the lab's source folder; locators are line numbers there); ABSTRACT =
abstract or record only. arXiv identifiers are given; journal data only where the opened text shows it.

- [BW] P. Bizoń, A. Wasserman, arXiv:gr-qc/0201046, l.14, l.35-40. OPENED.
- [Carr] B. J. Carr, A. A. Coley, M. Goliath, U. S. Nilsson, C. Uggla, arXiv:gr-qc/9902070, l.1083-1087. OPENED.
- [EC] C. R. Evans, J. S. Coleman, arXiv:gr-qc/9402041, l.13, l.20, l.36, l.134-141, l.303. OPENED.
- [GHJ] Y. Guo, M. Hadžić, J. Jang, arXiv:2112.10826, Thm 1.1 (l.215). OPENED.
- [GHJS] Y. Guo, M. Hadžić, J. Jang, M. Schrecker, arXiv:2509.12435, l.524-540. OPENED.
- [GHM] C. Gundlach, D. Hilditch, J. M. Martín-García, Living Review, arXiv:2507.07636 v1, §VII B l.2279-2282, §VII C
  l.2302-2306. OPENED.
- [GJ] K. Goldstein, V. Jejjala, arXiv:2609.21797 v2, abstract and pp. 1-2. OPENED.
- [Gon] E. A. Gonçalves, Zenodo record 10.5281/zenodo.23270506 (v2, 2026), description. ABSTRACT.
- [Gun] Gundlach, "Critical gravitational collapse of a perfect fluid: nonspherical perturbations", arXiv:gr-qc/9906124, abstract. ABSTRACT.
- [Har] T. Harada, arXiv:gr-qc/0109042 v2, abstract (l.15-24), l.470-478, l.785-795, l.900-908. OPENED.
- [HKA] T. Hara, T. Koike, S. Adachi, arXiv:gr-qc/9607010, l.54, l.108-111, l.2231-2234, l.2492-2496, l.2517, l.2545,
  l.2692-2697, l.2721-2725, l.2749. OPENED.
- [KHA] T. Koike, T. Hara, S. Adachi, arXiv:gr-qc/9503007, abstract l.10-17, eq. (16), eq. (18), l.210-218. OPENED.
- [Mai] D. Maison, arXiv:gr-qc/9504008, eq. (4), eq. (8), l.128-133, l.176. OPENED.
- [RT] M. Reiterer, E. Trubowitz, arXiv:1203.3766, abstract. OPENED.
- [San23] E. Sandine, arXiv:2310.18870, l.299-301. OPENED.
- [San25] E. Sandine, arXiv:2508.05817, l.264-268, l.359. OPENED.
- [SWYM] Smoller, Wasserman, Yau, McLeod, "Smooth static solutions of the Einstein-Yang/Mills equation",
  arXiv:math/9210226 (1992), abstract page. ABSTRACT.

---

## Appendix A. Far Wall bulk bootstrap (PROOF-B4-v4 sec. 3, AMEND1)

**Rows.** Row 1 uses the + equation from base x_c with h_+ = τ̂(p w_+ + z_−) + g_+/λ_+ and the propagator bound
|E_+(s,t)| e(s)/e(t) = 1. Row 2 is the sonic representation of z_−, with the N̄ part of R_− integrated by parts on the E_−
phase using N̄' = AĀ + H3 w (κ-free). Row 3 is K N̄(t) = −K∫_t^0 (f_w w_+ + f_z z_− + f_N N̄) with
f_w = H3(1 + p) − A(α̃_+ + α̃_− p)/(κ − A), f_z = H3 − Aα̃_−/(κ − A), f_N = −Ac_N/(κ − A); the term ∫ f_w w_+ is integrated by
parts with w_+ = (b_++w_+ + b_+−w_− + g_+ − λ_+w_+')/κ, so the factor K/|κ| = 1. All entries of 𝔄 are monotone in absolute
coefficient bounds.

**Discretisation.** For t in cell i, ∫_t^0 W(s,t)f(s)ds ≤ Ω_i(f_i h_i + G_i), G_{i+1} = Ω_i f_i h_i + E_i G_i, with
E_i = exp(up(h_i β_i)), Ω_i = exp(max(0, up(h_i β_i))) and f_i the cell supremum. For a K-cell [K_a, K_b] the program
integrates over every cell meeting [x_c(K_b), 0]; this is valid because the integrands are nonnegative and the
normalisation e(s) only rescales all bounds by a common factor (AMEND1 sec. 3).

**Enclosures.** Coefficient functions are degree-5 Taylor models on each bulk cell, built from a thin start hull
(width about 3e-15) with an order-11 Lagrange remainder on the certified a-priori box; steps are cut into 2^j sub-steps until
coef_11·(h/n)^11 ≤ 1e-13 (valid for any n). The series cap is 32 before anything is built, and any request for an unknown
coefficient raises. The background is re-run with the sealed E5 instrument on a box containing the E5 K box, with radii
rounded up. The τ/x switch sliver is closed by an explicit bootstrap: x_end ∈ −0.050000000 ± 7e-14, sliver length
ℓ ≤ 5.41e-14, sup|f_bg| = 0.727 on the last sub-cell box plus 1e-6, so every state on the sliver is within 3.93e-14 of
the landing state.

**Weights.** Any fixed positive weights are admissible; floating-point power iterations choose them and every candidate
is a complete rigorous certificate, so keeping the best is selection among valid proofs.

## Appendix B. The flat centre kernel (PROOF-B4-v4 sec. 4.2-4.3)

In ζ = tε (ε = e^{iθ}), with ρ = t_y/t_x ∈ [0, 1] and d = (1 − ρ)ζ_x, the 4×4 kernel K_f(x, y) = Φ(ζ_x)Φ(ζ_y)⁻¹ of
Y' = (C00 + κvC11)Y (rows and columns Ā, N̄, ω̄, V) is

  K = [[ρ, 0, 0, 0], [1−ρ, 1, 0, 0], [−3(1−ρ)sh1, 0, (1−ρ)sh1 + ρ ch, −(4/√3)ρ d sh1],
       [(3√3/4)ζ(ρ(1−ρ)²c2 + (1−ρ)³g3), 0, −(√3/4)ζ((1−ρ)³g3 + ρ(1−ρ)sh1), ρ²ch + ρ(1−ρ)³ζ²g3]],

with sh1 = sinh d/d, ch = cosh d, g3 = (d cosh d − sinh d)/d³, c2 = (cosh d − 1)/d², all entire. The flat regular
solutions are Y_ω^f = (0, 0, sh1(ζ), −(√3/4)ζ g3(ζ)) and Y_N^f = e_N̄. The majorant uses |sinh d|, |cosh d| ≤ cosh|Re d| and
|sh1(d)| ≤ sh1(|d|), |g3(d)| ≤ g3(|d|), |c2(d)| ≤ c2(|d|) (positive series), valid at θ = π/2. A referee checked the
kernel against a direct ODE solve (errors 4.8e-13 to 5.6e-11 at four points including ζ = 5i).

**Volterra supersolution and K-scaling.** With e = ρ_m − Y^f_m, |e| ≤ 𝒦(|Y^f| + |e|) for a positive Volterra operator
with kernel |K_f|Δ̂/K. If ε_a ≥ (1/K_a)𝒦̂(|Y^f| + ε_a) pointwise (checked on 150 log-cells, with the tail t ≤ 1e-3 as one
cell), then |e| ≤ ε_a at K = K_a, and (K_a/K)ε_a is again a supersolution for K ≥ K_a.

## Appendix C. The tail |κ| ≥ 2000 (AMEND2-TAIL)

**Outer region, monotone in u.** Every integrand is built by + and × from nonnegative cell bounds and the factors
iK = u, iK2 = u², iKA = u/(1 − A_max u), ρ_z ∝ u and Kρ_z = constant, all nondecreasing in u; every row-2 term carries at
least one u-factor. Hence the outer bounds for every K ≥ 2000 are at most their values at K = 2000.

**Deep region.** On the centre manifold x = log θ + χ(θ) + xc, write P = θP̃(θ) with P̃(0) = 1; then every bulk
coefficient is an exact Laurent quantity θ^m s(θ), giving envelopes on θ ∈ [0, θ_env], θ_env = 0.001455. θ_env covers the
deep region: the program asserts θ_env > F(θ_env) strictly (margin 1.0100000), so x(θ_env) > x_T, and by continuity and
θ → 0 as x → −∞, θ(x) < θ_env for all x ≤ x_T. The envelopes reproduce the flat limits (ξσ̂ = −0.357 = −c'/2,
β_+ = −1.000 ± 1.3e-3). Deep quantities are nonnegative sums Σ c ξ^a y^{−b} with y = Kξ. **Closed-form integral:** with
cc = a − b + e, ∫_t^{x_T}(y_s/y_t)^e ξ_s^a y_s^{−b} ds ≤ ξ_t^a y_t^{−b}/(−cc) if cc < 0, and
≤ ξ_T^{cc} ξ_t^{b−e} y_t^{−b}/cc if cc > 0; when |cc| < 1/20 the exponent e is raised to make cc = 1/20, which is valid since
(y_s/y_t)^e is nondecreasing in e. Junction terms carry an explicit ξ_T^{−α} step (written out after the referee).

## Appendix D. The two analytic arguments of the blind proof (b5/ARGUMENTS.md)

**D.1 Centre.** With Ω = z²σ, V = zU: g(0, σ, U) = −3U − 2/(1+k), independent of σ, with unique root U0 = −2/(3(1+k)) and
slope −3. A bounded solution must have U(0) = U0 (otherwise U ~ c log z). Writing g = −3(U − U0) + h gives
U = U0 + K3[h], K3 f(z) = z⁻³∫₀^z s²f(s)ds, of sup-norm 1/3 and causal, so the certified contraction restricts to every
[0, ε]. The Briot-Bouquet theorem (no positive-integer resonance, since −3 ∉ {1, 2, ...}) gives an analytic regular
solution, which lies in the contraction ball on small [0, ε] and therefore is the certified fixed point.

**D.2 Sonic crossing.** At the sonic equilibrium the spectrum is {0, λ₊, λ₋} with λ₊ > 0 > λ₋ certified and 0 belonging to
the line of equilibria. The spectral gap makes the strong-unstable manifold unique, one-dimensional and analytic; its
parameterisation coefficients solve (J − nλ₊)p_n = known terms, solvable since nλ₊ ∉ {0, λ₊, λ₋} for n ≥ 2. Harada states
the same uniqueness for this sonic-point problem (l.475). The certified curve satisfies F(P(θ)) = λ₊θP'(θ), so it is
that manifold. Along it dx/dθ = S(P(θ))Q(P(θ))/(λ₊θ) with S/θ certified nonzero on both sides (prove.py and C1) and
Q ≠ 0, so the solution is real-analytic in x through the sonic point, with no second sonic crossing on the validated disk.
The branch is forced: on the inner side dx/dτ < 0, so the solution reaches the equilibrium as τ → −∞, which is only
possible along λ₊.
