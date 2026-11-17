# B5: the two analytic arguments the referee asked to see written out (REFEREE-B5 sections 3 and 4)

## A. The centre: the K3 fixed point is the regular (analytic) centre solution

**Setting.** Write Ω = z²σ and V = zU. The centre-regular system is (system2.json, systems.centre_system):

    dσ/dz = Fσ(z, σ, U),        z dU/dz = g(z, σ, U),

with Fσ and g rational and analytic near z = 0. Their denominators are nonzero at z = 0: they reduce to k·(-1) and
k(k+1)·(-1).

From system2.json, g(0, σ, U) = -3U - 2/(1+k). This does not depend on σ. It has the unique root U0 = -2/(3(1+k)) and the
slope ∂g/∂U(0) = -3. Also Fσ = O(z) (dsig carries an explicit factor z).

**1. Every solution bounded at 0 has U(0) = U0.** If U is bounded and continuous on (0, ε] with a limit U(0), then
z dU/dz → g(0, σ(0), U(0)). A nonzero limit c forces U ~ c log z, which is unbounded. So g(0, σ(0), U(0)) = 0, and
therefore U(0) = U0.

**2. Integral form.** Write g = -3(U - U0) + h, where h = g + 3(U - U0) vanishes at (z = 0, U = U0) together with its
U-derivative. Then (z³(U - U0))' = z²h. For solutions bounded at 0, this is equivalent to

    U(z) = U0 + K3[h](z),      K3 f(z) = z⁻³ ∫₀^z s² f(s) ds.

The integration constant is fixed by boundedness: the homogeneous solutions C z⁻³ are unbounded. Also
σ(z) = σ0 + ∫₀^z Fσ. This is exactly the operator in rigode.py (ops 'int', 'K3').

**3. Causality.** K3 is causal: K3 f on [0, ε] uses only f on [0, ε]. Its sup-norm bound is 1/3, since
|K3 f(z)| ≤ z⁻³ ∫₀^z s² ds · sup|f| = sup|f|/3. So the contraction proved on [0, z1] (rigode.validate, Y + L r < r with
operator norm 1/3) restricts to every [0, ε], ε ≤ z1. On each [0, ε] the fixed point is unique in the sup-ball of
radius r around the numerical guess.

**4. The analytic solution exists and lies in that ball.**
- The system z dU/dz = g, dσ/dz = Fσ is of Briot-Bouquet type at z = 0. The linear part in U has eigenvalue -3, and there
  is no positive-integer resonance (a resonance would need -3 ∈ {1, 2, ...}).
- So for each σ0 there is a unique formal power-series solution with U(0) = U0, σ(0) = σ0. It converges, by the
  Briot-Bouquet theorem (classical; e.g. Hille, *Ordinary Differential Equations in the Complex Domain*, the chapter on
  Briot-Bouquet equations).
- This analytic solution is bounded, so by step 2 it solves the integral equation.
- On a small enough [0, ε] it lies within r of the guess, because both are continuous with the same value at z = 0.
- Hence it equals the certified fixed point on [0, ε], and by uniqueness of the regular initial-value problem for z > 0
  on the whole step.

**5. Physical regularity.** With μ = 2Ω(1 + kV² + (1+k)U)/(1-V²) = O(z²) and V = O(z): the centre has A → 1 and V → 0
(KHA condition (ii)). Ω = z²σ gives a finite central density. U is even in z (the right-hand sides are even in z for
even σ, U), so V is odd in r. The center is regular.

**6. Monotone z on the centre step.** C3 of BAR-LB-B5C certifies this: P/z < 0 on [0, 0.05], so z_x = zP/Q > 0 there,
with Q < 0.

## B. The sonic crossing: the validated curve is the analytic strong-unstable manifold

**Setting.** F = (zPS, ΩnO, (V²-1)nV/(k+1)) is the polynomial desingularised field, with dx/dτ = SQ.
- The sonic line e(z_s) consists of equilibria. sonic_sym.py proves S(e) = nO(e) = nV(e) = 0 and J t = 0 exactly, where
  t = de/dz_s.
- At the certified z_s, J = DF(e) has the spectrum {0, λ₊, λ₋}, with λ₊ = 0.16429349223752482… > 0 and
  λ₋ = -0.01480382456620133… < 0. Both are certified (result_ec.json), and λ₊λ₋ = c2 < 0 is certified. The eigenvalue 0
  belongs to the line of equilibria.

**1. Spectral gap.** λ₊ is simple and strictly larger than the rest of the spectrum {0, λ₋}.
- The strong-unstable (fast) manifold W^{uu}(e) is the set of points whose backward τ-orbit converges to e at rate
  O(e^{λ₊τ}). Faster than every rate in the centre-stable spectrum is what matters here.
- This manifold exists, is unique, is one-dimensional, is tangent to the λ₊ eigenvector, and is as smooth as F. For
  analytic F it is analytic.
- This is the strong (pseudo-)unstable manifold theorem for an equilibrium with a spectral gap. In the analytic case it
  is Poincaré's theorem on analytic invariant manifolds: the parameterization P(θ) with F(P(θ)) = λ₊θP'(θ) has
  coefficients fixed by (J - nλ₊)p_n = known terms. These are solvable because nλ₊ ∉ {0, λ₊, λ₋} for n ≥ 2 (nλ₊ > λ₊ > 0
  > λ₋). The Poincaré majorant argument gives convergence.
- Harada gr-qc/0109042 (l.475, opened by the referee) states the same for this sonic-point problem: "the analytic
  solution with each direction is unique for given ys".

**2. The certified curve is that manifold.**
- sonic.validate_tail proves that the curve P = H + ρ, with |ρ(θ)| ≤ C|θ|², satisfies F(P(θ)) = λ₊θP'(θ) on
  |θ| ≤ 0.05. It is the fixed point of the weighted Volterra operator, with K_i bounds 1/(2λ₊ - κ_i) > 0 asserted at
  sonic.py:207.
- By this conjugacy the τ-orbit through P(θ0) is θ(τ) = θ0 e^{λ₊τ}, so it converges to e as τ → -∞ at rate e^{λ₊τ},
  tangent to v₊.
- By uniqueness in step 1 the certified curve is (a piece of) W^{uu}(e), hence analytic in θ. The formal series of the
  analytic manifold also solves the same fixed-point equation, and the fixed point is unique in the ball.

**3. Analytic in x through the sonic point.** On the curve, dx/dθ = (dx/dτ)/(dτ/dθ) = S(P(θ)) Q(P(θ))/(λ₊θ).
- S(P(θ))/θ is analytic, since S(e) = 0 exactly.
- It is certified nonzero:
  - on [θ*, 0] in prove.py: S/θ ∈ [-0.25 ± 7e-3];
  - on [0, θ_o] by C1 of BAR-LB-B5C.
- Q = -z(1-V²)(1-μ) ≠ 0 there (μ < 1 and V² < 1 are certified).
- So x(θ) is analytic and strictly monotone. Its inverse θ(x) is analytic, and (Ω, V, z) are real-analytic functions of x
  across the sonic point x_s.
- Both halves (θ < 0, the inner side with S > 0; θ > 0, the outer side with S < 0) belong to the same analytic curve.
  There is no second sonic crossing on the validated disk, because S/θ ≠ 0 there.

**4. Which branch.** The inner solution comes from the centre with x increasing. On the inner side dx/dτ = SQ < 0, so it
reaches e as τ → -∞. That is only possible along the unstable eigen-direction λ₊; the λ₋ curve is reached as τ → +∞.
The saddle classification (c2 < 0) is certified, so the choice of branch is forced, not tuned.
