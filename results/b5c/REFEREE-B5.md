# REFEREE-B5: hostile referee on b5/ (blind second implementation, existence of the k = 1/3 CSS solution)

Seat: referee. I did not open engine/cert/driver_v32.py or PROPOSAL-V32.md. I did not re-run the proof. My checks are
in the <SCRATCH>: kha18c_fit.py, kha18_frw.py, ec_finiteN.py, recompute.py (all memcap, under 1 min each).

## In plain words
The proof is sound in its design and, as far as I can see, in its code. The equations are right. Two independent
derivations agree to 1e-38, and KHA's printed eq. (18) line 3 is a misprint: it fails even on the exact Friedmann
solution, which satisfies lines 1, 2 and 4. I found no float leak and no truncation trap. I did find three small
unchecked pieces: a short stretch of the outer sonic curve, the invertibility of the preconditioner, and z-monotonicity
on the centre step. Each is a few lines to fix and none should move a digit. The claim must stay "a solution", not
"the unique Evans-Coleman solution".

## Verdicts

**1. Equations: CORRECT (fix one reproducibility hole).**
- system.py:24-59 derives the constraint and asserts it equals Maison (8). It also asserts z_x/z equals Maison (4)
  and checks constraint propagation.
- At 4 random points, the B5 system and engine/derive.py (via eqs_gen.py) agree on ω_x, V_x, A_x and N_x to ≤ 5e-38.
- Both satisfy KHA (18) line 4 to about 1e-37. Both violate line 3 *identically*, with residuals 5.0814, 4.17666,
  2.97363 and 9.87289.
- The exact flat-FRW radiation solution (V = -Z/2, A = 1/(1-Z²/4), ω = (3/2)(Z²/4)A, N from engine/frw.py):
  - satisfies lines 1, 2 and 4 exactly;
  - leaves line 3 with residual (Z²+4)(3Z⁴-16Z²+48)/(12Z(Z-2)²(Z+2)).
- I read the page image (b5/img/kha_p3.png) exactly as system.py:70-71 transcribes it.
- Verdict: **KHA misprint, not a seat error.** The location is not identifiable from self-similar data, because A_x and
  N_x are algebraic. I tested three single-term repairs and all three fail.
- Fix: RESULT.md:34 cites frw_check.py, but that file does not exist in b5/. Restore it or remove the citation.

**2. Validated ODE enclosures: CORRECT WITH FIXES.**
- The criterion Y + opn·L·r < r (rigode.py:145-176) is a valid self-map plus contraction on a convex sup-ball.
- The remainder is enclosed:
  - truncation mass goes into err (cheb.py:56-63);
  - the product err is cheb.py:100;
  - the integral err is cheb.py:138;
  - the quotient defect is |R|/min|den| (rigode.py:132-137).
- K3 has norm 1/3 (rigode.py:143,145).
- No float leaks on the proof path:
  - floats appear only in guess, tolerances, step control and the preconditioner;
  - a rational guard is at rigode.py:19-22 and sonic.py:25-26;
  - there are no float literals in system2.json or sonic_sym.json.
- Fixes:
  - (i) The cache key (prove.py:189) has no code hash, and the final FRW control reused a cache from an earlier
    attempt (run_frw.out: "centre cache found").
  - (ii) Reload containment is tested only for G_mid (prove.py:217). G_lo, G_hi and G_box feed Miranda.

**3. Centre: CORRECT WITH FIXES (a write-up gap, not a math gap).**
- From system2.json dU, z·dU/dz → -3U - 2/(1+k) at z = 0. That expression does not depend on σ, has a unique root
  U0 and has slope -3.
- K3 is causal, so the fixed point restricted to [0,ε] is the unique fixed point in the restricted ball. The analytic
  (Briot-Bouquet-type) regular solution lies in that ball, so they coincide.
- Write these lines into RESULT.
- Also: P < 0 (z monotone) is not checked on [0, 0.05] (prove.py:141-152). It is analytically -z(1+O(z²)), but it
  should be checked.

**4. Sonic point: CORRECT WITH FIXES.**
- The weighted Volterra operator is valid:
  - the K_i bound (sonic.py:7-8,230) needs 2λ+ > κ_i, which is asserted (sonic.py:207);
  - d₀ and d₁ are exact zeros by the identities in sonic_sym.log;
  - DEG = 12N+4 ≥ 9N covers the degree of F.
- The strong-unstable argument is correct. λ+ > 0 and λ- < 0 are certified, and the line of equilibria contributes only
  the eigenvalue 0. The spectral gap (0, λ+) therefore makes the fast manifold unique and analytic, and nλ+ ≠ 0, λ-
  trivially.
- Harada gr-qc/0109042 (OPENED, l.475) supports this: "the analytic solution with each direction is unique for given
  ys". Name the theorem in the write-up.
- Fixes:
  - (i) The checks for S/θ ≠ 0, P < 0, μ < 1 and V > 0 cover only θ ∈ [θ*, 0] (prove.py:316-342). The outer piece
    θ ∈ [0, θ_o ≈ 0.0311] (run_ec.out) is unchecked, and "no second sonic point" needs it.
  - (ii) The `ok` flag is ignored at prove.py:314.

**5. Poincaré-Miranda: CORRECT WITH FIXES.**
- The faces are right. The σ faces are thin runs over 128 z_s sub-boxes. The z_s faces use the mean value with the
  scalar product taken first, which is valid for the scalar function C₂·Γ (prove.py:250-252). box_cover.py covers the
  σ faces.
- The preconditioner is applied to the function: H = C(Γ - c).
- Fix: det C ≠ 0 is never certified (prove.py:236), and H = 0 implies Γ = c only if it is. This is a one-line arb
  check.

**6. Library traps: CORRECT.**
- arb_poly is used, not arb_series, so ctx.cap does not apply. TP truncation is explicit.
- The ctx.prec switches (600/320) only coarsen rounding.
- There is no mpmath iv on the proof path. mpmath appears only in the non-rigorous cross-check.

**7. Controls: CORRECT WITH FIXES.**
- The FRW k = 1/2 control is a real positive control of the same pipeline. The box is centred on the exact data, so a
  wrong Γ or c would break Miranda. P3 contains the closed form to 1e-77.
- But cfg_frw.json has "outer": false. The t = 0 leg (the far N-system) has no positive control, only the mpmath
  cross-check, which shares system2.json.
- N1 is non-zero by construction whenever Γ is injective. N1 and N2 are sensitivity checks, which is acceptable under
  doctrine #23, but label them so.

**8. Their doubts.**
- **(a) Saddle: CORRECT.**
  - Harada (OPENED, l.905): EC "crosses a saddle sonic point for 0 < k ≲ 0.41".
  - Carr et al. gr-qc/9902070 (OPENED, l.1084-1087): the change from the saddle region happens at ≈0.41.
  - Maison's sentence (p.3) is about k_max ≈ 0.888, so there is no contradiction at k = 1/3.
  - I recomputed: c₂ = λ+λ- = -0.002432172036452878681 matches the stored value, and the ratio is
    -0.0901059705079427750186598448076.
- **(b) Identity: WEAK, needs a scoped claim.**
  - EC (gr-qc/9402041 l.138-139, OPENED) quotes a, m/r and Ω as asymptotic values for ξ → ∞, i.e. t = 0. Its Fig. 1
    caption (l.303) defines the plotted quantity as 4πr²ρ, the same as B5's Ω.
  - a = 1.0653 rounds to EC's 1.07.
  - The gaps are +0.25% (m/r) and +4.2% (Ω), recomputed.
  - It is not a finite-ξ read-off: Ω = 9.56e-3 occurs at N = 0.0225, where m/r = 0.0609.
  - EC's pair implies V(t=0) = 0.89274, against the certified 0.89684. Unresolved; most likely 1994 numerical
    accuracy.
  - The real identification is non-certified: the first implementation gives κ = 2.8105525488 against KHA's
    2.81055255 (LOG.md:7), and V0 agrees to 20 digits (runs3/RESULTS-E3.txt:33).
- **(c) Non-uniqueness: GAP, correctly disclosed.**
  - Harada l.390 says a node admits a one-parameter family along the primary direction.
  - The float scan finds sign-change candidates on the node segment (run_branch_scan2.out: strong at z_s 0.552-0.576;
    weak at 0.536-0.576).
  - None is excluded rigorously.

**Process note.** There is no sealed BAR for B5 under LIGHTBALL (only BAR-LB-E..E4). All agreement with E3 is
post hoc. That is acceptable for a blind replication, but label it POST HOC.

## 9. Strongest honest claim
Computer-assisted proof (arb ball arithmetic, plus three named classical facts written out): for k = 1/3, the
self-similar Einstein-Euler system has a solution with:
- a regular centre (σ0 in 8.9013232777556435191172160567756856 ± 5.31e-35);
- a real-analytic crossing of a saddle-type sonic point at z_s = 0.49718669611129068678126112853988237484 ± 2.61e-39,
  with V_s = 0.11243940138809833622376302760313511666 ± 6.18e-39;
- exactly one zero of V between the centre and the sonic point;
- a continuation to t = 0 with m/r = 0.059451404935530220048167641 ± 3.90e-28 (once fix 4(i) is done).

It agrees with an independent implementation to 20 digits. It is the solution KHA and Maison computed. Its identity with
EC's printed numbers holds only to 0.25% (m/r) and 4% (Ω), and uniqueness is not claimed.

## Blocks release
1. Outer sonic piece checks on [0, θ_o] (4(i)).
2. Certified det C ≠ 0 (5).
3. P < 0 on the centre step (3).
4. Write out the centre and strong-unstable arguments.
5. Restore frw_check.py.
6. Scope "the EC solution" as in 9.

## Single next strike
One new script (no edits to sealed or existing files) that loads the EC sonic branch and C, and certifies:
- S/θ ≠ 0, P < 0, μ < 1 and V > 0 on θ ∈ [0, θ_o];
- det C ≠ 0;
- P < 0 on [0, 0.05], from the cached centre step.

It needs only seconds of compute.
