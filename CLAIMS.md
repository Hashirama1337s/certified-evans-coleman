# CLAIMS: every claim with its status, and every number with its source

Status labels: **CERTIFIED** = computer-assisted enclosure (outward ball arithmetic) from a sealed run; **PROVED** =
written proof checked by at least one internal hostile referee seat; **CITED** = taken from an opened source (locator in
paper.md, References); **PROBE** = numerical value without an enclosure, or from an unsealed run; **OBSERVATION** = a
remark on a printed value whose cause we did not test. Seal ids are the lab ledger's claim hashes. Referee seats are named
by their report file. Source files are lab file names; the repository maps them under `results/` and `code/`.

## Claims

| # | Claim (paper section) | Status | Seals | Referee seats, blind seats | Novelty gate (seats/GATE.md, GATE-REFRESH-2026-11-17.md) |
|---|---|---|---|---|---|
| C1 | Existence of a k = 1/3 self-similar solution, regular centre, analytic saddle sonic crossing, one V sign change, no second sonic point, continues to t = 0; V0, W0, xc enclosures; m/r(t=0) (§3, Thm 1) | CERTIFIED | d0af75803f | REFEREE-E (on E3, repairs carried into E5), REFEREE-R8, REFEREE-E5CAP (E5 SAFE) | (a) OPEN; wording "first interval-arithmetic (computer-assisted) certificate", never "first existence proof" (HKA 1996 l.2231-2234) |
| C2 | Blind second existence proof (Chebyshev + Poincaré-Miranda), σ0, z_s, V_s, t = 0 values (§4, Thm 2) | CERTIFIED | 71904980c4 (+ bee365c702) | REFEREE-B5 (CORRECT WITH FIXES; fixes = closure test) | as C1 |
| C3 | E5 and B5 enclosures of V0 / V_s and W0 / σ0 overlap (19-digit agreement) (§4) | CERTIFIED comparison, POST HOC | d0af75803f, 71904980c4 | REFEREE-B5 ("agreement ... is POST HOC") | not a claim of identity (no local uniqueness) |
| C4 | Local uniqueness | NOT CLAIMED | - | REFEREE-E (no derivative bound on tail), REFEREE-B5 (node-segment GAP) | - |
| C5 | Third self-similar (type II) critical solution with an existence proof (after Bizoń-Wasserman CSS, Reiterer-Trubowitz DSS), first for a perfect fluid (§12) | CITED (list of the Living Review §VII C) + C1/C2 | - | - | GATE-REFRESH wording audit; Bartnik-McKinnon (type I, 1991) excluded from the count by "self-similar (type II)" |
| C6 | Constraint propagation dI/dx = c_I I; sonic exponents 0,0,0,ρ(κ); centre exponents {−1,0,0,−2}; gauge Lemma 5.1; Theorem 3 (D = 0 ⇔ physical eigenvalue; D_HKA = −(κ−κ_g)D/Ā_ss'(0)) (§5) | PROVED (PROOF-B4 §1, Thm 1.6; exact symbolic checks) | preconditions certified in b1b0ee1cde | REFEREE-B4FW2 ("what is counted: CORRECT WITH FIXES"; fixes = certify c_V(x_m), detA root) | gauge deflation parallels Gonçalves's rank-two deflation (credit) |
| C7 | Far Wall: no physical eigenvalue for Re κ ≥ 0, 40 ≤ \|κ\| ≤ 55 (M = 2.5) and 55 ≤ \|κ\| ≤ 2000 (M = 3.5) (§7) | CERTIFIED (fw_cert3), PROVED (PROOF-B4-v4 + AMEND1) | files pinned under b1b0ee1cde | REFEREE-B4FW (runs WRONG as run: cap leak; fixed), REFEREE-B4FW2 (repairs CORRECT) | method: a frequency split close to ours credited to GHJS 2509.12435 (energy-based at the extremes in theirs; Levinson slaving in ours) |
| C8 | Far Wall tail: no physical eigenvalue for Re κ ≥ 0, \|κ\| ≥ 2000 (§7, App. C) | CERTIFIED (fw_tail2), PROVED (AMEND2-TAIL) | wording condition of 4971bef2bf met | REFEREE-B4TAIL (claim holds; control replaced, rerun) | - |
| C9 | Exactly one growing non-gauge mode in Re κ ≥ 0 (spherical, linear, analytic through the sonic point, regular at the centre, modulo κ_g); simple, real; none on Re κ = 0 (§6, Thm 5) | CERTIFIED | b1b0ee1cde (+ b92e091b0c, 4971bef2bf) | REFEREE-B4FW2 (count logic sound; fixes applied before sealing) | (b) OPEN; first for a CSS / fluid critical solution; DSS Choptuik done by Gonçalves 2026 |
| C10 | κ* ∈ [2.810552549 ± 3.83e-10], γ = 1/κ* ∈ [0.3558019224 ± 6.97e-11] (§8, Cor. 6) | CERTIFIED | b1b0ee1cde; blind G2 8a1c27c7fb (+ b62c6aecc9): HOLD, overlap | blind seat b5g (b5g/result_g2.json, REPORT.md) | (c) OPEN as a certificate; digits TAKEN numerically (HKA 1996 2.8105525488) |
| C11 | Blind second enclosure κ* ∈ [2.8105525488271268040, 2.8105525488272177534] (γ ∈ [0.35580192244308 ± 8.03e-15]), inside C10's. Computed on the B5 solution; if E5 and B5 coincide (not claimed), HKA's table value 2.8105525488 is correct to its last digit and two further digits are certified (13 significant digits; certified, beyond the sealed bar). The first engine alone confirms HKA's 9-digit text value 2.81055255 (enclosure [2.8105525486170, 2.8105525493830]). The first engine's floor (~1e-10) is its stored background tube width (§8) | CERTIFIED (blind engine, two runs differing in every setting) + PROBE (floor diagnosis, kstar11, non-counting) | 8a1c27c7fb | blind seat b5g (brief quoted ≈ 0.3558, disclosed); release referee seats/REFEREE-RELEASE.md (READY WITH FIXES) | GATE-REFRESH wording trap 1 (now: conditional wording on the B5 solution; the first engine confirms HKA's 9-digit text value) |
| C12 | γ is the exponent of the KHA/Maison scaling argument; the mass law is not proved (§8) | CITED (KHA, Maison); NOT CLAIMED (mass law) | - | - | GATE-REFRESH wording audit |
| C13 | No kink-mode instability by Harada's criterion (saddle sonic point certified); "stable" not claimed (Remark 7) | CITED (Harada gr-qc/0109042 abstract, l.791-792) + CERTIFIED input (C1 saddle) | d0af75803f | - | (b) analytic piece PARTIAL (Harada) |
| C14 | κ_g ∈ [0.3556992037 ± 7.41e-11] = 2 − A0 + (2/3)ω0 (§6 controls) | CERTIFIED (C1c control) | b1b0ee1cde | - | - |
| C15 | HKA's printed gauge mode 0.35699 looks like 0.356992 with a dropped digit (§10) | OBSERVATION | - | - | - |
| C16 | KHA (1995) eq. (18) third line is a misprint (fails on exact Friedmann; lines 1, 2, 4 hold) (§10) | PROVED (symbolic check by REFEREE-B5) | - | REFEREE-B5 | fairness note |
| C17 | Evans-Coleman 1994 m/r 0.0596 and Ω 9.56e-3 sit 0.25% and about 4% above both implementations (§10) | OBSERVATION; cause UNTESTED | d0af75803f, 71904980c4 | REFEREE-B5 (same definitions; "most likely 1994 numerical accuracy") | fairness note |
| C18 | Frame change: TM remainder ratio 1.65e-2 → 1.43e-7 at 20i, r = 0.1; projection 21.5 h → 1.5-1.7 h; sealed count 3786 s (§6.1) | PROBE (smoke ratios) + run time of the sealed count | b1b0ee1cde (run time) | REFEREE-B4FW (LinearTM3 CORRECT) | method, not headline |
| C19 | Absence of continuous / residual spectrum of a linearised evolution operator; nonlinear stability; non-spherical modes | NOT CLAIMED | - | - | W1 (non-spherical) OPEN, WHERE-ELSE.md |
| C20 | Lean: the Collatz step of the Far Wall, i.e. rational rows on all 133 bulk K-cells of [40, 2000] (16cb6d873a) and on the uniform tail matrix (6a2a9e1d2f), and the first sentence of Corollary 3.1 (Collatz rows ⇒ I − M invertible, (I − M)⁻¹ ≥ 0, ρ < 1, v ≤ (I − M)⁻¹e₁) (750fcfda8b); teeth control fails as intended (3078ef4ec1). Supersedes the one-cell check ae7bc39344 | PROVED (Lean; rows + implication only) | 16cb6d873a, 6a2a9e1d2f, 750fcfda8b | statement match to the first sentence of Cor 3.1 confirmed by the release referee (seats/REFEREE-RELEASE.md §3) | the entries are arb enclosures; the Far Wall as a whole is not Lean-checked |

## Numbers and their sources

Direction: `ball` = arb ball as stored; `ends` = interval ends as stored; `exact` = exact value or value quoted verbatim;
`probe` = non-rigorous; `print` = printed in a source.

| Printed | Direction | Source file (field) |
|---|---|---|
| V0 [0.1124394013880983362 ± 3.31e-20]; W0 [8.90132327775564352 ± 7.92e-18]; xc [0.212365646758945385 ± 7.95e-19]; V0 width 1.86e-20 | ball | runs6/target.log (PE1 lines); ledger settle d0af75803f |
| contraction ratio 0.00703; Ebox [1.93e-24 ± 7.05e-28]; margin [1.385e-27 ± 3.78e-31]; 194 τ-boxes, 560 x-boxes; 731 outward steps | ball / count | runs6/target.log |
| λ [4.29174946397 ± 8.77e-12]; λ_b [−0.3867122506 ± 4.23e-11] | ball | runs6/target.log (PE2 saddle) |
| m/r(t=0) [0.0595 ± 8.24e-5] | ball | runs6/target.log (PE3) |
| warm-up ratio 0.0247; Ebox [1.76e-22 ± 6.63e-26]; V0 [−0.565685424949238020 ± 5.73e-19]; W0 [0.2962962962962963 ± 2.32e-17] | ball | runs6/warmup.log |
| −2√2/5, 8/27 (Friedmann k = 1/2 sonic data) | exact | AMEND-LB-E-1.md; LOG.md |
| content nulls V0, W0, xc refuse; offset 1e-13 | exact | runs6/nullV0.log, runs6/target.log; BAR-LB-E5.md |
| sonic tail eps [1.85525e-32 ± 1.91e-38]; centre eps [3.28636e-56 ± 2.50e-62] | ball | seats/REFEREE-E5CAP.md |
| 525 head coefficients; relative error 7.08e-23 | exact | seats/REFEREE-R8.md §1 |
| E3: claimed 3.7e-59, true about 3.7e-35 | exact | ledger correct row 8aec628fd2; seats/REFEREE-E5CAP.md |
| σ0, z_s (width 5.2e-39), V_s (width 1.2e-38), m/r, Ω, a at t = 0 | ends | b5/RESULT.md tables (result_ec.json, results_summary.json) |
| λ₊, λ₋, ratio λ₋/λ₊ | ends | b5/RESULT.md (3) |
| Miranda box σc ± 1e-35, zc ± 1e-39; z_m = 0.4945; about 55 min per run | exact | b5/RESULT.md |
| B5 Friedmann box [0.2962962962962962962962962963 ± 1.38e-29]; 9√2/14 | ball / exact | b5/RESULT.md (Controls) |
| θ_o [0.031099573504505126343 ± 4.02e-22]; 8 pieces, 0 ambiguous; det C [−81.709996897608312500 ± 6.34e-20]; P/z [−1 ± 0.0158]; 20 cache objects; Ω, V at z = 1 [0.96 ± 3e-30], [−0.8 ± 3e-30]; 131 steps | ball / count | b5/close_b5_results.json (c1-c5) |
| 24/25, −4/5 | exact | AMEND-LB-B5C-1.md |
| 5e-38 (agreement of the two derived systems) | exact | seats/REFEREE-B5.md §1 |
| c_I formula; exponents {−1, 0, 0, −2}; ρ(κ) = (d − κ)/λ' | exact | PROOF-B4.md §1.1-1.3 |
| λ_−' ∈ [1.0986, 1.0995]; s0 = 0.0036260; g = 0.0896 | ends | PROOF-B4-v4.md §3 (sonic piece) |
| c_V(x_m) [0.19288155165634 ± 1.41e-15] at x_m = −4.6; Ā_ss'(0) [−0.198928186801772 ± 4.57e-16]; detA root [−0.0990290820550 ± 2.91e-14]; 199 plan points | ball | runs7/count_R40.log; runs7/verdict_M.json (PM3) |
| N = 1; Δ/π [1.0000 ± 2.42e-5]; 99 discs; 0 refinements; 3786 s | ball / count | runs7/verdict_M.json (PM1, runtime_s) |
| R_p = 40.1935428950926; min chord distance² − 1600 = 3.200e-09; 8 chords | exact | runs7/count_R40.log (seal F1) |
| m = 16, J = 4, nsub = 2, prec 128 | exact | runs7/count_R40.log header |
| C1a winding 1; C1b winding 0; κ_g [0.3556992037 ± 7.41e-11]; R1 = 3: N = 1, Δ/π [1.0000000 ± 3.42e-8] | ball / count | runs7/verdict_M.json (controls) |
| κ* [2.810552549 ± 3.83e-10]; γ [0.3558019224 ± 6.97e-11]; 67 s | ball | runs7/kstar.log; runs7/verdict_M.json (PM2) |
| floor about 1e-10 (background tube width) | probe | LOG.md (kstar11 probe entry) |
| 1.65e-2 → 1.43e-7; r = 0.5; 21.5 h; 1039 discs; 1.5-1.7 h | probe | PROPOSAL-B4-SPEEDUP.md §0-1; SPEEDUP-LB-M.md |
| [40, 55]: M 2.5, 11 K-cells, Collatz 0.7566, lower/rhs 1.010825, margin 7.48e-4 | ends | PROOF-B4-v4-AMEND1.md §0; runs6/fw_cert3_M2.5_Ka40_K55.log |
| [55, 2000]: M 3.5, 122 K-cells, Collatz 0.5682, lower/rhs 1.314537 | ends | PROOF-B4-v4-AMEND1.md §0; runs6/fw_cert3_M3.5_Ka55_K2000.log |
| tail: Collatz 0.54072, e_b ≤ 0.06146, n̄ ≤ 0.48433, lower/rhs 14.77, 160 θ-cells, x_T = −6.3303, θ_env = 0.001455, margin 1.0100000 | ends | PROOF-B4-v4-AMEND2-TAIL.md §1-2; runs6/fw_tail2_KT2000.log |
| tail controls: 147 comparisons, 0 non-overlaps; 26 non-overlaps (power); [0.9789, 1.0738] | count / ends | PROOF-B4-v4-AMEND2-TAIL.md §2b; runs6/fw_tail2_ctl.log |
| ξσ̂ = −0.357; β_+ = −1.000 ± 1.3e-3 | ends | PROOF-B4-v4-AMEND2-TAIL.md §1.2 |
| x_end −0.050000000 ± 7e-14; ℓ ≤ 5.41e-14; F = 0.727; δ ≤ 3.93e-14 | ends | PROOF-B4-v4-AMEND1.md §2 |
| 128 dropped-coefficient checks; capchk 7/7 | count | PROOF-B4-v4-AMEND1.md §1 |
| kernel check errors 4.8e-13 to 5.6e-11 | exact | seats/REFEREE-B4FW.md §A |
| Lean ae7bc39344, K-cell [55, 56.65], Collatz 0.568212 | exact | PROOF-B4-v4-AMEND1.md §4; ledger lean row |
| 4.9e-15 (arc dip before the fix) | exact | seats/REFEREE-B4FW2.md §2 |
| 2.81055255 (β ≃ 0.35580192; "exact to the last digit", l.2492-2496) in HKA's text; 2.8105525488 (HKA table, l.2545); 0.35699 | print | seats/src/gr-qc_9607010.txt l.2492-2496, l.2545, l.2517 |
| β ≃ 0.3558019; 0.35699 | print | seats/src/gr-qc_9503007.txt l.16-17, l.215 |
| γ = 0.3558 | print | seats/src/gr-qc_9504008.txt l.176 |
| m/r 0.0596, Ω 9.56e-3, a 1.07, β ≈ 0.36 | print | seats/src/gr-qc_9402041.txt l.20, l.138-139 |
| −0.25% / about 4% (EC gaps) | exact | b5/RESULT.md (comparison); seats/REFEREE-B5.md §8(b) |
| λ ∈ [2.674040, 2.674115] (Gonçalves) | print | seats/src/refresh-2026-11-17-abstracts.txt (Zenodo record) |
| k ≲ 0.41 (saddle), k ≲ 0.89 (kink-stable) | print | seats/src/gr-qc_0109042.txt abstract, l.900-908; gr-qc_9902070.txt l.1083-1087 |
| 7 referee reports; 9 seals (4 HOLD, 3 FAILED, 2 VOID; the duplicate 9592b015f0 is VOID, duplicate, never run) | count | seats/REFEREE-*.md; ledger list (LIGHTBALL) |

Tokens allowed without a result-file source (definitions and exact constants): k = 1/3, γ = 4/3 (KHA's adiabatic index),
c = 1/√3, √k; 2, 3, 2/3, −2/(3(1+k)); R1 = 40, 3, 55, 2000 (sealed contour and band limits); 1e-15 ... 1e-9 (sealed ladder);
1e-12 (vertex inflation); 3 halvings (sealed refinement); 0.40, 2.70, 1e-6 (BAR-LB-G2 constants); 2026-11-16 (search
floor).
