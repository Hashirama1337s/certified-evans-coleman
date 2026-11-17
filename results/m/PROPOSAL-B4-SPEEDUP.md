# PROPOSAL-B4-SPEEDUP: a cheap certified mode count (speedup seat, input for BAR-LB-M)

Speedup seat, 2026-11-16. Every TM number is a NON-COUNTING smoke (logs in `b4speed/`). No count was run. Only new files: `engine/cert/lin_tm_v3.py`, `engine/cert/fw_cert2.py`, `b4speed/*`. Sealed, pinned and build-seat files are untouched (count_b4.py, lin_*_v2.py, fw_cert.py are imported, not edited).

## 0. Bottom line
| item | before | now | log |
|---|---|---|---|
| TM remainder, 20i, r = 0.1 (same settings: nsub 2, m 16, J 4) | 1.65e-2 (v2 frame) | **1.43e-7** (characteristic frame) | b4speed/rt3_a.log, rt3_probe.log |
| largest usable disc on the imaginary axis, 40 < abs(kappa) <= 55 | r = 0.025-0.03 | **r = 0.5** (ratio <= 2.1e-3 at 50i, 55i) | rt3_b.log |
| Far Wall start K_a | 55 | **40, CERTIFIED** (40 <= abs(kappa) <= 55, plus the existing 55-2000 certificate) | runs6/fw_cert2_Ka40_K55.log |
| projected count, half-disc R1 = 40 | 13.8 h (v2 rule, R1 40) / 21.5 h (R1 55) | **1.5 h on 1 process, 0.76 h on 2** (99 discs) | b4speed/plan_v3_55s.log |
| with route 3 (Rouche on the arc, design only) | | 0.82 h on 1 process (54 discs, axis only) | plan_v3_55s.log |

Recommended pair: **R1 = 40, frame 'A1' (lin_tm_v3), nsub 2, m 16, J 4, prec 128, 8-chord circumscribed arc: 99 discs, about 1.5 h on 1 process, 0.76 h on 2 (plus about 10% for the point-value arg sharpening).** A fallback that needs no new Far Wall: R1 = 55, 134 discs, 2.05 h / 1.02 h.

## 1. Route 1: characteristic frame for the kappa-disc remainder (lin_tm_v3.py)
**Method.** v2 encloses the disc transition matrix as P_0 + sum_j |P_j| r^j + tail, with the absolute values taken entrywise in the original variables (Abar, Nbar, omega, V). The d-dependence of P is, to first order, d times the integral of Phi A1 Phi. In the eigenbasis of the large-kappa principal part A1 this is diagonal. So v3 does this per step s:
- T_s = [[I, 0], [-B^{-1}C, U]] at the step midpoint, where A1 = [[0, 0], [C, B]] and U holds the eigenvectors of B. Columns 3-4 are the characteristic directions mu_+-; columns 1-2 span the kernel of A1. This is the slaving structure of PROOF-B4-v4 sec. 2 (w_+-, z_-), frozen per step. T_s is a float matrix used as an exact matrix; any invertible T is rigorous.
- The remainder set is written in T_s coordinates. One step encloses N(d) = T_{s+1}^{-1} P(kc+d) T_s as N_0 + sum_j csym(|T_{s+1}^{-1} P_j T_s| r^j) + T_{s+1}^{-1} csym(tail) T_s. The absolute values are taken AFTER the change of variables; that is the whole change. The forcing box enters as T_{s+1}^{-1} f. Then the usual weighted-QR Lohner step. At the end T := I.
- Everything else is v2 unchanged: exact d-hierarchy polynomial part, C1 derivative matrices, Cauchy tail, t-truncation, omitted d^(m+1) forcing, and rounding radii (lin_tm_v3.py docstring).

**Result (remainder ratio R/abs(c0), lower is better; nsub 2, m 16, J 4).** c0 agrees with v2 to every printed digit at every point.
| centre | v2 (BUILD2 sec. 3 + r2_rt55.log) | v3 frame A1 | log |
|---|---|---|---|
| 1i | | r 0.25: 2.5e-7 (v2 1.2e-3); r 0.5: 2.7 (refuse) | rt3_d_small.log, rt3_b.log |
| 2.5i | | r 0.4: 9.6e-5 | rt3_d_small.log |
| 5i | r 0.2: 3.4e-4 | r 0.5: 1.6e-3 | rt3_b.log |
| 10i | r 0.2: 2.6e-2 | r 0.5: 5.4e-4 | rt3_b.log |
| 20i | r 0.1: 1.65e-2 (rerun here: identical) | r 0.1: 1.4e-7; 0.3: 5.5e-6; 0.5: 4.7e-4; 0.7: 5.5e-2; 1.0: 310 | rt3_probe/a/b.log |
| 30i | | r 0.5: 6.5e-4 | rt3_b.log |
| 35i | | r 0.6: 6.8e-3 | rt3_e_arc.log |
| 40i | r 0.05: 1.3e-2; 0.1: 2.4 | r 0.05: 3.6e-7; 0.1: 5.9e-7; 0.3: 1.7e-5; 0.5: 9.9e-4 | rt3_a.log |
| 50i | r 0.03: 1.7e-3; 0.05: 0.29 | r 0.5: 1.7e-3 | rt3_b.log |
| 55i | r 0.03: 4.2e-4; 0.05: 0.23 | r 0.5: 2.1e-3; 0.7: 0.15 | rt3_b.log |
| 5+35i | | r 0.7: 4.5e-7 | rt3_e_arc.log |
| 10+38i | | r 1.0: 2.2e-6 | rt3_e_arc.log |
| 20+40i | r 0.3: 2.2e-3 | r 1.0: 9.0e-6; 1.5: 2.5e-2 | rt3_b/e.log |
| 30+30i | | r 1.5: 2.6e-2; 2.0 (m 24): 8.6e4 | rt3_e_arc/f.log |
| 40+20i | r 0.3: 8.2e-7 | r 1.5: 4.4e-3; 2.0 (nsub 1): 3.0e3 | rt3_e/c.log |

What the table says:
- On the imaginary axis the remainder is now **flat in abs(kappa)** from 10i to 55i at r = 0.5 (5e-4 to 2e-3).
- It still grows with r: at 20i it goes 1.4e-7, 5.5e-6, 4.7e-4, 5.5e-2, 310 for r = 0.1, 0.3, 0.5, 0.7, 1.0. The slope dropped from about 24 decades per unit r (v2) to about 10 between r = 0.3 and 0.7, and the level dropped by 5 decades.
- m = 24 does not change r = 0.7 at 20i (5.53e-2 vs 5.46e-2, rt3_f_m24.log), so the polynomial degree is not the limit. The rest is the remaining frame growth, not truncation.
- The 'full' frame (eigenvectors of A0 + kc A1) is worse than 'A1' at 1i: 1.3e-6 vs 2.5e-7 at r 0.25 (rt3_d_small.log). 'A1' is the default.

**Speed.** prec 128 instead of 200 gives the same numbers to every printed digit at 20i and 55i, r 0.5. It takes 54-56 s per TM instead of 68-75 s (rt3_g_p128.log; ball arithmetic is rigorous at any precision). J = 3 is about 50 s with ratios 14-29% larger (rt3_g_p128J3.log): kept J = 4. nsub 1 is slower (111-115 s, mean order 54; rt3_c_nsub1.log) although tighter at 55i (7.7e-5): kept nsub 2.

**Exact step cache (b4speed/cacheL.py).** Linear2 setup (384-466 s) is done once and stored. Mids are bit-identical (man/exp integers). Radii come back rounded UP by a relative amount below 1e-8, so the reload is an outward enclosure. The script's round-trip check counts these rad differences as "mismatches": 862 (nsub 2), 485 (nsub 1); a unit test showed equal mids and larger rads. Point D(0) and D(20i) agree to every printed digit with the freshly built object (cacheL_build.log, cacheL_build1.log). The v2 TM on the cache reproduces BUILD2's 20i r 0.1 number exactly (1.65e-2, rt3_a.log). Reload takes 7 s. The cache is for speed only: a sealed count may rebuild instead (+6 min).

## 2. Route 2: lower the Far Wall start (fw_cert2.py)
fw_cert2 reuses every fw_cert object (bulk_ia, centre_theta, centre_test_ia, cells, centre data) and changes only discretisation knobs:
- M per band;
- an ns ladder (log-cells of the centre Volterra majorant): each failing theta-cell is retried with finer ns;
- optional theta bisection and K-cell ratio.

Every accepted (theta-cell, K-cell) pair passes exactly fw_cert.centre_test_ia.

| probe (NON-COUNTING unless marked) | worst lower/rhs at the theta = pi/2 cell (needs > 1) | log |
|---|---|---|
| Ka 35, M 3.5, ns 150 | -0.27 (8 of 9 sampled cells fail) | fw2_probe_Ka35_ns150.log |
| Ka 45 / 50, M 3.5, ns 150 | -0.07 / 0.41 | fw2_probe_Ka45/50_ns150.log |
| Ka 50, M 3.5: ns 300 / NB1 6000 / M 3 / M 4 | 1.40 / 0.44 / 2.00 / -0.09 | fw2_probe_Ka50_*.log |
| Ka 40, M 3, ns 300 (9 cells) | 0.84 (only the axis cell fails) | fw2_probe_Ka40_M3_ns300.log |
| Ka 40, axis cell: M 2.5 ns 300 / M 3 ns 600 / M 3 ns 300 kratio 1.01 | 1.064 / 1.052 / 0.88 | fw2_ax_Ka40_*.log |
| Ka 35, M 3, ns 300 (9 cells) | 0.21 (only the axis cell fails) | fw2_probe_Ka35_M3_ns300.log |
| Ka 35, axis cell, M 2.5 ns 600 | 0.52 | fw2_ax_Ka35_M2.5_ns600_kr1.03.log |
| Ka 37, axis cell, M 2.5, ns ladder 300/600/1200 | 0.80 (fails at 1200) | fw2_ax_Ka37_M2.5_ladder1200.log |
| **Ka 40 -> 55, M 2.5, 160 theta-cells, ns ladder 150/300/600** | **CERTIFIED**: 11 K-cells, worst Collatz 0.7566, e_b <= 0.3709, nbar <= 5.635; 157 cells at ns 150, 3 at ns 300, no bisection, no failure; worst lower/rhs 1.011, worst margin 7.5e-4; 756 s | **runs6/fw_cert2_Ka40_K55.log** |

With runs6/fw_cert_Ka55_K2000.log this gives **no physical eigenvalue for Re kappa >= 0, 40 <= abs(kappa) <= 2000 (CERTIFIED, interval)**. M is a free matching parameter, so bands with different M (2.5 below 55, 3.5 above) combine.
- Caveats: the margin is thin (lower/rhs 1.011 at theta in [pi/2 - pi/320, pi/2], K in [40, 41.2]). The new runs need the same referee as fw_cert (sections 2-5 are unrefereed). The tail abs(kappa) > 2000 stays OPEN.
- Below 40: discretisation stops paying (Ka 37: ns 1200 still 0.80). The binding cell is always theta = pi/2, K at the band floor. The **25 route** (exact centre series) is not uniform in K. For a finite band [K_lo, 40] it would not need uniformity, but it does need the centre solutions transported from xbar = -4.6 to x_c(K) at kappa-cells: the same kappa-disc transport problem as the count. So it is no cheaper than counting the band directly. Not pursued.

## 3. Route 3 (design only): Rouche on the large-abs(kappa) arc
Goal: on abs(kappa) = R1 >= K_a, show D(kappa) = D_mod(kappa)(1 + eta) with abs(eta) < 1, where D_mod is explicit. Then Delta arg D along the arc = Delta arg D_mod + (a term pinned by the certified endpoint values, abs < pi), with no discs. Ingredients, all present in fw_cert:
1. Bulk: row 1 of the bootstrap is based at x_c, so w_+(0) = E_+(0, x_c) w_+(x_c)(1 + eta_1) with abs(eta_1) <= v_1 - 1, where v = (I - A)^{-1} e_1 (already computed per K-cell). Need v_1 < 2 (to check).
2. E_+(0, x_c) = exp(-kappa I_1(K) + I_2(K)), with I_1 = int dx/lambda_+ and I_2 = int b_++/lambda_+ over [x_c(K), 0]. These are kappa-free cell sums, so the arg change along the arc is -R1 Delta(sin theta) I_1 plus an enclosure.
3. Sonic data w_+(0) for the (Abar, Nbar)(0) = (0, 1) solution: polynomial in kappa (the pole-free Cramer form, lin_tm.py:196-198).
4. Centre: rho in Y^f(Z eps) + B1/K +- E2/K^2 per theta-cell (closed-form flat kernel).
5. The mu-disc margin lower > e_b |w_+| is exactly abs(det_model) > abs(error) at x_c, so it bounds eta_2.
6. Transport x_c -> xbar for the TM normalisation: the Wronskian factor exp(int tr A0 + kappa int tr A1), explicit, plus the change of variables (Abar elimination factor kappa - A).

Cost: one extra pass of fw_cert data per theta-cell (seconds). Gain: the arc goes to zero discs. At R1 = 40 that is 99 -> 54 discs (0.82 h on 1 process). Risk: v_1 < 2 is not yet checked. A thin arg enclosure per theta-cell is needed so that consecutive cells pin the winding, so 160 cells may be too coarse.

## 4. Cost table (upper half-disc path, frame A1, nsub 2, m 16, J 4, prec 128, 55 s per TM, b4speed/plan_v3_55s.log)
Radius rule (Re kappa < 5): r = 0.25 for abs(kappa) <= 2, r = 0.4 for <= 5, r = 0.5 above. Re kappa >= 5: r = 1.0. The arc is a circumscribed 8-chord polygon (vertices at 1.005 R1). count_b4.N_ARC = 64 forces at least 64 arc discs: with 64 chords R1 = 40 needs 125 discs (1.91 h).

| R1 | 20 | 30 | 35 | **40** | 55 |
|---|---|---|---|---|---|
| Far Wall available at R1? | no | no | no | **yes (fw_cert2)** | yes (fw_cert) |
| discs (axis part) | 57 (29) | 78 (42) | 91 (48) | **99 (54)** | 134 (73) |
| 1 process | 0.87 h | 1.19 h | 1.39 h | **1.51 h** | 2.05 h |
| 2 processes | 0.44 h | 0.60 h | 0.70 h | **0.76 h** | 1.02 h |
| route 3 (axis only), 1 / 2 processes | 0.44 / 0.22 h | 0.64 / 0.32 h | 0.73 / 0.37 h | 0.82 / 0.41 h | 1.12 / 0.56 h |

Add about 10% for count_b4's point-value arg sharpening (one point value per disc) and for refinements. The v2 projection at R1 = 55 was 1039 discs, 21.5 h (PROPOSAL-B4-BUILD2 sec. 4). The same plan code with the v2 rule reproduces 1039 here.

## 5. What the count needs (for the build seat; nothing edited)
- count_b4.count builds LT2.LinearTM2. Swap in lin_tm_v3.LinearTM3(L, m=16, J=4, frame='A1'), and set ctx.prec = 128 after Linear2 is built.
- Make N_ARC = 8 and the radius rule above configurable (R_RULE = [(2, .25), (5, .4), (90, .5)], R_RIGHT = 1.0).
- R1 = 40 with the far-wall note "40 <= abs(kappa) <= 2000: fw_cert2 (40-55) + fw_cert (55-2000)".
- Optional: CachedLinear2 from b4speed/cacheL.py to skip the 6-min setup (outward-rounded reload).

## 6. Open risks
1. **Radius-rule coverage.** Measured on the axis at 1, 2.5, 5, 10, 20, 30, 35, 40, 50 and 55i, and off-axis at 5+35i, 10+38i, 20+40i, 30+30i (1.5) and 40+20i (1.5). NOT measured: kappa in [0, 0.5i]; the arc points near the real axis (40 + 0i at r 1); and anything on the R1 = 55 arc. 1i jumps from 2.5e-7 (r 0.25) to 2.7 (r 0.5), a sharp cliff, so small abs(kappa) must keep r <= 0.25-0.4. Pre-registered refinement (count_b4 REFINE 3) covers misses at extra cost.
2. **Far Wall 40-55 has a thin margin** (1.1%, at the imaginary axis), on unrefereed v4 mathematics. A referee finding that loosens any constant would put K_a back to 55 (cost 2.05 h, still fine).
3. The r-growth is not gone: 20i goes from 4.7e-4 at r 0.5 to 310 at r 1.0. Larger discs need a better frame. Two candidates: a kappa-dependent slaving frame (p = sigma/kappa + chi/kappa^2), or T_s varying within the step (the T' term). Neither is built.
4. The cache radii are rounded up. Mids are exact, so this is rigorous, but cacheL's "mismatch" count reports rad differences and should be relabelled before pinning.
5. Controls (C1-C3, N1-N3) have not been rerun with LinearTM3. Before sealing, at least C1a/C1b, N2 and the k=1/2 C3 must pass with frame A1 (frame A1 at k=1/2 is untested).
6. Route 3 is design only.

## 7. Commands (from <LAB>, all under memcap)
```
py -3.11 <MEMCAP> 8 py -3.11 -u b4speed/cacheL.py build 1/3 2 1e-20 b4speed/L13_nsub2.pkl     # 474 s
py -3.11 <MEMCAP> 8 py -3.11 -u b4speed/rt3.py b4speed/L13_nsub2.pkl 2 "A1,0,55,0.5;v2,0,20,0.1"   # env LB_PREC, LB_M, LB_J
FW2_NS_LADDER=150,300,600 py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/fw_cert2.py run 2.5 40 55 160 150 1500 2   # 756 s
LB_NARC=8 py -3.11 <MEMCAP> 2 py -3.11 b4speed/plan3.py 55 "2:0.25,5:0.4,90:0.5" 1.0 20,30,35,40,55
```

## Round 2 (lead requests: count_b4v3, radius gaps, controls with TM3, cap-leak audit, relabel, pins)
All NON-COUNTING. Logs in `b4speed/r2/`. Only new files (count_b4v3.py, b4speed/r2/*) plus the relabel inside b4speed/cacheL.py. lin_tm_v3.py and fw_cert2.py are untouched in this round, and PROOF-B4-v4 was not edited.

### R2.1 `engine/cert/count_b4v3.py`
- Imports count_b4 unchanged (it stays as pinned v2). It substitutes two engines: `Linear2` becomes `CachedLinear2` (k = 1/3, nsub 2) and `LinearTM2` becomes `LinearTM3(frame 'A1')`. It also sets `ctx.prec = 128` after setup and **`ctx.cap = 32` before any import** (see R2.4).
- Configurable through env vars: LB_R1 (default 40), LB_NARC (8), LB_RULE (`2:0.25,5:0.4,90:0.5`), LB_RRIGHT (1.0), LB_FRAME (A1), LB_PREC (128), LB_CACHE13.
- Modes:
  - `plan`: 99 discs at R1 = 40, matching plan3.
  - `count`: count_b4's counting code, called as is.
  - `kstar`.
  - `smoke <k> <parts>`: runs lin_b4_smoke3.py with the v3 engine.
- Caveat: the smoke banner "prec 200" comes from lin_b4_smoke3's own PREC variable. The effective precision is the `[v3] ... prec set to 128` line.

### R2.2 Radius-rule gaps closed (b4speed/r2/gaps.log; frame A1, nsub 2, m 16, J 4, prec 128). Every ratio <= 1.6e-3.
| centre, r | ratio R/abs(c0) |
|---|---|
| 0, 0.25 | 4.9e-9 |
| 0.25i, 0.25 | 1.0e-7 |
| 0.5i, 0.25 | 1.9e-7 |
| 2i, 0.25 | 5.1e-8 |
| 2.2i, 0.4 | 2.1e-4 |
| 4.8i, 0.4 | 5.0e-5 |
| 5.2i, 0.5 | 1.6e-3 |
| 2.5+40.1i (last chord near the axis), 0.5 | 2.3e-6 |
| R1 = 40 arc vertices (radius 40.194, angles k*pi/16, k = 0..7), 1.0 | 1.7e-7, 2.9e-7, 6.7e-7, 4.0e-6, 1.0e-5, 9.9e-6, 4.6e-6, 3.0e-6 |

The 1i cliff (2.5e-7 at r 0.25, 2.7 at r 0.5) lies inside the |kappa| <= 2 zone, which uses r 0.25. The rule stands unchanged: r = 0.25 for |kappa| <= 2, 0.4 for <= 5, 0.5 above, and 1.0 where Re kappa >= 5.

### R2.3 Controls with TM3 through count_b4v3: all PASS
| control | result | v2 (build seat) | log |
|---|---|---|---|
| C1a D_HKA around kappa_g, r 0.08 | winding 1 PASS, 4 discs, 155 s | PASS, 222 s | r2/c1c2.log |
| C1b D_new | winding 0 PASS, 154 s | PASS | r2/c1c2.log |
| C1c real bisection | 0.3556992037 +- 7.41e-11, overlaps closed-form kappa_g | same | r2/c1c2.log |
| C2 deflation identity at 0.5, 1.5, 2+i, 4-2i | residual contains 0 at all four | same | r2/c1c2.log |
| C3 k = 1/2 (fresh Linear2, nsub 1, 66 splits) | kappa_g contains 15/17; C3a winding 1 PASS (4 discs, 502 s); C3b 0.88235 +- 4.99e-6 contains 15/17 | PASS | r2/c3.log |
| N1 new, box [3.5,8]x[-4,4], r 0.4, all edges | Delta/pi = 0 +- 6.5e-8 (unsharpened +- 7.4e-6), winding 0 PASS, within 0.25, parity True, 22 discs, 986 s | only partial edges run | r2/n1.log |
| N2 D_HKA/(kappa+kappa_g), [0.1,3.5]x[-1,1], r 0.3 | 2 +- 2.2e-8 (unsharpened +- 0.011), PASS, parity True, 14 discs | 2 +- 0.67 | r2/n2.log |
| N3 e_Abar, [0.5,3.5]x[-1,1], r 0.3 | 2 +- 3.7e-9 (unsharpened +- 5.8e-7), PASS, parity True, 13 discs | 2 +- 5.7e-3 | r2/n3abar.log |
| N3-ref e_V, same box | 1 +- 1.7e-8 (unsharpened +- 5.8e-7), PASS, parity True, 13 discs | 1 +- 2.2e-3 | r2/n3ref.log |
| Dry run: count_b4v3 count, R1 = 3 box (LB_SHAPE=box), the real count code with ledger + sharpening | **N = 1**, Delta/pi = 1 +- 3.2e-8, parity (-1)^N = sign D(0) sign D(R1): True, 17 discs, no refinement, 669 s (37-40 s per disc incl. point value) | v2 dry run paused after 11 discs | r2/dry_R3.log, r2/dry_R3.jsonl |

C1/C2 ran at the default cap, before the cap line was added. The containment test in R2.4 shows that every object they use is identical at cap 32. Everything else ran with cap 32.

### R2.4 Cap-leak audit of the count stack (referee B1 class)
**Static check.** I grepped for flint series objects (arb_series, acb_series, *_poly, ctx.cap, integral, derivative, compose) in lin_stages.py, lin_stages_v2.py, lin_tm.py, lin_tm_v2.py, lin_tm_v3.py, lin_mc.py, count_b4.py and count_b4v3.py. **None of them builds a flint series**: the jets are Python lists (core.Tape.jet_order, smul/sdiv/horner) with explicit orders. The only flint series in the chain sit in the tube build:
- param.tail_bound, rebound to driver_v31.arb_series_exact (prec 2^16; R7);
- core.Tape.series_eval and driver_v32, which both set ctx.cap explicitly before building series and restore it after.

**Dynamic containment test** (`b4speed/r2/capchk_count.py`, ctx.cap set before every import). Each object rebuilt at cap 32 is compared with the as-run object:
| object | result | log |
|---|---|---|
| tube k = 1/3 (lin_tube_v2.build; the pinned pkl is not written, hash unchanged f6bff9e2...) | 89677 serialized arbs + 728 other fields: **all identical** | r2/capchk_tube13.log |
| tube k = 1/2 (hash unchanged 17dadbfc...) | 146446 arbs + 1165 fields: **all identical** | r2/capchk_tube12.log |
| Linear2 stage L, k = 1/3 nsub 2, ALL 1440 steps, every coefficient | 4,128,717 identical, 1,498,803 radius rounded up (the cache's outward rounding), 0 mid differences, 0 not contained | r2/capchk_linear13.log |
| stages S + C (4495 arbs, exact mid/rad hash) | cap 10 and cap 32: same sha256 da9a457e... | r2/capchk_sc.log |
| D_tm(20i, r 0.5), frame A1, prec 128, at cap 32 | c0 and R identical to the as-run value (R 0.00101) | r2/capchk_linear13.log |

**Verdict.** The count stack does not leak through ctx.cap; the as-run objects equal the cap-32 objects. The fix is `ctx.cap = 32` at the top of count_b4v3.py: valid (code that sets the cap restores it to 32) and number-neutral by the test above.

**Not covered:** fw_cert2.py imports fw_cert.py, so **runs6/fw_cert2_Ka40_K55.log inherits B1 and is VOID as a certificate** until the proof seat's cap-fixed rerun of `run 2.5 40 55` passes. R1 = 40 depends on that rerun. If it fails, the fallback is R1 = 55, which needs fw_cert_Ka55_K2000 itself rerun under the fix (also in the proof seat's queue).

### R2.5 Relabel
cacheL.py now reports "identical / radius rounded up (outward, allowed) / mid differs (ERROR) / radius smaller (ERROR)" through `roundtrip()`, plus a `check` mode (fresh build vs an existing cache). The full-stage-L comparison above uses the same categories. The round-1 lines "862 / 485 mismatches" were all radius-rounded-up: mids identical, consistent with the full check.

### R2.6 Final plan
R1 = 40 half-disc, 8-chord circumscribed arc, frame A1, nsub 2, m 16, J 4, prec 128, cap 32, rule as in R2.2: **99 discs**. With 55 s per TM plus about 10% for point values and refinements: **about 1.7 h on 1 process, 0.85 h on 2** (1.51 / 0.76 h without the margin). The R1 = 3 dry run measured 37-40 s per disc including the point value, so this is conservative. Fallback R1 = 55: 134 discs, 2.05 h / 1.02 h (+10%). Precondition: a cap-fixed Far Wall on [40, 55] (or [55, 2000] for the fallback).
Command: `py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/count_b4v3.py count 40 runs7/count_R40.jsonl`. For 2 processes, split the plan by edge (count_b4 LB_MAXDISCS/ledger resume, or two ledgers by edge). The split is not built.

### R2.7 Files to pin (sha256, b4speed/r2/pin_sha256.txt)
New:
```
54c1e2f7f2defad39483f3de6af9774212b2b78da2d9f416234598255708f4a3  engine/cert/lin_tm_v3.py
29abac07b0ea07c27a4ebf79303e1018bf690261f8be94109d52401e7c92ee07  engine/cert/count_b4v3.py
a0489e264e31f61a9934b0089644025048adcc304c7e2401031c61c9e8b98533  b4speed/cacheL.py
c239cbe5347059080ea3c4144b8c5a65133cb4c4812a280096c8f5116dfbb782  b4speed/L13_nsub2.pkl
efa020736bcafd903e34031077a2d50016ec8a7bf71c4fb4f2bdeb46fe47b382  engine/cert/fw_cert2.py        (Far Wall [40,55]; pending the cap-fixed rerun)
51b5d2be5cf33a71372373321629fa8149a06994c4fb06d1d9e0668b2e8cfe81  b4speed/r2/capchk_count.py     (control)
8a8d24a996f57c4f1f145009859113c067e9ac9e10cfab3f9313d64728e3f59a  b4speed/rt3.py                 (RT smoke tool)
8a0b27c6b9aedbd2dd084ca1b3da2e98222272147677ab2ab332d4d8c1fd6509  b4speed/plan3.py               (plan only)
```
Carried (imported unchanged):
```
7f17adc9c0b10ec529fba28ea356296b7b3bd77dd95bd312e3bb791c4d2c7204  engine/cert/count_b4.py
e477ec1f6fc3f36ec16e776c26c0516e1fa5b964ee769286d8e0c290f20bf2d5  engine/cert/lin_stages_v2.py
582e24236c8ff81594cb56a4c2afdb0edc9a083dfdc261cf61623f4679deb821  engine/cert/lin_stages.py
80b6a52355e28605a8f6d72784fa025ddac354733c413191cb086f84f60985d1  engine/cert/lin_tm_v2.py
98a06289ee5db85f16c7610968754b4a9c5a38417ec4efe1cf2a5273bf6411c9  engine/cert/lin_tm.py
a39312d4a76008e638e850ccb8417a0261b5c6ea590619c973b143a15d73a368  engine/cert/lin_tube_v2.py
f6bff9e203c8a4e054a35bff0453d179d570dca6bd520665b900cb2b2c3071df  engine/cert/lin_tube2_k1_3.pkl
17dadbfcfb78e235010610e0b7b627dcb76a35ceedc6601e4fc1b845bad3badf  engine/cert/lin_tube2_k1_2.pkl
dce8eba80ec94e1dff72d4187f17029cd3a0e62ae65f3e71e41658ed7761fc26  engine/cert/lin_fields2.pkl
ee336f0e6185a620c9d995466af69b486372716dc1f593ed56d0fc13bbcbf64d  engine/cert/lin_mc.py
45657bd361d2c18b7a2476e5f232f8b3c16f73835e8dc0618242fac9bc5f68cd  engine/cert/lin_b4_smoke3.py   (controls)
932a32b5f6b23f416147ce112b3a51e2d840b855fea79b339fdc397d3ca88f85  engine/cert/core.py
5d3321a150c4b55db1aa4052f4197e6583d6145c2b2ad257905d24dff3752fc4  engine/cert/fw_cert.py          (B1 leak; superseded by the proof seat's fix)
```
Also carried: the E5 pins (driver_v32.py, driver_v31.py, ... per BAR-LB-E5).

### R2.8 Remaining risks
1. The Far Wall [40, 55] (and [55, 2000]) is void until the cap-fixed rerun. The R1 = 40 plan depends on it.
2. The 2-process split is not built.
3. Untested: the A1 frame at k = 1/2 for large |kappa| (C3 only exercises small |kappa|).
4. The sharpening relies on point values: cheap, already in the measured dry-run time.

## Round 3: REFEREE-B4FW2 sealing fixes in the count wrapper (b4speed/r3/, NON-COUNTING)
New file `engine/cert/count_b4v3_seal.py`, wired in by `count_b4v3.py`. count_b4.py is unedited; its upper_path, disc_plan and count are replaced in the running process only.

| referee fix | implementation | evidence |
|---|---|---|
| 1 arc vertices | Rp = R1/cos(dphi/2)·(1 + 1e-12). Exact Fraction checks on the float vertices actually used: min chord distance² − R1² ≥ 0, max vertex ≤ 55. Plan rebuilt with each edge's last piece snapped to the edge end. Exact asserts: plan[i].b == plan[i+1].a, first a == (Rp, 0), last b == 0, every piece ≤ 2r. | R1 = 40: min dist² − 1600 = +3.2e-9, vertex 40.1935428951, 99 pieces continuous (r3/pre_R40.log). R1 = 3: +1.8e-11, 15 pieces (r3/dry_hd_R3.log). |
| 2 parity | `(-1)**N != sign D(0) sign D(Rp)` raises ArithmeticError | |
| 3 preconditions (ball evaluations, printed) | c_V(x_m) = Irow[3] on the certified tube-end box at x = xbar; Ā_ss'(0) = L.abar_x0; detA root = −dA(0)/(dA(1) − dA(0)), with the second difference printed; iw == 2 at all plan points, plus run-time guards in sonic_tm, detA and sonic_start | **c_V(x_m) = 0.19288155165634 ± 1.41e-15**. Ā_ss'(0) = −0.198928186801772 ± 4.6e-16. detA root = **−0.0990290820550 ± 2.9e-14** (referee: −0.09902908206 ± 5e-12); second difference ± 1.5e-13. iw == 2 at all 199 R1 = 40 plan points and in all 49 run-time calls of the dry run. |
| 5 ledger / env | Fresh ledger starts with a header line (R1, shape, N_ARC, rule, M, J, nsub, tol, REFINE, plan size, frame, prec, cap, sha256 of 13 engine files incl. cache and tube); resume asserts it. LB_R_FIXED and LB_MAXDISCS refused in count mode. Cache sha256 verified at load (c239cbe5...). All settings printed. | header refusal (LB_RRIGHT changed): refused (r3/refuse_header.log); LB_MAXDISCS: refused (r3/refuse_maxdiscs.log); resume with matching header: 15 discs from record, N = 1 again (r3/resume_check.log) |

**Positive control (next strike): half-disc dry run R1 = 3**, LB_SHAPE=halfdisc, the exact sealed code path. Result: **N = 1**, Delta/pi = 1.0000000 ± 3.4e-8, parity True, 15 discs, 0 refinements, 834 s (r3/dry_hd_R3.log, ledger r3/dry_hd_R3.jsonl). D(0) = −108.527 ± 2.8e-4, D(Rp = 3.0145) = 0.653570 ± 3.7e-7, so κ* ≈ 2.81 lies inside.

**Scope wording for the sealed count:** "zeros of D_new with Re κ ≥ 0, |κ| ≤ 2000". It requires c_V(x_m) ≠ 0 (now certified), the [40, 55] and [55, 2000] Far Wall (fw_cert3 pins per REFEREE-B4FW2 Blocks-sealing 4), and the half-disc count at R1 = 40. The tail |κ| > 2000 stays OPEN.

**Final count-side pins (sha256, b4speed/r3/pin_sha256_count.txt).** These supersede the count rows of R2.7. count_b4v3.py changed in this round; the others are unchanged.
```
9d6db027f0ca12a2fefc6fcb716bfd36d3a5d422c971611403bcec04dd13b28f  engine/cert/count_b4v3.py
6e7d49f3e1af8e8f0a3caf32476adfd3557cee49cf7d5a067f528ab7d02f470c  engine/cert/count_b4v3_seal.py
54c1e2f7f2defad39483f3de6af9774212b2b78da2d9f416234598255708f4a3  engine/cert/lin_tm_v3.py
a0489e264e31f61a9934b0089644025048adcc304c7e2401031c61c9e8b98533  b4speed/cacheL.py
c239cbe5347059080ea3c4144b8c5a65133cb4c4812a280096c8f5116dfbb782  b4speed/L13_nsub2.pkl
7f17adc9c0b10ec529fba28ea356296b7b3bd77dd95bd312e3bb791c4d2c7204  engine/cert/count_b4.py
80b6a52355e28605a8f6d72784fa025ddac354733c413191cb086f84f60985d1  engine/cert/lin_tm_v2.py
98a06289ee5db85f16c7610968754b4a9c5a38417ec4efe1cf2a5273bf6411c9  engine/cert/lin_tm.py
e477ec1f6fc3f36ec16e776c26c0516e1fa5b964ee769286d8e0c290f20bf2d5  engine/cert/lin_stages_v2.py
582e24236c8ff81594cb56a4c2afdb0edc9a083dfdc261cf61623f4679deb821  engine/cert/lin_stages.py
932a32b5f6b23f416147ce112b3a51e2d840b855fea79b339fdc397d3ca88f85  engine/cert/core.py
ee336f0e6185a620c9d995466af69b486372716dc1f593ed56d0fc13bbcbf64d  engine/cert/lin_mc.py
dce8eba80ec94e1dff72d4187f17029cd3a0e62ae65f3e71e41658ed7761fc26  engine/cert/lin_fields2.pkl
a39312d4a76008e638e850ccb8417a0261b5c6ea590619c973b143a15d73a368  engine/cert/lin_tube_v2.py
f6bff9e203c8a4e054a35bff0453d179d570dca6bd520665b900cb2b2c3071df  engine/cert/lin_tube2_k1_3.pkl
17dadbfcfb78e235010610e0b7b627dcb76a35ceedc6601e4fc1b845bad3badf  engine/cert/lin_tube2_k1_2.pkl
45657bd361d2c18b7a2476e5f232f8b3c16f73835e8dc0618242fac9bc5f68cd  engine/cert/lin_b4_smoke3.py   (controls)
51b5d2be5cf33a71372373321629fa8149a06994c4fb06d1d9e0668b2e8cfe81  b4speed/r2/capchk_count.py     (control)
```
Plus the E5 pins (driver_v32/v31, ... per BAR-LB-E5) and the Far Wall pins from REFEREE-B4FW2 item 4 (fw_cert3.py 8f71e486..., lin_fw4_bg3.pkl 53008a7d..., ...).

Sealed command (NOT started):
```
py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/count_b4v3.py count 40 runs7/count_R40.jsonl
```
Use a fresh path. Projected 99 discs, about 1.5-1.7 h on one process (the dry runs measured 37-56 s per disc).
