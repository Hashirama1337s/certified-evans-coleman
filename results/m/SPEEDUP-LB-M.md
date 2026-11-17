# SPEEDUP-LB-M: speedup design pass for the B4 certified mode count (doctrine #33)

The full record is in PROPOSAL-B4-SPEEDUP.md (rounds 1-3). This file is the summary that doctrine #33 asks for.

- **Projection before:** 21.5 h on one process. That was 1039 discs at R1 = 55 with the v2 Lohner frame (PROPOSAL-B4-BUILD2.md).
- **Projection after:** 1.5-1.7 h on one process. That is 99 discs at R1 = 40 (b4speed/plan_v3_55s.log; r3/pre_R40.log).
- **Dominant cost driver:** the number of discs on the imaginary axis. The v2 frame's remainder grows like exp(c r), which forces r = 0.025-0.05 for |kappa| > 40.
- **Lever 1, change of frame:** a characteristic (A1, slaving-direction) frame, lin_tm_v3.py. Ratio at 20i, r 0.1: 1.65e-2 -> 1.43e-7. Axis discs grow to r = 0.5.
- **Lever 2, cheap exclusion bound:** the Far Wall (fw_cert3) certifies there is no eigenvalue for 40 <= |kappa| <= 2000 (refereed, REFEREE-B4FW2). This shrinks the counted region from R1 = 55 to 40.
- **Lever 3, certified asymptotics (Rouche on the arc):** designed but not built (PROPOSAL-B4-SPEEDUP sec. 3). It would give about 0.82 h. Not needed under the 2 h line.
