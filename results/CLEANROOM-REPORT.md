# Clean-room replay: certified-evans-coleman (code / results / tools side)

Build seat, 2026-11-16/17. Nothing was published. The repository was copied (plain `cp -r`, so file times were not
kept, as in a fresh clone) from `release/repo/certified-evans-coleman/` to `release/cleanroom/run1/`, and only the
documented commands were run there. Every Python process ran under the lab memory cap
(`memcap.py 8 python ...`; the cap covers child processes). One heavy process at a time.

Machine: Windows 11, Python 3.11.9, python-flint 0.9.0, sympy 1.14.0, mpmath 1.3.0, numpy 2.4.4.

## 1. Pins

Command (in `run1/`): `python tools/check_pins.py`

Result: `check_pins: ALL MATCH (84 pins, 81 patched files; 108 and 100 after G2 was added, section 4)`; `manifest: 183 files checked`. Log: `run1_check_pins.log`.
The 84 pins are every file sha256 in the ledger rows of the three seals and their amendments (seal-file rows: BAR /
AMEND files, instrument and formula-reference pins; pin rows; settle-row result files). No pinned file is missing
from the repository (`tools/pins.json`, `pinned_but_not_shipped` is empty). Path-patched copies are checked after
in-memory restoration.

## 2. Quick tier

Command (in `run1/`): `python tools/reproduce.py --tier quick` (work folder `run1/work/`).
Started 2026-11-16 23:51:47, finished 2026-11-17 00:08:14 (local time): **16.4 min**. Log: `run1_quick.log`; report
`run1/work/_out/reproduce_report.json`.

Data files: the seven pickles the quick checks need were copied from `results/pickles/` only after their sha256
matched (five ledger pins; `lin_fields.pkl` and `lin_fw4_bg.pkl` against the lab file sha256, they have no ledger pin).

| Step | Time | Status | Canonical sha256 (produced = expected) |
|---|---|---|---|
| e5_tails13 | 0.5 min | IDENTICAL-CANONICAL | `c3694218104365df2dfaf312b8ec94678b330c3da5edcbf378343d5976a73a0f` |
| e5_tails12 | 0.5 min | IDENTICAL-CANONICAL | `160e25ecddcf36bb036de3a5ec2a35494a1719a6ce7874e8d9991a694767e650` |
| e5_tails31 | 0.4 min | IDENTICAL-CANONICAL | `91907ef2fa47aa2c35e482e42d6fd8b4b602b1623e60dab9019263642f9f1f86` |
| e5_struct | 0.0 min | IDENTICAL-CANONICAL | `3c5a403b8ce6281cb83d568396b35173aad89df28646827b88f9692cc5ae7de6` |
| b5c (log) | 5.2 min | IDENTICAL-CANONICAL | `c65d32dc800e5391cdbf27e3c446a58e8bdc33ec9a104ccaed53ffa5a43e68f9` |
| b5c (`close_b5_results.json`) | | IDENTICAL-CANONICAL | `59491378a07010e0198d993baf9aaef84db202f1b180398766a0380085c0d1d4` |
| fw_tail_ctl | 0.7 min | IDENTICAL-CANONICAL | `140ae475017d66444762578f91c137f8f66492c42f907038571614c015d6ce67` |
| fw_tail_run | 7.9 min | IDENTICAL-CANONICAL | `8a4a44e0008159ffeb68e2aea278258d4044895e4e92c9b45100520a5280b0a2` |
| m_kstar | 1.2 min | IDENTICAL-CANONICAL | `2013b5a4f5260ef6102018b8a21c237cf56315193604d336f87bf2d8301cef10` |

Last line: `[reproduce] ALL CHECKS PASS`, exit 0.

## 3. Full tier

Command (in `run1/`): `python tools/reproduce.py --tier full --work work_full`. Started 2026-11-17 00:08:22, finished
03:40:04 (local time): **211.7 min**, peak 1.40 GB. Log: `run1_full.log`. The full tier never loads a shipped pickle:
all nine data files were regenerated from their generators in the work folder and compared byte for byte.

| Data file | Time | Status |
|---|---|---|
| gen_eqs, gen_fields, gen_lin_fields, gen_lin_fields2 | 1.9 min | IDENTICAL (sha256) |
| gen_tube13, gen_tube12 | 15.7 min | IDENTICAL |
| gen_cacheL, gen_fw4_bg, gen_fw4_bg3 | 9.7 min | IDENTICAL |

| Check (in addition to the nine quick checks, all IDENTICAL-CANONICAL again) | Time | Status |
|---|---|---|
| e5_lohner_thin, e5_lohner_box | 9.4 min | IDENTICAL-CANONICAL |
| e5_warmup (exact Friedmann data recovered) | 13.3 min | IDENTICAL-CANONICAL |
| e5_target (existence, Theorem 1) | 16.8 min | IDENTICAL-CANONICAL |
| e5_nullV0, e5_null2 (content nulls refuse) | 14.3 min | IDENTICAL-CANONICAL |
| m_C1 (gauge/winding controls), m_R3 (R1 = 3 count) | 23.5 min | IDENTICAL-CANONICAL (+ ledger jsonl) |
| m_count40 (the count, N = 1) | 71.9 min | IDENTICAL-CANONICAL (+ ledger jsonl) |
| fw_band1 [40, 55], fw_band2 [55, 2000] | 19.0 min | IDENTICAL-CANONICAL |

Last line: `[reproduce] ALL CHECKS PASS`, exit 0. (The G2 steps below ran alongside, in a separate folder.)

## 4. Blind G2 (added by the lead after the build seat stopped)

The build seat stopped on a network error after the quick tier; its full-tier process kept running and finished
(section 3). The lead then added the blind second kappa* computation to the repository (`code/b5g/`, 41 files;
`results/g2/` BAR + amendment; 24 ledger pins of seal 8a1c27c7fb / b62c6aecc9, all matching) and five full-tier steps.

- `run2/` (first attempt): every JSON output IDENTICAL; the two certificate logs DIFFERENT only in a run-time field
  (`"time_s"`) printed inside a record line. `tools/canon.py` was given that run-time rule (as it already had for
  stamps and JSON files); no certified character is ignored.
- `run3/` (fresh copy after the fix): `python tools/check_pins.py` ALL MATCH (108 pins, 100 patched files);
  `python tools/reproduce.py --only g2_certA g2_certB g2_refineA g2_refineB g2_result --work work_g2`:
  g2_certA 5.3 min, g2_certB 3.1 min, g2_refineA 0.8 min, g2_refineB 1.5 min, g2_result 0.0 min, all
  IDENTICAL-CANONICAL (logs and cert_A/B, sign_intervals_A, refine_A/B, result_g2 JSON). `ALL CHECKS PASS`, exit 0.
  Log: `run3_g2.log`.

The canon.py changes only add run-time rules; runs in sections 2-3 used the earlier, stricter version.

## Verdict

ALL CHECKS PASS: pins (108), quick tier (9), full tier (9 data files + 21 checks), G2 (5). Nothing published.
