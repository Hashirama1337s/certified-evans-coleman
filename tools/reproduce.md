# Reproduce: re-run the certificates

Requirements: Python 3.11 (tested with 3.11.9), python-flint 0.9.0, sympy 1.14.0, mpmath 1.3.0, numpy 2.4.4.
Every step is a single process; the largest peak in the lab was about 1.1 GB (the mode count). The lab ran every step
under an 8 GB memory cap. Lean 4 + Mathlib only for the optional Lean checks (see the end).

## Commands

From the repository root:

```
python tools/check_pins.py                   # every sealed file against its ledger pin; every file against MANIFEST.sha256
python tools/reproduce.py --tier quick       # about 20 minutes on one core
python tools/reproduce.py --tier full        # about 5 hours on one core
```

Expected last line: `[reproduce] ALL CHECKS PASS`. The runner works in `./work` (or `--work DIR`, which must be
empty or absent); `--python "CMD"` picks another interpreter; `--list` prints every step; `--only NAME ...` runs named
steps (the data files they need are prepared first); `--regen` makes the quick tier regenerate its data files instead
of using the shipped ones. A JSON report is written to `work/_out/reproduce_report.json`; every step's output is in
`work/_out/<step>.log`, and the launcher's substitution lines in `work/_out/<step>.launcher.log`.

## Steps

Times are the lab's run times on one core (the clean-room times are in the clean-room report). "Canonical" = the
output compared with the shipped log by `tools/canon.py`: elapsed-time values, peak memory, the shell's wrapper lines
and the project folder are removed; every other character (each enclosure, count, parity flag, sha256 provenance line
and verdict line) must agree. Data files are compared byte for byte.

### Quick tier

| Step | Claim | Command (from the work folder) | Lab time | Compared with |
|---|---|---|---|---|
| pins | all | `tools/check_pins.py` | seconds | `tools/pins.json` (84 ledger pins), `tools/path_patch.json` |
| e5_tails13 | E5 control | `engine/cert/test_tails_v32.py 13 1e-12` | 0.6 min | `results/e5/runs6/ctl_test_tails_v32.py_13_1e-12.log` |
| e5_tails12 | E5 control | `engine/cert/test_tails_v32.py 12 1e-12` | 0.8 min | `results/e5/runs6/ctl_test_tails_v32.py_12_1e-12.log` |
| e5_tails31 | E5 control | `engine/cert/test_tails_v31.py v31` | 0.7 min | `results/e5/runs6/ctl_test_tails_v31.py_v31.log` |
| e5_struct | E5 control | `engine/cert/check_struct.py` | seconds | `results/e5/runs6/ctl_check_struct.py.log` |
| b5c | B5C (all of C1-C5, NEG) | `b5/close_b5.py all` | 7.4 min | `results/b5c/close_sealed.log`; `b5/close_b5_results.json` vs `results/b5c/close_b5_results.json` (canonical JSON without the `seconds` fields) |
| fw_tail_ctl | M, Far Wall tail controls C1-C3 | `engine/cert/fw_tail2.py ctl` | 0.7 min | `results/m/runs6/fw_tail2_ctl.log` |
| fw_tail_run | M, Far Wall tail, every K >= 2000 | `engine/cert/fw_tail2.py run` | 8.1 min | `results/m/runs6/fw_tail2_KT2000.log` |
| m_kstar | M, kappa* and gamma | `engine/cert/count_b4v3.py kstar` | 1.1 min | `results/m/runs7/kstar.log` |

In the quick tier the data files are the shipped ones (`results/pickles/`), copied into the work folder only after
their sha256 matches the ledger pin (see "Data files" below).

### Full tier (adds)

| Step | Claim | Command (from the work folder) | Lab time | Compared with |
|---|---|---|---|---|
| gen_eqs | engine | `engine/derive.py` (writes `engine/eqs_gen.py`) | under a minute | ledger pin of `engine/eqs_gen.py` (byte for byte) |
| gen_fields ... gen_fw4_bg3 | data | the eight generators in the table below | about 40 min | ledger pin / lab sha256 (byte for byte) |
| e5_lohner_thin | E5 control | `engine/cert/test_lohner2.py thin` | 7.7 min | `results/e5/runs6/ctl_test_lohner2.py_thin.log` |
| e5_lohner_box | E5 control | `engine/cert/test_lohner2.py box` | 7.8 min | `results/e5/runs6/ctl_test_lohner2.py_box.log` |
| e5_warmup | E5 PE4 (k = 1/2 positive control) | `engine/cert/run_v32.py warmup` | 21.7 min | `results/e5/runs6/warmup.log` |
| e5_target | E5 PE1-PE3 (the existence proof) | `engine/cert/run_v32.py target` | 24.2 min | `results/e5/runs6/target.log` |
| e5_nullV0 | E5 content nulls | `engine/cert/run_v32.py nullV0` | 15.7 min | `results/e5/runs6/nullV0.log` |
| e5_null2 | E5 smoke check | `engine/cert/run_v32.py null2` | 7.5 min | `results/e5/runs6/null2_smoke.log` |
| m_C1 | M controls C1a-c | `engine/cert/count_b4v3.py smoke 1/3 C1` | 14.2 min | `results/m/runs7/ctl_C1.log` |
| m_R3 | M control R1 = 3 | `engine/cert/count_b4v3.py count 3 runs7/ctl_R3.jsonl` | 9.6 min | `results/m/runs7/ctl_R3.log`; the disc ledger vs `results/m/runs7/ctl_R3.jsonl` (without `sec`) |
| m_count40 | M, N(R1 = 40) = 1 | `engine/cert/count_b4v3.py count 40 runs7/count_R40.jsonl` | 63.1 min | `results/m/runs7/count_R40.log`; the disc ledger vs `results/m/runs7/count_R40.jsonl` |
| fw_band1 | M, Far Wall 40 <= abs(kappa) <= 55 | `FW2_NS_LADDER=150,300,600 engine/cert/fw_cert3.py run 2.5 40 55 160 150 1500 2` | 14.3 min | `results/m/runs6/fw_cert3_M2.5_Ka40_K55.log` |
| fw_band2 | M, Far Wall 55 <= abs(kappa) <= 2000 | `engine/cert/fw_cert3.py run 3.5 55 2000 160 150 1500 0` | 15.4 min | `results/m/runs6/fw_cert3_M3.5_Ka55_K2000.log` |
| g2_certA | G2 (blind kappa*), configuration A | `b5g/cert_g2.py A` | 4.6 min | `code/b5g/cert_A.log`, `cert_A.json`, `sign_intervals_A.json` |
| g2_certB | G2, configuration B (every setting different) | `b5g/cert_g2.py B` | 2.7 min | `code/b5g/cert_B.log`, `cert_B.json` |
| g2_refineA | G2, certified bisection of kappa* (A) | `b5g/refine_g2.py A` | 0.7 min | `code/b5g/refine_A.log`, `refine_A.json` |
| g2_refineB | G2, certified bisection of kappa* (B) | `b5g/refine_g2.py B` | 1.4 min | `code/b5g/refine_B.log`, `refine_B.json` |
| g2_result | G2, assemble the verdict file | `b5g/make_result_g2.py` | seconds | `code/b5g/make_result.log`, `result_g2.json` |

The G2 steps read only JSON files written by the blind seat's own symbolic scripts (`b5g/derive_g2.py`,
`centre_sym.py`, `sonic_sym_g2.py`, `bgjac_g2.py`; see `code/b5g/RUN.md`), which are shipped and pinned. Its sealed
test is `results/g2/BAR-LB-G2.md` with the amendment `AMEND-LB-G2-1.md`, sealed before the run.

`results/m/runs7/verdict_M.json` was assembled by the lab from these logs (no program writes it): every enclosure in
it appears verbatim in `ctl_C1.log`, `ctl_R3.log`, `count_R40.log` or `kstar.log`, and its two bands
(`2.8106 +- 0.0005`, `0.3558 +- 0.0001`) are the sealed bands of `results/m/BAR-LB-M.md`.

## Data files (pickles)

Loading a pickle can execute code. The full tier therefore never loads a shipped pickle: it regenerates all eight in
the work folder, in this order, and checks each against its sha256 before the next step reads it.

| File (work folder) | Generator | Reads | Lab time | Expected sha256 |
|---|---|---|---|---|
| `engine/cert/fields.pkl` | `engine/cert/fields.py` | `engine/eqs_gen.py` | seconds | `736818c0...` (ledger pin, E5) |
| `engine/cert/lin_fields.pkl` | `engine/cert/lin_sym.py` | `engine/eqs_gen.py` | 2 min | `6e645282...` (lab file; not pinned) |
| `engine/cert/lin_fields2.pkl` | `engine/cert/lin_mc.py` | `lin_fields.pkl` | seconds | `dce8eba8...` (ledger pin, M) |
| `engine/cert/lin_tube2_k1_3.pkl` | `engine/cert/lin_tube_v2.py 1/3` | `fields.pkl` (E5 boxes are literals in the file) | 8-9 min | `f6bff9e2...` (ledger pin, M) |
| `engine/cert/lin_tube2_k1_2.pkl` | `engine/cert/lin_tube_v2.py 1/2` | `fields.pkl` | 12-15 min | `17dadbfc...` (ledger pin, M) |
| `b4speed/L13_nsub2.pkl` | `b4speed/cacheL.py build 1/3 2 1e-20 b4speed/L13_nsub2.pkl` | `lin_fields2.pkl`, `lin_tube2_k1_3.pkl` | 8 min | `c239cbe5...` (ledger pin, M; also checked by `count_b4v3.py` itself) |
| `engine/cert/lin_fw4_bg.pkl` | `engine/cert/lin_fw4_bg.py` | `fields.pkl` | 7 min | `a5d3f462...` (lab file; not pinned) |
| `engine/cert/lin_fw4_bg3.pkl` | `engine/cert/lin_fw4_bg3_make.py` | `lin_fw4_bg.pkl` | seconds | `53008a7d...` (ledger pin, M) |

The quick tier uses the shipped copies (`results/pickles/`) after the same sha256 check; if you do not want to load
shipped pickles at all, add `--regen` (about 40 extra minutes). The generation logs of the lab are shipped as
`results/pickles/gen_*.log`.

## What the runner does (and why)

- **Work folder = the lab tree.** The instruments find each other and their data by paths under one project folder
  (`engine/`, `engine/cert/`, `b4speed/`, `b5/`), and the lab ran every step from that folder. `code/` is that tree;
  the runner copies it to `work/` and runs every step with `work/` as the current folder.
- **Sealed bytes.** The shipped copies are path-patched (PATCHES.md); the runner restores the originals in `work/`
  first, so every file has exactly the bytes the ledger pinned.
- **Lab folder literals.** The instruments name the lab folder in string literals (`sys.path.insert`, data file
  names, a `runpy` target, a default cache location). `tools/run_sealed.py` compiles every module loaded from `work/`
  with only that literal pointed at `work/`; the files on disk stay unchanged, so every self-hash (the engine sha256s
  in the count ledger header, fw_tail2.py's provenance line, close_b5.py's cache key) is the sealed value. One recorded
  config file names a data file by its lab path (`b5/cfg_ec.json`); the launcher redirects such an open() to `work/`.
  Every substitution is printed (`work/_out/<step>.launcher.log`), and an audit hook refuses any open() outside
  `work/` under the lab folder name, so a missed literal fails instead of reading elsewhere.
- **File times.** `close_b5.py` (C4) checks that the cached proof objects are newer than the files that produced
  them. Version control does not keep file times, so the runner gives the `b5/` files their lab modification times
  (`tools/mtimes.json`). The recorded outcome includes the FRW cache predating the last `prove.py` edit, as in the lab.
- **Environment.** The `LB_*` and `FW*` flags are cleared before each step (`run_v32.py` refuses to start with any set);
  only `FW2_NS_LADDER=150,300,600` is set for the first Far Wall band, as in the lab.

## Running one step by hand

```
python tools/reproduce.py --only e5_struct --work run     # sets up run/ and runs one step
cd run
python ../tools/run_sealed.py --root . engine/cert/check_struct.py
```

## Lean (optional; not run by tools/reproduce.py)

`code/lean_fw/` holds the Lean 4 files of the Far Wall Collatz step and their generators (`lean_all.py` extracts the
rational rows of the 133 bulk cells from `fw_cert3.py`, `lean_tail.py` those of the tail from `fw_tail2.py`;
`LEAN-NOTE.md` says exactly what is and is not kernel-checked). To check them, put the four `.lean` files into a Lake
project that depends on Mathlib (the lab used toolchain v4.36.0-rc2) and run `lake env lean <file>` on each:
`fw_collatz_all.lean`, `fw_collatz_tail.lean` and `fw_collatz_logic.lean` must compile with no `sorry` and only the
axioms `propext`, `Classical.choice`, `Quot.sound`; `fw_collatz_logic_control.lean` is a teeth control and must FAIL.
The lab's check logs are in `results/lean/`.

## Figures (optional; not run by tools/reproduce.py)

The three README figures are drawn from shipped files only: `python figures/make_figures.py .` (Python 3.11,
matplotlib 3.11; seconds; reads only shipped files). It rewrites `figures/*.svg`; two runs give byte-identical SVGs.
