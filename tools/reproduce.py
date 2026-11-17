#!/usr/bin/env python
"""reproduce.py -- re-run the certificates of this repository in a clean work folder and compare with results/.

Tiers (tools/reproduce.md has the table with commands and lab run times):
    python tools/reproduce.py --tier quick     pins, E5 tail controls, check_struct, B5 closure, Far Wall tail (controls
                                               and run), kappa*                         (about 20 minutes on one core)
    python tools/reproduce.py --tier full      everything: the quick checks, all eight data files regenerated from their
                                               generators and checked against their sha256, the E5 Lohner controls,
                                               warm-up, target and both content nulls, the M controls C1 and R1 = 3,
                                               the count to R1 = 40, both Far Wall bands         (about 5 hours)
    python tools/reproduce.py --only b5c kstar  named checks (the data they need are prepared first)
Options: --work DIR (default ./work; must not exist or be empty); --python "CMD" (default: this interpreter);
         --regen (quick tier / --only: regenerate the data files instead of using the shipped, sha256-checked ones);
         --list (print the checks and exit).

What it does
  1. tools/check_pins.py: every sealed file against its ledger pin (restored bytes for path-patched copies).
  2. Work folder = a mirror of the lab tree: code/ is copied to WORK/, and every path-patched copy is restored to its
     sealed bytes (PATCHES.md). The b5/ files get their lab modification times (tools/mtimes.json; close_b5.py C4
     compares them).
  3. Data files (pickles). Loading a pickle can execute code, so the full tier does not load the shipped ones at all:
     it regenerates all eight from their generators inside WORK (in dependency order) and checks each result against
     its pinned sha256. The quick tier copies the shipped pickles into WORK only after their sha256 matches the pin
     (tools/layout.json, tools/pins.json); pass --regen to regenerate there too.
  4. Every step runs as `python tools/run_sealed.py --root WORK <instrument> <args>` from WORK (the lab ran every
     step from its project root), with the LB_* / FW* environment flags cleared and only the sealed ones set. The
     launcher compiles every module with the lab-root literal pointed at WORK and refuses any open outside WORK.
  5. Output (stdout + stderr) goes to WORK/_out/<check>.log; launcher lines go to WORK/_out/<check>.launcher.log.
     Each output is compared with the shipped file by canonical content (tools/canon.py: elapsed times, peak memory,
     shell wrapper lines and the project folder removed; every other character must agree). Result files written by
     an instrument (close_b5_results.json, the count ledgers) are compared the same way; pickles byte for byte.
Last line: [reproduce] ALL CHECKS PASS (or the list of failures). Exit code 0 on success. A JSON report is written to
WORK/_out/reproduce_report.json.
"""
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOLS = os.path.join(ROOT, 'tools')
sys.path.insert(0, TOOLS)
import pathpatch  # noqa: E402
import canon  # noqa: E402

LADDER = {'FW2_NS_LADDER': '150,300,600'}

# data files: (name, lab path, generator argv (relative to WORK), minutes in the lab)
GENS = [
    ('gen_eqs', 'engine/eqs_gen.py', ['engine/derive.py'], 1),
    ('gen_fields', 'engine/cert/fields.pkl', ['engine/cert/fields.py'], 1),
    ('gen_lin_fields', 'engine/cert/lin_fields.pkl', ['engine/cert/lin_sym.py'], 1),
    ('gen_lin_fields2', 'engine/cert/lin_fields2.pkl', ['engine/cert/lin_mc.py'], 1),
    ('gen_tube13', 'engine/cert/lin_tube2_k1_3.pkl', ['engine/cert/lin_tube_v2.py', '1/3'], 9),
    ('gen_tube12', 'engine/cert/lin_tube2_k1_2.pkl', ['engine/cert/lin_tube_v2.py', '1/2'], 15),
    ('gen_cacheL', 'b4speed/L13_nsub2.pkl', ['b4speed/cacheL.py', 'build', '1/3', '2', '1e-20', 'b4speed/L13_nsub2.pkl'], 8),
    ('gen_fw4_bg', 'engine/cert/lin_fw4_bg.pkl', ['engine/cert/lin_fw4_bg.py'], 7),
    ('gen_fw4_bg3', 'engine/cert/lin_fw4_bg3.pkl', ['engine/cert/lin_fw4_bg3_make.py'], 0.1),
]

# checks: name -> (tier, argv relative to WORK, extra env, expected log (results/), [(produced file in WORK, shipped
#                  file, kind)], files that must be absent before the run, minutes in the lab, data files needed)
R = 'results/'
CHECKS = [
    ('e5_tails13', 'quick', ['engine/cert/test_tails_v32.py', '13', '1e-12'], {},
     R + 'e5/runs6/ctl_test_tails_v32.py_13_1e-12.log', [], [], 0.6, ['fields']),
    ('e5_tails12', 'quick', ['engine/cert/test_tails_v32.py', '12', '1e-12'], {},
     R + 'e5/runs6/ctl_test_tails_v32.py_12_1e-12.log', [], [], 0.8, ['fields']),
    ('e5_tails31', 'quick', ['engine/cert/test_tails_v31.py', 'v31'], {},
     R + 'e5/runs6/ctl_test_tails_v31.py_v31.log', [], [], 0.7, ['fields']),
    ('e5_struct', 'quick', ['engine/cert/check_struct.py'], {}, R + 'e5/runs6/ctl_check_struct.py.log', [], [], 0.1,
     ['fields']),
    ('b5c', 'quick', ['b5/close_b5.py', 'all'], {}, R + 'b5c/close_sealed.log',
     [('b5/close_b5_results.json', R + 'b5c/close_b5_results.json', 'json')], ['b5/close_b5_results.json'], 7.4, []),
    ('fw_tail_ctl', 'quick', ['engine/cert/fw_tail2.py', 'ctl'], {}, R + 'm/runs6/fw_tail2_ctl.log', [], [], 0.7,
     ['lin_fields', 'fw4_bg', 'fw4_bg3']),
    ('fw_tail_run', 'quick', ['engine/cert/fw_tail2.py', 'run'], {}, R + 'm/runs6/fw_tail2_KT2000.log', [], [], 8.1,
     ['lin_fields', 'fw4_bg', 'fw4_bg3']),
    ('m_kstar', 'quick', ['engine/cert/count_b4v3.py', 'kstar'], {}, R + 'm/runs7/kstar.log', [], [], 1.1,
     ['fields', 'lin_fields2', 'tube13', 'cacheL']),
    ('e5_lohner_thin', 'full', ['engine/cert/test_lohner2.py', 'thin'], {}, R + 'e5/runs6/ctl_test_lohner2.py_thin.log',
     [], [], 7.7, ['fields']),
    ('e5_lohner_box', 'full', ['engine/cert/test_lohner2.py', 'box'], {}, R + 'e5/runs6/ctl_test_lohner2.py_box.log',
     [], [], 7.8, ['fields']),
    ('e5_warmup', 'full', ['engine/cert/run_v32.py', 'warmup'], {}, R + 'e5/runs6/warmup.log', [], [], 21.7, ['fields']),
    ('e5_target', 'full', ['engine/cert/run_v32.py', 'target'], {}, R + 'e5/runs6/target.log', [], [], 24.2, ['fields']),
    ('e5_nullV0', 'full', ['engine/cert/run_v32.py', 'nullV0'], {}, R + 'e5/runs6/nullV0.log', [], [], 15.7, ['fields']),
    ('e5_null2', 'full', ['engine/cert/run_v32.py', 'null2'], {}, R + 'e5/runs6/null2_smoke.log', [], [], 7.5, ['fields']),
    ('m_C1', 'full', ['engine/cert/count_b4v3.py', 'smoke', '1/3', 'C1'], {}, R + 'm/runs7/ctl_C1.log', [], [], 14.2,
     ['fields', 'lin_fields2', 'tube13', 'tube12', 'cacheL']),
    ('m_R3', 'full', ['engine/cert/count_b4v3.py', 'count', '3', 'runs7/ctl_R3.jsonl'], {}, R + 'm/runs7/ctl_R3.log',
     [('runs7/ctl_R3.jsonl', R + 'm/runs7/ctl_R3.jsonl', 'jsonl')], ['runs7/ctl_R3.jsonl'], 9.6,
     ['fields', 'lin_fields2', 'tube13', 'tube12', 'cacheL']),
    ('m_count40', 'full', ['engine/cert/count_b4v3.py', 'count', '40', 'runs7/count_R40.jsonl'], {},
     R + 'm/runs7/count_R40.log', [('runs7/count_R40.jsonl', R + 'm/runs7/count_R40.jsonl', 'jsonl')],
     ['runs7/count_R40.jsonl'], 63.1, ['fields', 'lin_fields2', 'tube13', 'tube12', 'cacheL']),
    ('fw_band1', 'full', ['engine/cert/fw_cert3.py', 'run', '2.5', '40', '55', '160', '150', '1500', '2'], LADDER,
     R + 'm/runs6/fw_cert3_M2.5_Ka40_K55.log', [], [], 14.3, ['lin_fields', 'fw4_bg', 'fw4_bg3']),
    ('fw_band2', 'full', ['engine/cert/fw_cert3.py', 'run', '3.5', '55', '2000', '160', '150', '1500', '0'], {},
     R + 'm/runs6/fw_cert3_M3.5_Ka55_K2000.log', [], [], 15.4, ['lin_fields', 'fw4_bg', 'fw4_bg3']),
    # blind G2 (seal 8a1c27c7fb): own equations and gauge, reads only the JSON files its own derivation scripts wrote
    ('g2_certA', 'full', ['b5g/cert_g2.py', 'A'], {}, 'code/b5g/cert_A.log',
     [('b5g/cert_A.json', 'code/b5g/cert_A.json', 'json'),
      ('b5g/sign_intervals_A.json', 'code/b5g/sign_intervals_A.json', 'json')],
     ['b5g/cert_A.json', 'b5g/sign_intervals_A.json'], 4.6, []),
    ('g2_certB', 'full', ['b5g/cert_g2.py', 'B'], {}, 'code/b5g/cert_B.log',
     [('b5g/cert_B.json', 'code/b5g/cert_B.json', 'json')], ['b5g/cert_B.json'], 2.7, []),
    ('g2_refineA', 'full', ['b5g/refine_g2.py', 'A'], {}, 'code/b5g/refine_A.log',
     [('b5g/refine_A.json', 'code/b5g/refine_A.json', 'json')], ['b5g/refine_A.json'], 0.7, []),
    ('g2_refineB', 'full', ['b5g/refine_g2.py', 'B'], {}, 'code/b5g/refine_B.log',
     [('b5g/refine_B.json', 'code/b5g/refine_B.json', 'json')], ['b5g/refine_B.json'], 1.4, []),
    ('g2_result', 'full', ['b5g/make_result_g2.py'], {}, 'code/b5g/make_result.log',
     [('b5g/result_g2.json', 'code/b5g/result_g2.json', 'json')], ['b5g/result_g2.json'], 0.1, []),
]
# which data files each data file needs to be generated
GEN_DEPS = {'eqs': [], 'fields': [], 'lin_fields': [], 'lin_fields2': ['lin_fields'], 'tube13': ['fields'], 'tube12': ['fields'],
            'cacheL': ['fields', 'lin_fields2', 'tube13'], 'fw4_bg': ['fields'], 'fw4_bg3': ['fw4_bg']}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load_json(rel):
    return json.load(open(os.path.join(ROOT, rel), encoding='utf-8'))


def lab_root():
    t = load_json('tools/path_patch.json')
    return {e['id']: bytes.fromhex(e['original_hex']).decode() for e in t['table']}['R3']


def data_pins():
    """lab path of each data file -> (shipped repo path, expected sha256, source of the expected value)"""
    pins = {p['path']: p['sha256'] for p in load_json('tools/pins.json')['pins']}
    out = {}
    out['engine/eqs_gen.py'] = ('code/engine/eqs_gen.py', pins['code/engine/eqs_gen.py'], 'lab ledger pin')
    for g in load_json('tools/layout.json')['pickles']:
        repo = g['repo']
        if repo in pins:
            if pins[repo] != g['sha256']:
                sys.exit('[reproduce] tools/layout.json and tools/pins.json disagree on %s' % repo)
            out[g['lab']] = (repo, pins[repo], 'lab ledger pin')
        else:
            out[g['lab']] = (repo, g['sha256'], 'lab file as assembled; not pinned by the ledger')
    return out


def setup(work):
    os.makedirs(work, exist_ok=True)
    if os.listdir(work):
        sys.exit('[reproduce] work folder %s is not empty' % work)
    shutil.copytree(os.path.join(ROOT, 'code'), work, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('__pycache__'))
    patched = pathpatch.patched_files()
    n = 0
    for rel in patched:
        if rel.startswith('code/'):
            p = os.path.join(work, rel[5:])
            pathpatch.restore_file(p)
            if sha(p) != patched[rel]['original_sha256']:
                sys.exit('[reproduce] restore failed for %s' % rel)
            n += 1
    for d in ('runs6', 'runs7', '_out'):
        os.makedirs(os.path.join(work, d), exist_ok=True)
    mt = load_json('tools/mtimes.json')['mtimes']
    for rel, t in mt.items():
        p = os.path.join(work, rel)
        if os.path.exists(p):
            os.utime(p, (t, t))
    return n


def run(work, py, name, argv, env_extra, log):
    env = {k: v for k, v in os.environ.items() if not (k.startswith('LB_') or k.startswith('FW'))}
    env.update(env_extra)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['PYTHONIOENCODING'] = 'utf-8'
    out = os.path.join(work, '_out', name + '.raw')
    cmd = py + ['-u', os.path.join(TOOLS, 'run_sealed.py'), '--root', work, '--trace',
                os.path.join(work, '_out', name + '.trace'), os.path.join(work, argv[0])] + argv[1:]
    t0 = time.time()
    with open(out, 'wb') as fo:
        rc = subprocess.run(cmd, cwd=work, env=env, stdout=fo, stderr=subprocess.STDOUT).returncode
    raw = open(out, 'rb').read().decode('utf-8', errors='replace').splitlines(True)
    with open(os.path.join(work, '_out', name + '.log'), 'w', encoding='utf-8', newline='') as f:
        f.writelines(l for l in raw if not l.startswith('[run_sealed]'))
    with open(os.path.join(work, '_out', name + '.launcher.log'), 'w', encoding='utf-8', newline='') as f:
        f.writelines(l for l in raw if l.startswith('[run_sealed]'))
    os.remove(out)
    return rc, (time.time() - t0) / 60


def main(argv):
    work, py, tier, only, regen, lst = os.path.join(os.getcwd(), 'work'), [sys.executable], 'quick', None, False, False
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--work':
            work = os.path.abspath(argv[i + 1]); i += 2
        elif a == '--python':
            py = shlex.split(argv[i + 1]); i += 2
        elif a == '--tier':
            tier = argv[i + 1]; i += 2
        elif a == '--regen':
            regen = True; i += 1
        elif a == '--list':
            lst = True; i += 1
        elif a == '--only':
            only = []; i += 1
            while i < len(argv) and not argv[i].startswith('--'):
                only.append(argv[i]); i += 1
        else:
            sys.exit('[reproduce] unknown argument %s\n%s' % (a, __doc__))
    if tier not in ('quick', 'full'):
        sys.exit('[reproduce] --tier is quick or full')
    if lst:
        for g in GENS:
            print('%-15s full   %5.1f min  %s' % (g[0], g[3], ' '.join(g[2])))
        for c in CHECKS:
            print('%-15s %-6s %5.1f min  %s%s' % (c[0], c[1], c[7], ' '.join('%s=%s' % kv for kv in c[3].items()) + ' ' if c[3] else '', ' '.join(c[2])))
        return 0
    names = [c[0] for c in CHECKS] + [g[0] for g in GENS]
    for n in only or []:
        if n not in names:
            sys.exit('[reproduce] unknown check %s (known: %s)' % (n, ', '.join(names)))
    todo = [c for c in CHECKS if (c[0] in only if only is not None else (tier == 'full' or c[1] == 'quick'))]
    regen = regen or (tier == 'full' and only is None)

    def log(s):
        print('[reproduce] ' + s, flush=True)

    T0 = time.time()
    rc = subprocess.run(py + [os.path.join(TOOLS, 'check_pins.py')]).returncode
    if rc != 0:
        log('pins FAILED: a shipped file differs from its pinned sha256; stopping')
        return 1
    log('pins: every sealed file matches its pinned sha256')
    n = setup(work)
    log('work folder: %s (%d path-patched copies restored to their sealed bytes)' % (work, n))
    roots = [lab_root(), work]
    results = []
    dp = data_pins()
    need = set()
    gen_only = [g[0][4:] for g in GENS if only and g[0] in only]
    if gen_only:
        regen = True
    for c in todo + [(None,) * 8 + (gen_only,)]:
        for d in c[8]:
            stack = [d]
            while stack:
                x = stack.pop()
                if x not in need:
                    need.add(x)
                    stack += GEN_DEPS[x]
    if tier == 'full' and only is None:
        need = {g[0][4:] for g in GENS}
    for g in GENS:
        key = g[0][4:]
        if key not in need or (key == 'eqs' and not regen):
            continue
        repo, want, src = dp[g[1]]
        dst = os.path.join(work, g[1])
        if regen:
            log('regenerating %s with %s (about %s min)' % (g[1], ' '.join(g[2]), g[3]))
            if os.path.exists(dst):
                os.remove(dst)
            rc, mins = run(work, py, g[0], g[2], {}, log)
            got = sha(dst) if os.path.exists(dst) else None
            st = 'IDENTICAL' if rc == 0 and got == want else ('RUN FAILED (exit %d)' % rc if rc else 'DIFFERENT')
            results.append(dict(name=g[0], status=st, produced=got, expected=want, expected_source=src, minutes=mins,
                                output=g[1]))
            log('%s: %s in %.1f min (sha256 %s, expected %s, %s)' % (g[0], st, mins, got, want, src))
            if st != 'IDENTICAL':
                log('stopping: a data file could not be regenerated exactly')
                break
        else:
            src_p = os.path.join(ROOT, repo)
            if key == 'eqs':
                continue
            got = sha(src_p)
            if got != want:
                log('%s: shipped file sha256 %s differs from %s (%s); not loading it' % (repo, got, want, src))
                results.append(dict(name='data_' + key, status='DIFFERENT', produced=got, expected=want))
                break
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src_p, dst)
            log('data %s: shipped copy sha256 matches (%s); copied to the work folder' % (g[1], src))
    if all(r['status'].startswith('IDENTICAL') for r in results):
        for c in todo:
            name, _, args, env, expected, extra, absent, minutes, _ = c
            for a in absent:
                p = os.path.join(work, a)
                if os.path.exists(p):
                    os.remove(p)
            log('running %s: %s%s (about %s min in the lab)' % (
                name, ' '.join('%s=%s ' % kv for kv in env.items()), ' '.join(args), minutes))
            rc, mins = run(work, py, name, args, env, log)
            r = dict(name=name, minutes=mins, rc=rc, output='_out/%s.log' % name, shipped=expected)
            if rc != 0:
                r['status'] = 'RUN FAILED (exit %d)' % rc
            else:
                st, hp, he, diffs = canon.compare('log', open(os.path.join(work, '_out', name + '.log'), 'rb').read(),
                                                  pathpatch.restore_bytes(open(os.path.join(ROOT, expected), 'rb').read()),
                                                  roots)
                r.update(status=st, produced=hp, expected=he, differences=diffs, files=[])
                for prod, ship, kind in extra:
                    pp = os.path.join(work, prod)
                    if not os.path.exists(pp):
                        s2 = ('MISSING', None, None, [])
                    else:
                        s2 = canon.compare(kind, open(pp, 'rb').read(),
                                           pathpatch.restore_bytes(open(os.path.join(ROOT, ship), 'rb').read()), roots)
                    r['files'].append(dict(produced=prod, shipped=ship, status=s2[0], sha_produced=s2[1],
                                           sha_expected=s2[2], differences=s2[3]))
                    if not s2[0].startswith('IDENTICAL'):
                        r['status'] = 'DIFFERENT (%s)' % prod
            results.append(r)
            log('%s: %s in %.1f min' % (name, r['status'], mins))
            if 'produced' in r:
                log('  canonical sha256 produced %s expected %s' % (r['produced'], r['expected']))
            for d in r.get('differences') or []:
                log('  ' + d)
            for f in r.get('files', []):
                log('  %s: %s (canonical sha256 %s, expected %s)' % (f['produced'], f['status'], f['sha_produced'],
                                                                   f['sha_expected']))
                for d in f['differences']:
                    log('    ' + d)
    with open(os.path.join(work, '_out', 'reproduce_report.json'), 'w', encoding='utf-8') as f:
        json.dump(dict(tier=tier, only=only, regen=regen, minutes=(time.time() - T0) / 60, results=results), f, indent=1)
    bad = [r['name'] for r in results if not r['status'].startswith('IDENTICAL')]
    ran = {r['name'] for r in results}
    missing = [c[0] for c in todo if c[0] not in ran]
    log('total %.1f min' % ((time.time() - T0) / 60))
    log('ALL CHECKS PASS' if not bad and not missing else 'FAILED: ' + ', '.join(bad + ['%s (not run)' % m for m in missing]))
    return 0 if not bad and not missing else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
