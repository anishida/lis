#!/usr/bin/env python3
"""Audit native LIS ILUTP-full against ORIGINAL Matrix Market A,b snapshots.

No SuperLU dependency, no manufactured RHS and no modification of A or b.
Bench executable is the LIS ./test/benchilutp built from this branch.

Exit status: 0 if each full-fill case passes its original-system residual,
1 for rejected solver results/timeouts/limits, 2 for invalid CLI/input.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def integer_nonnegative(value: str) -> int:
    if not value.isascii() or not value.isdecimal():
        raise argparse.ArgumentTypeError('expected a nonnegative decimal integer')
    try:
        n = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError('invalid decimal integer') from exc
    # The benchmark itself validates against the build's LIS_INT width.
    # Do not reject valid 64-bit LIS_INT limits in this Python front end.
    if n > 2 ** 63 - 1:
        raise argparse.ArgumentTypeError('factor limit exceeds 64-bit LIS_INT range')
    return n


def positive_number(value: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError('expected positive finite number')
    return number


def nonnegative_number(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError('expected nonnegative finite number')
    return number


def gnu_time() -> str | None:
    path = Path('/usr/bin/time')
    if not path.is_file():
        return None
    try:
        test = subprocess.run([str(path), '-f', '%M', '-o', '/dev/null', 'true'],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                              timeout=2, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return str(path) if test.returncode == 0 else None


def memory_limit(mib: int):
    """Bound address space (not RSS), only on POSIX with RLIMIT_AS."""
    import resource
    byte_limit = mib * 1024 * 1024

    def set_limit() -> None:
        resource.setrlimit(resource.RLIMIT_AS, (byte_limit, byte_limit))
    return set_limit



def run_process_tree(command: list[str], timeout: float, preexec,
                     env: dict[str, str]) -> subprocess.CompletedProcess:
    """Reap the benchmark AND /usr/bin/time's descendant process tree.

    subprocess.run(timeout=...) only kills its immediate child.  When that
    child is GNU time, the underlying solver could otherwise outlive the
    audit's timeout and keep consuming CPU or memory.
    """
    proc = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, preexec_fn=preexec, env=env,
        start_new_session=(os.name == 'posix'))
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        if os.name == 'posix':
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            proc.kill()
        proc.communicate()  # Always reap the direct child.
        raise
    return subprocess.CompletedProcess(command, proc.returncode, stdout, stderr)



def parse_time_file(path: Path) -> dict:
    metrics = {}
    if not path.exists():
        return metrics
    for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
        if line.startswith('MAX_RSS_KIB='):
            try:
                metrics['peak_rss_kib'] = int(line.split('=', 1)[1])
            except ValueError:
                pass
        if line.startswith('WALL_SECONDS='):
            try:
                metrics['gnu_wall_seconds'] = float(line.split('=', 1)[1])
            except ValueError:
                pass
    return metrics


def parse_full_row(output: str) -> tuple[dict | None, str | None]:
    # The native C benchmark CSV is kept in a separate file.  This method
    # validates its schema rather than interpreting human-facing LIS logs.
    try:
        rows = list(csv.DictReader(io.StringIO(output)))
    except csv.Error as exc:
        return None, str(exc)
    if not rows:
        return None, 'empty or invalid CSV'
    full_rows = [row for row in rows if row.get('variant') == 'ilutp_full']
    if len(full_rows) != 1:
        return None, f'expected one ilutp_full row, got {len(full_rows)}'
    return full_rows[0], None


def run_one(args: argparse.Namespace, name: str, matrix: Path, rhs: Path,
            log_dir: Path, timing_tool: str | None) -> dict:
    result: dict = {'case': name, 'matrix': str(matrix), 'rhs': str(rhs),
                    'accepted': False, 'reasons': [], 'metrics': {}}
    hashes = {}
    for which, path in (('matrix', matrix), ('rhs', rhs)):
        if not path.is_file():
            result['reasons'].append(f'missing {which}: {path}')
        else:
            hashes[which] = sha256_file(path)
    result['sha256_before'] = hashes
    if len(hashes) != 2:
        return result
    safe_name = re.sub(r'[^a-zA-Z0-9._-]', '_', name)[:80]
    csv_path = log_dir / f'{safe_name}.csv'
    log_path = log_dir / f'{safe_name}.log'
    # Never reuse a CSV from a prior run after an early benchmark failure.
    csv_path.unlink(missing_ok=True)
    result['benchmark_csv'] = str(csv_path)
    result['benchmark_log'] = str(log_path)
    command = [str(args.bench), str(matrix), str(rhs),
               '--full', '1', '--only-full', '0' if args.compare else '1',
               '--full-limit', str(args.full_limit),
               '--full-tau', str(args.full_tau),
               '--maxiter', str(args.maxiter), '--restart', str(args.restart),
               '--rtol', str(args.rtol), '--drop', str(args.drop),
               '--rate', str(args.rate), '--rcm', str(args.rcm),
               '--row-scale', str(args.row_scale), '--csv', str(csv_path)]
    result['command'] = command
    with tempfile.TemporaryDirectory(prefix='lis-native-timing-') as d:
        time_file = Path(d) / 'time.txt'
        launch = ([timing_tool, '-f',
                   'MAX_RSS_KIB=%M\nWALL_SECONDS=%e', '-o', str(time_file)]
                  + command) if timing_tool else command
        preexec = memory_limit(args.max_vmem_mib) if args.max_vmem_mib else None
        start = time.monotonic()
        try:
            completed = run_process_tree(
                launch, timeout=args.timeout, preexec=preexec,
                env={**os.environ, 'OMP_NUM_THREADS': '1',
                     'OMP_DYNAMIC': 'FALSE'})
            result['exit_code'] = completed.returncode
            log_path.write_text(completed.stdout + '\n--- STDERR ---\n' + completed.stderr,
                                encoding='utf-8')
        except subprocess.TimeoutExpired as exc:
            result['exit_code'] = None
            result['reasons'].append(f'timeout after {args.timeout}s')
            with log_path.open('w', encoding='utf-8') as f:
                f.write(f'TIMEOUT after {args.timeout}s\n{exc}')
        except OSError as exc:
            result['exit_code'] = None
            result['reasons'].append(f'could not run benchmark: {exc}')
            log_path.write_text(str(exc), encoding='utf-8')
        result['metrics']['wall_seconds'] = time.monotonic() - start
        result['metrics'].update(parse_time_file(time_file))
    if result['exit_code'] != 0:
        result['reasons'].append(f'benchmark exited with {result["exit_code"]}')
    if not csv_path.is_file():
        result['reasons'].append('benchmark did not write CSV')
    else:
        row, error = parse_full_row(csv_path.read_text(encoding='utf-8'))
        if error:
            result['reasons'].append(error)
        else:
            result['solver_row'] = row
            for key in ('iterations', 'nnz_L', 'nnz_U', 'moved_columns',
                        'api_error', 'solver_status'):
                try:
                    result['metrics'][key] = int(row[key])
                except (KeyError, ValueError, TypeError):
                    result['reasons'].append(f'invalid integer: {key}')
            for key in ('setup_cpu_s', 'solve_cpu_s', 'physical_relres'):
                try:
                    value = float(row[key])
                    if not math.isfinite(value) or (key == 'physical_relres' and value < 0):
                        raise ValueError(key)
                    result['metrics'][key] = value
                except (KeyError, ValueError, TypeError):
                    result['reasons'].append(f'invalid finite number: {key}')
            if row.get('status') != 'converged':
                result['reasons'].append('native full-fill did not converge')
            if result['metrics'].get('api_error') != 0:
                result['reasons'].append('LIS API returned an error')
            if result['metrics'].get('solver_status') != 0:
                result['reasons'].append('LIS solver reported a failure')
            if result['metrics'].get('physical_relres', float('inf')) > args.accept_residual:
                result['reasons'].append('original-system physical residual exceeds acceptance')
            if result['metrics'].get('nnz_L', 0) + result['metrics'].get('nnz_U', 0) <= 0:
                result['reasons'].append('empty factors')
    if args.max_rss_mib is not None:
        rss = result['metrics'].get('peak_rss_kib')
        if rss is None:
            result['reasons'].append('peak RSS unavailable for requested RSS limit')
        elif rss > args.max_rss_mib * 1024:
            result['reasons'].append('peak RSS exceeds limit')
    result['sha256_after'] = {
        which: sha256_file(path) for which, path in (('matrix', matrix), ('rhs', rhs))
        if path.is_file()}
    if result['sha256_after'] != result['sha256_before']:
        result['reasons'].append('original A/b files changed during benchmark')
    result['accepted'] = not result['reasons']
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bench', type=Path, default=Path('./test/benchilutp'))
    parser.add_argument('--case', nargs=3, action='append', required=True,
                        metavar=('LABEL', 'A.mtx', 'b.mtx'))
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--csv', type=Path, default=None,
                        help='optional one-line-per-case summary')
    parser.add_argument('--compare', action='store_true',
                        help='also record three incomplete ILUT/P controls')
    parser.add_argument('--full-limit', type=integer_nonnegative, default=20000000)
    parser.add_argument('--full-tau', type=nonnegative_number, default=1.0,
                        help='native complete-fill ILUTP pivot threshold [0,1]')
    parser.add_argument('--rtol', type=positive_number, default=1e-6)
    parser.add_argument('--accept-residual', type=positive_number, default=1e-5)
    parser.add_argument('--maxiter', type=int, default=300)
    parser.add_argument('--restart', type=int, default=80)
    parser.add_argument('--drop', type=nonnegative_number, default=0.005)
    parser.add_argument('--rate', type=nonnegative_number, default=12)
    parser.add_argument('--rcm', type=int, choices=(0, 1), default=1)
    parser.add_argument('--row-scale', type=int, choices=(0, 1), default=1)
    parser.add_argument('--timeout', type=positive_number, default=600,
                        help='max seconds per benchmark case')
    parser.add_argument('--max-rss-mib', type=positive_number, default=None,
                        help='reject if measured max RSS is above this threshold')
    parser.add_argument('--max-vmem-mib', type=int, default=None,
                        help='POSIX RLIMIT_AS hard virtual-memory cap per process')
    args = parser.parse_args(argv)
    if args.full_tau > 1.0:
        parser.error('--full-tau must be <= 1')
    if args.maxiter <= 0 or args.restart <= 0:
        parser.error('--maxiter and --restart must be positive')
    if args.max_vmem_mib is not None and args.max_vmem_mib <= 0:
        parser.error('--max-vmem-mib must be positive')
    if args.max_vmem_mib is not None and os.name != 'posix':
        parser.error('--max-vmem-mib requires POSIX RLIMIT_AS')
    labels = [triple[0] for triple in args.case]
    if len(set(labels)) != len(labels):
        parser.error('case labels must be unique')
    safe_labels = [re.sub(r'[^a-zA-Z0-9._-]', '_', x)[:80] for x in labels]
    if len(set(safe_labels)) != len(safe_labels):
        parser.error('case labels must have distinct safe filenames')
    args.bench = args.bench.resolve()
    if not args.bench.is_file() or not os.access(args.bench, os.X_OK):
        parser.error(f'not an executable LIS benchmark: {args.bench}')
    args.report = args.report.resolve()
    if args.csv is not None:
        args.csv = args.csv.resolve()
    logs = args.report.parent / (args.report.stem + '_logs')
    logs.mkdir(parents=True, exist_ok=True)
    cases = []
    timing_tool = gnu_time()
    for name, a_path, b_path in args.case:
        result = run_one(args, name, Path(a_path).resolve(), Path(b_path).resolve(),
                         logs, timing_tool)
        cases.append(result)
        metric = result['metrics']
        print(f'{name}: {"ACCEPTED" if result["accepted"] else "REJECTED"} '
              f'relres={metric.get("physical_relres", "n/a")} '
              f'fill={metric.get("nnz_L", 0)+metric.get("nnz_U", 0)} '
              f'rss_kib={metric.get("peak_rss_kib", "n/a")} '
              f'reasons={result["reasons"]}', flush=True)
    report = {'schema': 'native-lis-contact-v1',
              'implementation': 'native LIS ILUTP full fill; no SuperLU',
              'benchmark': str(args.bench),
              'requirements': {'rtol': args.rtol,
                               'max_physical_relres': args.accept_residual,
                               'full_factor_nnz_limit': args.full_limit,
                               'max_rss_mib': args.max_rss_mib,
                               'max_vmem_mib': args.max_vmem_mib,
                               'timeout_seconds_per_case': args.timeout},
              'cases': cases, 'accepted': bool(cases) and all(x['accepted'] for x in cases)}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n',
                           encoding='utf-8')
    if args.csv is not None:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=(
                'case', 'accepted', 'reason', 'physical_relres', 'iterations',
                'nnz_L', 'nnz_U', 'moved_columns', 'setup_cpu_s',
                'solve_cpu_s', 'wall_seconds', 'peak_rss_kib',
                'matrix_sha256', 'rhs_sha256'))
            writer.writeheader()
            for result in cases:
                metrics = result['metrics']
                writer.writerow({'case': result['case'],
                                 'accepted': result['accepted'],
                                 'reason': '; '.join(result['reasons']),
                                 **{key: metrics.get(key) for key in (
                                     'physical_relres', 'iterations', 'nnz_L',
                                     'nnz_U', 'moved_columns', 'setup_cpu_s',
                                     'solve_cpu_s', 'wall_seconds', 'peak_rss_kib')},
                                 'matrix_sha256': result['sha256_before'].get('matrix'),
                                 'rhs_sha256': result['sha256_before'].get('rhs')})
    return 0 if report['accepted'] else 1


if __name__ == '__main__':
    sys.exit(main())
