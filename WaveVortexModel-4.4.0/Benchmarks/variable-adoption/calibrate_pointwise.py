#!/usr/bin/env python3
"""Bounded integrated pointwise calibration: serial and the prior 4/8-worker finalists."""
import argparse
import json
from pathlib import Path
import os
import platform

import numpy as np
from run_derivative_qualification import compare, digest, save, source_receipt
from run_compact_qualification import run_child


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('manifest', 'source', 'build', 'receipt', 'output'):
        p.add_argument(name, type=Path)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    fixture = json.loads(args.manifest.read_text())
    receipt = json.loads(args.receipt.read_text())
    source = source_receipt(args.source.resolve())
    worker = args.build.resolve() / 'wv-variable-adoption'
    worker_hash = digest(worker)
    if (source['status'] or source['commit'] != receipt['source']['commit']
            or source['files'] != receipt['source']['files'] or worker_hash != receipt['worker']['sha256']):
        raise ValueError('Frozen build receipt mismatch')
    cache = (args.build / 'CMakeCache.txt').read_text()
    configured = next(line.split('=',1)[1] for line in cache.splitlines() if line.startswith('WVM_SOURCE:PATH='))
    if Path(configured).resolve() != args.source.resolve() or not receipt.get('harness'):
        raise ValueError('Build configuration or harness identity mismatch')
    provider = receipt['verifiedRuntimeProvider']
    counts = (1, 4, 8)
    save(args.output / 'protocol.json', dict(schema='wvm-variable-pointwise-calibration-v1', fixture=fixture,
         manifestSHA256=digest(args.manifest), source=source, buildReceipt=receipt,
         runnerSHA256=digest(__file__), host=platform.platform(),
         environment={k:v for k,v in os.environ.items() if k.startswith(('VECLIB_','OMP_','MKL_','OPENBLAS_','DYLD_'))}, pointwiseWorkers=counts, horizontalWorkers=12,
         fftInternalWorkers=1, generalVerticalOuterWorkers=1, blocks=3, warmups=1, samples=3,
         rule='One count across all declared profiles; smallest count within 1% of fastest geometric complete-flux time',
         rationale='Serial control plus spectral-kernel-benchmarks #23 Lyra interleaved/compact provisional pointwise finalists 4 and 8.',
         retention='First payload per profile/count and every failure; later verified duplicates hashed then removed'))
    runs, retained = [], {}
    for case in fixture['cases']:
        for path, expected in (('sourcePath', 'sourceSHA256'), ('matlabFluxPath', 'matlabFluxSHA256')):
            if digest(case[path]) != case[expected]:
                raise ValueError('Fixture identity mismatch')
        for block in range(3):
            order = counts[block:] + counts[:block]
            for count in order:
                prefix = args.output / f"{case['id']}-{block}-p{count}"
                payload = prefix.with_suffix('.bin')
                command = [str(worker), case['sourcePath'], 'compact', '12', '1', '3', str(payload), str(count)]
                row = dict(profile=case['id'], block=block, pointwiseWorkers=count, command=command)
                runs.append(row); save(args.output / 'runs.json', runs)
                row.update(run_child(command, prefix)); save(args.output / 'runs.json', runs)
                if row['returncode']:
                    raise RuntimeError('Worker failed; evidence retained')
                report = json.loads(prefix.with_suffix('.stdout').read_text()); row['report'] = report
                family={'stratified-qg':'WVTransformStratifiedQG','hydrostatic':'WVTransformHydrostatic','boussinesq':'WVTransformBoussinesq'}[case['family']]
                save(args.output / 'runs.json', runs)
                if (report['schema'] != 'wvm-variable-screening-v1' or report['family'] != family
                        or report['warmups'] != 1 or not report['streamedNonlinear']
                        or report['horizontalSchedule'] != 'fftw-streaming-pruned-tile16'
                        or any(report['provider'][key] != provider[key] for key in ('version','baseLibrary','threadLibrary','fftThreads'))
                        or report['grid'] != case['grid'] or report['selection'] != 'compact'
                        or report['pointwiseWorkers'] != count or not report['compactSplitViews']
                        or report['matrixBackendIdentifier'] != 'accelerate-gemm-v1'
                        or report['horizontalWorkers'] != 12 or report['provider']['fftThreads'] != 1
                        or len(report['samplesSeconds']) != 3
                        or not all(np.isfinite(x) and x > 0 for x in report['samplesSeconds'])):
                    raise ValueError('Worker metadata mismatch')
                for key in ('baseLibrary', 'threadLibrary'):
                    report['provider'][key+'SHA256'] = digest(report['provider'][key])
                    if report['provider'][key+'SHA256'] != provider[key+'SHA256']:
                        raise ValueError('Provider library identity changed')
                row['matlabComparison'] = compare(payload, Path(case['matlabFluxPath']), True)
                row['payloadSHA256'] = digest(payload); save(args.output / 'runs.json', runs)
                if not row['matlabComparison']['passed']:
                    raise ValueError('Scientific comparison failed')
                key = (case['id'], count)
                row['payloadRetained'] = key not in retained
                if key in retained: payload.unlink()
                else: retained[key] = payload
                save(args.output / 'runs.json', runs)
                print(f"{case['id']} block={block} pointwise={count} passed", flush=True)
    profile_scores = {case['id']: {str(count): float(np.exp(np.mean(np.log([
        np.median(row['report']['samplesSeconds']) for row in runs
        if row['profile'] == case['id'] and row['pointwiseWorkers'] == count])))) for count in counts} for case in fixture['cases']}
    scores = {count: float(np.exp(np.mean(np.log([profile[str(count)] for profile in profile_scores.values()])))) for count in counts}
    selected = min(count for count in counts if scores[count] <= min(scores.values()) * 1.01)
    unchanged = source == source_receipt(args.source.resolve()) and worker_hash == digest(worker)
    save(args.output / 'selection.json', dict(selectedPointwiseWorkers=selected, scores=scores, profileScores=profile_scores,
         sourcesAndWorkerUnchanged=unchanged, defaultAdoption=False))
    if not unchanged: raise ValueError('Source or worker changed during calibration')


if __name__ == '__main__': main()
