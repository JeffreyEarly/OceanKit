#!/usr/bin/env python3
"""Independent-source compact flux qualification; model adoption is separate."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np

from run_derivative_qualification import compare, digest, save, source_receipt


def run_child(command, prefix):
    with prefix.with_suffix('.stdout').open('wb') as out, prefix.with_suffix('.stderr').open('wb') as err:
        child = subprocess.Popen(command, stdout=out, stderr=err)
        _, status, usage = os.wait4(child.pid, 0)
        child.returncode = os.waitstatus_to_exitcode(status)
    return {'returncode': child.returncode,
            'completeLifetimePeakRSSBytes': int(usage.ru_maxrss * (1 if sys.platform == 'darwin' else 1024))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('manifest', 'control_source', 'candidate_source', 'control_build', 'candidate_build', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--control-receipt', type=Path, required=True)
    parser.add_argument('--candidate-receipt', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=12)
    parser.add_argument('--pointwise-workers', type=int, default=1)
    parser.add_argument('--mode', choices=('smoke', 'pilot', 'final'), default='pilot')
    args = parser.parse_args()
    if min(args.workers, args.pointwise_workers) < 1:
        parser.error('Worker counts must be positive')
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads(args.manifest.read_text())
    blocks, warmups, samples = {'smoke': (1, 0, 1), 'pilot': (4, 2, 4), 'final': (8, 2, 4)}[args.mode]
    sources = {'control': args.control_source.resolve(), 'candidate': args.candidate_source.resolve()}
    builds = {'control': args.control_build.resolve(), 'candidate': args.candidate_build.resolve()}
    workers = {side: build / 'wv-variable-adoption' for side, build in builds.items()}
    receipts = {side: source_receipt(source) for side, source in sources.items()}
    if args.mode != 'smoke' and any(receipt['status'] for receipt in receipts.values()):
        raise ValueError('Timing requires clean, frozen source trees')
    worker_hashes = {side: digest(worker) for side, worker in workers.items()}
    build_receipts = {side: json.loads(path.read_text()) for side, path in
                      (('control', args.control_receipt), ('candidate', args.candidate_receipt))}
    for side, receipt in build_receipts.items():
        cache = (builds[side] / 'CMakeCache.txt').read_text()
        configured_source = next(line.split('=', 1)[1] for line in cache.splitlines() if line.startswith('WVM_SOURCE:PATH='))
        if (Path(configured_source).resolve() != sources[side]
                or receipt['source']['commit'] != receipts[side]['commit']
                or receipt['source']['files'] != receipts[side]['files']
                or receipt['worker']['sha256'] != worker_hashes[side]
                or not receipt.get('harness')):
            raise ValueError('Source/build/worker identity mismatch: ' + side)
    provenance = dict(schema='wvm-variable-compact-qualification-v1', mode=args.mode,
                      boundary='complete-direct-flux', host=platform.platform(), sources=receipts, buildReceipts=build_receipts,
                      workers={side: dict(path=str(worker), sha256=worker_hashes[side]) for side, worker in workers.items()},
                      fixtureManifest=manifest, manifestSHA256=digest(args.manifest), runnerSHA256=digest(__file__),
                      helperSHA256=digest(Path(__file__).with_name('run_derivative_qualification.py')),
                      buildCaches={side: (build / 'CMakeCache.txt').read_text() for side, build in builds.items()},
                      blocks=blocks, warmups=warmups, samples=samples, horizontalWorkers=args.workers,
                      pointwiseWorkers=args.pointwise_workers,
                      environment={k: v for k, v in os.environ.items() if k.startswith(('VECLIB_', 'OMP_', 'OPENBLAS_', 'MKL_', 'DYLD_'))},
                      payloadRetention='First successful payload per profile/selection and every failure; verified duplicates hashed then removed.')
    save(args.output / 'provenance.json', provenance)
    selections = [('production', 'control', 'frozen'), ('prior-pruned', 'control', 'pruned-streamed'),
                  ('interleaved', 'candidate', 'pruned-streamed'), ('compact', 'candidate', 'compact')]
    runs, retained, profiles = [], {}, {}
    for case in manifest['cases']:
        for path, sha in (('sourcePath', 'sourceSHA256'), ('matlabFluxPath', 'matlabFluxSHA256')):
            if digest(case[path]) != case[sha]:
                raise ValueError('Fixture identity changed: ' + case[path])
        for block in range(blocks):
            for label, side, selection in selections if block % 2 == 0 else reversed(selections):
                name = f"{case['id']}-{block}-{label}"
                prefix = args.output / name
                payload = prefix.with_suffix('.bin')
                command = [str(workers[side]), case['sourcePath'], selection, str(args.workers), str(warmups), str(samples), str(payload)]
                if side == 'candidate':
                    command.append(str(args.pointwise_workers))
                row = dict(id=name, profile=case['id'], block=block, selection=label, command=command)
                runs.append(row)
                save(args.output / 'runs.json', runs)
                row.update(run_child(command, prefix))
                save(args.output / 'runs.json', runs)
                if row['returncode']:
                    raise RuntimeError('Worker failed; evidence retained: ' + name)
                r = json.loads(prefix.with_suffix('.stdout').read_text())
                row['report'] = r
                family = {'stratified-qg': 'WVTransformStratifiedQG', 'hydrostatic': 'WVTransformHydrostatic', 'boussinesq': 'WVTransformBoussinesq'}[case['family']]
                expected_bytes = case['Nj'] * case['Nkl'] * (1 if case['family'] == 'stratified-qg' else 3) * 16
                if (payload.stat().st_size != expected_bytes or r['schema'] != 'wvm-variable-screening-v1'
                        or r['family'] != family or r['grid'] != case['grid'] or r['selection'] != selection
                        or r['horizontalSchedule'] != ('full-fft-gather' if selection == 'frozen' else 'fftw-streaming-pruned-tile16')
                        or r['horizontalWorkers'] != (1 if selection == 'frozen' else args.workers)
                        or r['matrixBackend'] != 'accelerate' or r['provider']['fftThreads'] != 1
                        or r['warmups'] != warmups or len(r['samplesSeconds']) != samples
                        or not all(np.isfinite(t) and t > 0 for t in r['samplesSeconds'])):
                    raise ValueError('Worker metadata mismatch: ' + name)
                if side == 'candidate' and (r['compactSplitViews'] != (selection == 'compact')
                        or not r['streamedNonlinear'] or r['pointwiseWorkers'] != args.pointwise_workers
                        or r['matrixBackendIdentifier'] != 'accelerate-gemm-v1'):
                    raise ValueError('Candidate execution policy mismatch: ' + name)
                for key in ('baseLibrary', 'threadLibrary'):
                    r['provider'][key + 'SHA256'] = digest(r['provider'][key])
                row['payloadSHA256'] = digest(payload)
                production = retained.get((case['id'], 'production'), payload)
                if production == payload and label != 'production':
                    raise ValueError('Production comparison missing')
                row['matlabComparison'] = compare(payload, Path(case['matlabFluxPath']), True)
                row['productionComparison'] = compare(payload, production, True)
                save(args.output / 'runs.json', runs)
                if not all(row[key]['passed'] for key in ('matlabComparison', 'productionComparison')):
                    raise ValueError('Scientific mismatch: ' + name)
                key = (case['id'], label)
                row['payloadRetained'] = key not in retained
                if key in retained:
                    payload.unlink()
                else:
                    retained[key] = payload
                save(args.output / 'runs.json', runs)
                print(name + ' passed', flush=True)
        selected = {label: [r for r in runs if r['profile'] == case['id'] and r['selection'] == label] for label, _, _ in selections}
        metrics = {label: dict(processMediansSeconds=[float(np.median(r['report']['samplesSeconds'])) for r in rows],
                              maximumCompleteLifetimePeakRSSBytes=max(r['completeLifetimePeakRSSBytes'] for r in rows),
                              persistentBytes=[r['report']['persistentBytes'] for r in rows]) for label, rows in selected.items()}
        ratios = {label: np.array(metrics['compact']['processMediansSeconds']) / metrics[label]['processMediansSeconds'] for label in metrics if label != 'compact'}
        profiles[case['id']] = dict(metrics=metrics, compactTimeRatios={label: float(np.exp(np.mean(np.log(values)))) for label, values in ratios.items()})
        save(args.output / 'summary.json', dict(mode=args.mode, decision='Flux evidence only; complete model gates remain separate', profiles=profiles))
    unchanged = (receipts == {side: source_receipt(source) for side, source in sources.items()}
                 and worker_hashes == {side: digest(worker) for side, worker in workers.items()})
    save(args.output / 'postflight.json', dict(sourcesAndWorkersUnchanged=unchanged))
    if not unchanged:
        raise ValueError('Sources or binaries changed during campaign')


if __name__ == '__main__':
    main()
