#!/usr/bin/env python3
"""Apply the frozen variable adoption gates; retain all comparator decisions."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

CONTROL_COMMIT = 'b157988355d2fabb6cbc62d56de18b203456897f'
CANDIDATE_COMMIT = 'a3acbf4dc9305c4d656399e9f98bc8360ddfa66f'
FLUX_PROFILES = {
    f'{family}-{n}x{n}x{nz}': (family, [n, n, nz])
    for family in ('stratified-qg', 'hydrostatic', 'boussinesq')
    for n, nz in ((256, 129), (512, 257))
    if family != 'boussinesq' or n == 256}
MODEL_PROFILES = {'variable-stratified-qg-composite-large', 'variable-hydrostatic-composite-large'}
BOUSS_MODEL_PROFILES = {'variable-boussinesq-composite-256'}


def positive(value):
    return np.isfinite(value) and value > 0


def paired_summary(profiles):
    """Resample paired process observations within each equally weighted profile."""
    rng = np.random.default_rng(455)
    draws = np.zeros(10000)
    rows = {}
    for name, ratios in sorted(profiles.items()):
        values = np.log(np.asarray(ratios, dtype=float))
        if not len(values) or not np.all(np.isfinite(values)):
            raise ValueError('Invalid paired observations: ' + name)
        draws += values[rng.integers(0, len(values), (10000, len(values)))].mean(axis=1) / len(profiles)
        rows[name] = {'geometricRatio': float(np.exp(values.mean())), 'pairedRatios': list(ratios)}
    return {'geometricRatio': float(np.exp(np.mean([np.log(r['geometricRatio']) for r in rows.values()]))),
            'confidence95': np.exp(np.quantile(draws, [.025, .975])).tolist(),
            'profiles': rows, 'resamples': 10000, 'seed': 455}


def flux_summary(folder):
    provenance = json.loads((folder / 'provenance.json').read_text())
    runs = json.loads((folder / 'runs.json').read_text())
    post = json.loads((folder / 'postflight.json').read_text())
    profiles = {c['id'] for c in provenance['fixtureManifest']['cases']}
    if (profiles != set(FLUX_PROFILES) or len(provenance['fixtureManifest']['cases']) != 5
            or any((c['family'], c['grid']) != FLUX_PROFILES[c['id']]
                   for c in provenance['fixtureManifest']['cases'])
            or provenance['schema'] != 'wvm-variable-compact-qualification-v1'
            or provenance['horizontalWorkers'] != 12 or provenance['pointwiseWorkers'] != 8):
        raise ValueError('Direct-flux profile or policy differs from frozen contract')
    for side, commit in (('control', CONTROL_COMMIT), ('candidate', CANDIDATE_COMMIT)):
        source, build = provenance['sources'][side], provenance['buildReceipts'][side]
        if (source['commit'] != commit or source['status']
                or build['source']['commit'] != commit or build['source']['files'] != source['files']
                or build['worker']['sha256'] != provenance['workers'][side]['sha256']
                or not build.get('harness')):
            raise ValueError('Direct-flux frozen source/build receipt mismatch: ' + side)
    labels = ('production', 'prior-pruned', 'interleaved', 'compact')
    indexed = {(r['profile'], r['block'], r['selection']): r for r in runs}
    expected = {(p, b, s) for p in profiles for b in range(8) for s in labels}
    complete = (provenance['mode'] == 'final' and provenance['blocks'] == 8
                and provenance['warmups'] == 2 and provenance['samples'] == 4
                and set(indexed) == expected and len(indexed) == len(runs)
                and post['sourcesAndWorkersUnchanged'])
    if not complete:
        raise ValueError('Incomplete final direct-flux protocol')
    provider_identities = {tuple(r['report']['provider'][key] for key in
                                ('version', 'baseLibrary', 'baseLibrarySHA256', 'threadLibrary', 'threadLibrarySHA256'))
                           for r in runs}
    if len(provider_identities) != 1:
        raise ValueError('Direct-flux provider changed within or across source selections')
    for row in runs:
        report = row['report']
        selection = {'production': 'frozen', 'prior-pruned': 'pruned-streamed',
                     'interleaved': 'pruned-streamed', 'compact': 'compact'}[row['selection']]
        if (report['selection'] != selection or report['grid'] != FLUX_PROFILES[row['profile']][1]
                or report['warmups'] != 2 or len(report['samplesSeconds']) != 4
                or not all(positive(t) for t in report['samplesSeconds'])
                or not positive(row['completeLifetimePeakRSSBytes']) or not positive(report['persistentBytes'])
                or report['matrixBackend'] != 'accelerate' or report['provider']['fftThreads'] != 1
                or report['horizontalWorkers'] != (1 if selection == 'frozen' else 12)
                or (row['selection'] in ('compact', 'interleaved') and
                    (report['pointwiseWorkers'] != 8 or not report['streamedNonlinear']
                     or report['compactSplitViews'] != (selection == 'compact')
                     or report['matrixBackendIdentifier'] != 'accelerate-gemm-v1'))):
            raise ValueError('Direct-flux runtime observation mismatch: ' + row['id'])
    scientific = all(r['returncode'] == 0 and r['matlabComparison']['passed']
                     and r['productionComparison']['passed'] for r in runs)
    result = {'scientificChecksPassed': scientific, 'completeProtocol': complete,
              'comparators': {}}
    for comparator in labels[:-1]:
        time, rss, owned = {}, {}, {}
        for profile in sorted(profiles):
            time[profile], rss[profile], owned[profile] = [], [], []
            for block in range(8):
                candidate = indexed[profile, block, 'compact']
                control = indexed[profile, block, comparator]
                time[profile].append(float(np.median(candidate['report']['samplesSeconds']) /
                                           np.median(control['report']['samplesSeconds'])))
                rss[profile].append(candidate['completeLifetimePeakRSSBytes'] / control['completeLifetimePeakRSSBytes'])
                owned[profile].append(candidate['report']['persistentBytes'] / control['report']['persistentBytes'])
        timing, memory = paired_summary(time), paired_summary(rss)
        gates = {'scientific': scientific,
                 'noProfileTimeRegression': all(r['geometricRatio'] <= 1.03 for r in timing['profiles'].values()),
                 'noProfileRSSGrowth': all(r['geometricRatio'] <= 1 for r in memory['profiles'].values()),
                 'noOwnedMemoryGrowth': all(r <= 1 for values in owned.values() for r in values)}
        if comparator == 'production':
            gates.update(geometricFluxImprovement=timing['geometricRatio'] <= .9,
                         pairedIntervalExcludesTie=timing['confidence95'][1] < 1)
        result['comparators'][comparator] = {'timing': timing, 'peakRSS': memory,
                                            'ownedCapacityRatios': owned, 'gates': gates,
                                            'passed': all(gates.values())}
    result['productionGatesPassed'] = result['comparators']['production']['passed']
    result['priorPrunedGatesPassed'] = result['comparators']['prior-pruned']['passed']
    return result


def model_summary(path):
    campaign = json.loads(path.read_text())
    manifest_bytes = Path(campaign['manifestPath']).read_bytes()
    manifest = json.loads(manifest_bytes)
    expected_profiles = (BOUSS_MODEL_PROFILES if manifest['schema'] == 'wvm-variable-boussinesq-model-fixture-256-v1'
                         else MODEL_PROFILES)
    if (hashlib.sha256(manifest_bytes).hexdigest() != campaign['manifestSHA256']
            or manifest['schema'] not in ('wvm-variable-model-fixtures-large-v1', 'wvm-variable-boussinesq-model-fixture-256-v1')
            or manifest['parameters']['grid'] != [256, 256, 129]
            or {p['id'] for p in manifest['profiles']} != expected_profiles):
        raise ValueError('Model fixture contract differs from the final freeze')
    blocks = campaign['blocks']
    profiles = {b['profile'] for b in blocks}
    selected_profiles = set(campaign.get('selectedProfiles', expected_profiles))
    expected = {(p, b) for p in profiles for b in range(10)}
    if (campaign['mode'] != 'final' or not selected_profiles or not selected_profiles <= expected_profiles
            or profiles != selected_profiles or len(blocks) != len(expected)
            or {(b['profile'], b['block']) for b in blocks} != expected
            or any(b['warmup'] != (b['block'] < 2) for b in blocks)
            or not campaign['summary']['passed'] or not campaign['summary']['sourceBuildAndProviderIdentityStable']):
        raise ValueError('Incomplete or scientifically failing final model protocol')
    for block in blocks:
        labels = {'control-frozen', 'candidate-frozen', 'candidate-interleaved', 'candidate-compact'}
        comparisons = {'matlab-vs-' + label for label in labels}
        if campaign['retentionPolicy'] == 'first-block' and block['block'] > 0:
            comparisons |= {'retained-control-vs-' + label for label in labels}
        else:
            comparisons |= {'control-vs-' + label for label in labels if label != 'control-frozen'}
        if (not block['passed'] or set(block['runs']) !=
                labels or set(block['comparisons']) != comparisons
                or not all(c['passed'] for c in block['comparisons'].values())):
            raise ValueError('Model scientific comparison inventory mismatch')
        for key, run in block['runs'].items():
            report = run['report']
            expected_commit = CONTROL_COMMIT if key == 'control-frozen' else CANDIDATE_COMMIT
            if (not run['passed'] or run['process']['exitCode'] != 0
                    or report['source']['commit'] != expected_commit or report['source']['compiledSourceDirty']
                    or report['work']['acceptedSteps'] != 1 or report['work']['rhsEvaluations'] != 4
                    or report['work']['scalarAdvections'] != 4 or report['work']['denseOutputEvaluations'] != 1
                    or report['selection']['fftwInternalThreads'] != 1
                    or not positive(report['timingSeconds']['integrate'])
                    or not positive(run['process']['completeLifetimePeakRSSBytes'])
                    or not positive(report['storageBytes']['fullModelRetainedBeforeClose'])):
                raise ValueError('Model observed execution/source mismatch: ' + key)
            if key in ('candidate-compact', 'candidate-interleaved') and (
                    report['selection']['horizontalWorkers'] != 12 or report['selection']['pointwiseWorkers'] != 8
                    or report['selection']['matrixBackendIdentifier'] != 'accelerate-gemm-v1'):
                raise ValueError('Model observed optimized topology mismatch')
    measured = [b for b in blocks if not b['warmup']]
    result = {'scientificChecksPassed': True, 'completeProtocol': True,
              'profiles': sorted(profiles), 'selections': {}}
    for selection in ('frozen', 'interleaved', 'compact'):
        timing, rss, owned = {}, {}, {}
        for profile in sorted(profiles):
            timing[profile], rss[profile], owned[profile] = [], [], []
            for block in measured:
                if block['profile'] != profile:
                    continue
                control, candidate = block['runs']['control-frozen'], block['runs']['candidate-' + selection]
                timing[profile].append(candidate['report']['timingSeconds']['integrate'] /
                                       control['report']['timingSeconds']['integrate'])
                rss[profile].append(candidate['process']['completeLifetimePeakRSSBytes'] /
                                    control['process']['completeLifetimePeakRSSBytes'])
                owned[profile].append(candidate['report']['storageBytes']['fullModelRetainedBeforeClose'] /
                                      control['report']['storageBytes']['fullModelRetainedBeforeClose'])
        time, memory = paired_summary(timing), paired_summary(rss)
        gates = {'noProfileTimeRegression': all(r['geometricRatio'] <= 1.03 for r in time['profiles'].values()),
                 'noProfileRSSGrowth': all(r['geometricRatio'] <= 1 for r in memory['profiles'].values()),
                 'noOwnedMemoryGrowth': all(r <= 1 for values in owned.values() for r in values)}
        result['selections'][selection] = {'integration': time, 'peakRSS': memory,
                                           'ownedCapacityRatios': owned, 'gates': gates,
                                           'passed': all(gates.values())}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--flux', type=Path, required=True)
    parser.add_argument('--model', type=Path, action='append')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = {'schema': 'wvm-variable-final-gates-v1', 'flux': flux_summary(args.flux),
              'evaluatorSHA256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'defaultAdoption': False,
              'limitations': ['Boussinesq 512 direct flux excluded by capacity.',
                              'Model qualification uses short 256x256x129 composite workloads; 512 model workloads are not qualified.',
                              'Public MATLAB-loaded variable service injection is not exposed.']}
    if args.model:
        models = [model_summary(path) for path in args.model]
        observed = [p for model in models for p in model['profiles']]
        complete = set(observed) == MODEL_PROFILES | BOUSS_MODEL_PROFILES and len(observed) == 3
        result['model'] = {'campaigns': models, 'allThreeFamiliesQualified': complete}
        result['measuredGatesPassed'] = (result['flux']['productionGatesPassed']
                                        and result['flux']['priorPrunedGatesPassed']
                                        and complete
                                        and all(model['selections']['compact']['passed'] for model in models))
    result['decision'] = 'Measurement evidence; final source/default/CI handoff remains separate.'
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('flux', 'model')}, indent=2))


if __name__ == '__main__':
    main()
