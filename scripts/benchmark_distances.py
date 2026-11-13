#!/usr/bin/env python3
"""Compare prepared FASTA variants using the installed Virtitta distance engine.

Manifest JSON rows: sample_run_id, method, mode, subtype, fasta, mask.
FASTA and JSON boolean mask must already have the SAME documented trimming applied.
Each input is one record; masks are positional evidence, never inferred from letters.
"""
import argparse
import csv
from dataclasses import asdict
import itertools
import json
import re
from pathlib import Path
import subprocess
import statistics
import sys

from benchmark_consensus import digest, write_json
from consensus_metadata import sequence


def aligned_records(path):
    result = {}
    identifier = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            identifier = line[1:].split()[0]
            result[identifier] = ''
        elif line.strip():
            if identifier is None:
                raise ValueError('Missing FASTA header')
            result[identifier] += line.strip()
    return result


def nearest(cells):
    result = {}
    for (first, second), cell in cells.items():
        if cell['compared_bases']:
            for a, b in ((first, second), (second, first)):
                result.setdefault(a, []).append((cell['total_events'], b))
    return {a: sorted(b for value, b in values if value == min(x[0] for x in values))
            for a, values in result.items()}


def compare_cells(old, new):
    changes = []
    for pair in sorted(old.keys() & new.keys()):
        a, b = old[pair], new[pair]
        delta = {key: b[key]-a[key] for key in a}
        changes.append(dict(first=pair[0], second=pair[1], baseline=a, candidate=b,
                            delta=delta, comparable=bool(a['compared_bases'] and b['compared_bases']),
                            baseline_substitution_density=a['substitutions']/a['compared_bases'] if a['compared_bases'] else None,
                            candidate_substitution_density=b['substitutions']/b['compared_bases'] if b['compared_bases'] else None))
    shared = old.keys() & new.keys()
    old_nn = nearest({key: old[key] for key in shared})
    new_nn = nearest({key: new[key] for key in shared})
    return dict(pairs=changes, missing_pairs=dict(baseline_only=len(old.keys()-new.keys()),
                                                 candidate_only=len(new.keys()-old.keys())),
                nearest_neighbour_changes={
        k: dict(baseline=old_nn.get(k, []), candidate=new_nn.get(k, []))
        for k in sorted(old_nn.keys() | new_nn.keys()) if old_nn.get(k) != new_nn.get(k)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--virtitta-source', type=Path, required=True)
    parser.add_argument('--mafft', default='mafft')
    parser.add_argument('--baseline', default='baseline')
    parser.add_argument('--plots', action='store_true', help='Write distance-difference heatmaps (requires matplotlib)')
    args = parser.parse_args()
    sys.path.insert(0, str(args.virtitta_source))
    from virtitta.distance import pair_distance, project_mask, project_sequence

    rows = json.loads(args.manifest.read_text())
    keys = [(r['sample_run_id'], r['method'], r['mode']) for r in rows]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate sample/method/mode')
    args.output.mkdir(parents=True, exist_ok=False)
    groups = {}
    for index, row in enumerate(rows):
        if row['mode'] not in ('iupac', 'majority'):
            raise ValueError('Mode must be iupac or majority')
        if not re.fullmatch(r'[A-Za-z0-9_-]+', row['method']) or row['method'] == 'pooled':
            raise ValueError('Method names must be safe identifiers other than pooled')
        _, seq = sequence(row['fasta'])
        mask = json.loads(Path(row['mask']).read_text())
        if len(mask) != len(seq) or any(type(x) is not bool for x in mask):
            raise ValueError('Mask must be boolean and match FASTA coordinates')
        if '-' in seq:
            raise ValueError('Manifest requires unaligned FASTAs')
        row = dict(row, token=f's{index}', sequence=seq, coverage=mask)
        groups.setdefault((row['subtype'] or 'uncertain', row['mode']), []).append(row)
    commands, comparisons, incompatible = [], [], []
    row_lookup = {(r['sample_run_id'], r['method'], r['mode']): r for g in groups.values() for r in g}
    for index, ((subtype, mode), group) in enumerate(sorted(groups.items())):
        directory = args.output/f'group_{index}'
        directory.mkdir()
        methods = sorted({r['method'] for r in group})
        # Pooled alignment gives all methods a shared gap pattern. Operational
        # alignments are calculated independently, without forcing equal lengths.
        matrices = {}
        for method in ['pooled']+methods:
            selected = group if method == 'pooled' else [r for r in group if r['method'] == method]
            templates = {r['token']: r['sequence'] for r in selected}
            if method != 'pooled' and mode == 'majority':
                partners = [row_lookup.get((r['sample_run_id'], method, 'iupac')) for r in selected]
                if any(p is None or len(p['sequence']) != len(r['sequence'])
                       for r, p in zip(selected, partners)):
                    incompatible.append(dict(subtype=subtype, method=method, mode=mode,
                                             reason='Missing or length-incompatible IUPAC/majority pair'))
                    continue
                templates = {r['token']: p['sequence'] for r, p in zip(selected, partners)}
            source = directory/f'{method}.input.fasta'
            source.write_text(''.join('>'+r['token']+'\n'+templates[r['token']]+'\n' for r in selected))
            destination = directory/f'{method}.aligned.fasta'
            command = [args.mafft, '--auto', '--thread', '1', str(source)]
            with destination.open('w') as output, (directory/f'{method}.log').open('w') as log:
                subprocess.run(command, stdout=output, stderr=log, check=True)
            commands.append(command)
            aligned = aligned_records(destination)
            if set(aligned) != {r['token'] for r in selected}:
                raise ValueError('MAFFT identifiers differ from inputs')
            for r in selected:
                if aligned[r['token']].replace('-', '').upper() != templates[r['token']]:
                    raise ValueError('MAFFT changed sequence payload')
                if templates[r['token']] != r['sequence']:
                    aligned[r['token']] = project_sequence(aligned[r['token']], r['sequence'], validate_bases=False)
            masks = {r['token']: project_mask(aligned[r['token']], r['coverage']) for r in selected}
            for target in methods if method == 'pooled' else [method]:
                members = {r['sample_run_id']: r for r in selected if r['method'] == target}
                cells = {}
                for a, b in itertools.combinations(sorted(members), 2):
                    x, y = members[a]['token'], members[b]['token']
                    cells[a, b] = asdict(pair_distance(aligned[x], aligned[y], masks[x], masks[y]))
                matrices[method, target] = cells
                identifiers = sorted(members)
                with (directory/f'{method}-{target}.matrix.tsv').open('w') as stream:
                    writer = csv.writer(stream, delimiter='\t')
                    writer.writerow(['sample_run_id']+identifiers)
                    for a in identifiers:
                        values = []
                        for b in identifiers:
                            cell = cells.get(tuple(sorted((a, b))))
                            values.append(0 if a == b else cell['total_events'] if cell['compared_bases'] else -1)
                        writer.writerow([a]+values)
            if method == 'pooled':
                pooled, pooled_masks = aligned, masks
        for candidate in methods:
            if candidate == args.baseline:
                continue
            for alignment in ('operational', 'pooled'):
                if alignment == 'operational' and ((candidate, candidate) not in matrices
                        or (args.baseline, args.baseline) not in matrices):
                    continue
                old = matrices[(args.baseline if alignment == 'operational' else 'pooled'), args.baseline]
                new = matrices[(candidate if alignment == 'operational' else 'pooled'), candidate]
                comparisons.append(dict(subtype=subtype, mode=mode, candidate=candidate,
                                        alignment=alignment, **compare_cells(old, new)))
            # Hold callable nucleotide positions constant for each sample pair.
            members = {(r['sample_run_id'], r['method']): r['token'] for r in group}
            common_old, common_new = {}, {}
            samples = sorted({r['sample_run_id'] for r in group})
            for a, b in itertools.combinations(samples, 2):
                identities = [(s, m) for m in (args.baseline, candidate) for s in (a, b)]
                if not all(key in members for key in identities):
                    continue
                tokens = [members[key] for key in identities]
                common = [all(pooled_masks[t][i] and pooled[t][i].upper() in 'ACGTRYSWKMBDHV'
                              for t in tokens) for i in range(len(pooled[tokens[0]]))]
                for result, x, y in ((common_old, tokens[0], tokens[1]), (common_new, tokens[2], tokens[3])):
                    # Common-base comparison intentionally excludes gaps; native
                    # event counts above retain indels and their flank semantics.
                    result[a, b] = asdict(pair_distance(pooled[x], pooled[y], common, common))
            comparisons.append(dict(subtype=subtype, mode=mode, candidate=candidate,
                                    alignment='common_callable_bases', **compare_cells(common_old, common_new)))
    write_json(args.output/'comparisons.json', comparisons)
    write_json(args.output/'incompatible.json', incompatible)
    summaries = []
    with (args.output/'pair-differences.tsv').open('w') as stream:
        fields = ['subtype', 'mode', 'candidate', 'alignment', 'first', 'second', 'comparable',
                  'substitutions', 'indel_events', 'total_events', 'compared_bases']
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for result in comparisons:
            tag = {k: result[k] for k in fields[:4]}
            for pair in sorted(result['pairs'], key=lambda p: abs(p['delta']['total_events']), reverse=True):
                writer.writerow(tag | {k: pair[k] for k in ('first', 'second', 'comparable')} | pair['delta'])
            metrics = {}
            for metric in fields[7:]:
                values = sorted(abs(p['delta'][metric]) for p in result['pairs'] if p['comparable'])
                metrics[metric] = dict(median=statistics.median(values) if values else None,
                    p95=values[max(0, (95*len(values)+99)//100-1)] if values else None,
                    maximum=max(values) if values else None,
                    changed_pairs=sum(v != 0 for v in values))
            summaries.append(tag | dict(metrics=metrics, nearest_neighbour_changes=len(result['nearest_neighbour_changes'])))
    write_json(args.output/'summary.json', summaries)
    if args.plots:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        for index, result in enumerate(comparisons):
            ids = sorted({p[k] for p in result['pairs'] for k in ('first', 'second')})
            if not ids:
                continue
            lookup = {identifier: n for n, identifier in enumerate(ids)}
            fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
            for axis, metric in zip(axes, ('total_events', 'compared_bases')):
                values = np.full((len(ids), len(ids)), np.nan)
                np.fill_diagonal(values, 0)
                for pair in result['pairs']:
                    if not pair['comparable']:
                        continue
                    a, b = lookup[pair['first']], lookup[pair['second']]
                    values[a, b] = values[b, a] = pair['delta'][metric]
                bound = max(1, np.nanmax(np.abs(values)))
                artist = axis.imshow(values, cmap='coolwarm', vmin=-bound, vmax=bound)
                axis.set_title('Candidate − baseline: '+metric.replace('_', ' '))
                if len(ids) <= 40:
                    axis.set_xticks(range(len(ids)), ids, rotation=90, fontsize=6)
                    axis.set_yticks(range(len(ids)), ids, fontsize=6)
                fig.colorbar(artist, ax=axis, shrink=0.7)
            fig.suptitle(f"{result['subtype']} / {result['mode']} / {result['candidate']} / {result['alignment']}")
            fig.savefig(args.output/f'difference-{index}.png', dpi=150)
            plt.close(fig)
    write_json(args.output/'commands.json', commands)
    write_json(args.output/'provenance.json', dict(manifest_sha256=digest(args.manifest),
        inputs=[{k: r[k] for k in ('sample_run_id', 'method', 'mode', 'subtype')} |
                {k+'_sha256': digest(r[k]) for k in ('fasta', 'mask')} for r in rows],
        virtitta_commit=subprocess.check_output(['git', '-C', str(args.virtitta_source), 'rev-parse', 'HEAD']).decode().strip(),
        distance_source_sha256=digest(args.virtitta_source/'virtitta/distance.py'),
        mafft_version=subprocess.check_output([args.mafft, '--version'], stderr=subprocess.STDOUT).decode()))


if __name__ == '__main__':
    main()
