"""Run the portable lab wrapper against the engine's synthetic genome fixture."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from test_engine_synthetic import synthetic, input_paths

ROOT = Path(__file__).resolve().parents[1]


def test_presplit_wrapper_end_to_end(synthetic, tmp_path):
    _, genome_root = synthetic
    paths = input_paths(tmp_path)
    # Exercise chr-prefix reconciliation against the existing fixture's .fai.
    paths['up'].write_text(paths['up'].read_text().replace('chr1', '1'))
    motifs = tmp_path / 'motifs.tsv'
    motifs.write_text('Protein_name\tregularExpression\nTEST\tAA\nOTHER\tATAT\n')
    output = tmp_path / 'output'
    config = {
        'arm': 'synthetic',
        'inputs': {key: path.name for key, path in paths.items()},
        'genome': {'root': str(genome_root), 'build': 'synthetic'},
        'motifs': {'known': motifs.name},
        'engine': {'stat_method': 'fisher', 'workers': 1, 'intron': 20,
                   'exon': 10, 'window': 4, 'step': 2},
        'blas_threads': 1,
        'summary': {'perms': 20, 'seed': 149, 'workers': 1},
        'output_root': output.name,
    }
    source = tmp_path / 'config.json'
    source.write_text(json.dumps(config))
    command = [sys.executable, str(ROOT / 'tools' / 'rmaps3_lab_run.py'),
               '--config', str(source), '--perms', '20']
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                            env=dict(os.environ, RMAPS_FORCE_MOTIF_FALLBACK='1'),
                            timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    for name in ('run_manifest.json', 'engine/run_manifest.json',
                 'engine/pVal.up.vs.bg.RNAmap.txt',
                 'engine/pVal.dn.vs.bg.RNAmap.txt',
                 'summary/condensed_per_rbp.tsv'):
        assert (output / name).is_file(), name
    assert list((output / 'engine' / 'positional').glob('*.hits.npz'))
    with (output / 'summary' / 'condensed_per_rbp.tsv').open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
        assert {'arm', 'RBP', 'direction', 'pooled_region', 'selected_motif_key',
                'calib_p_pooled', 'calib_q', 'native_p_pooled',
                'fg_proportion', 'bg_proportion', 'enrichment_ratio',
                'n_fg_hit', 'n_fg_elig', 'n_bg_hit', 'n_bg_elig'} <= set(reader.fieldnames)
    assert rows and {row['RBP'] for row in rows} == {'TEST', 'OTHER'}
    assert {row['arm'] for row in rows} == {'synthetic'}
    manifest = json.loads((output / 'run_manifest.json').read_text())
    assert manifest['status'] == 'complete'
    assert manifest['config']['sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert manifest['perms_override'] == 20
    assert manifest['environment']['OPENBLAS_NUM_THREADS'] == '1'
    steps = {step['name']: step for step in manifest['steps']}
    assert {'scaffold_mapping', 'motif_map', 'summarize'} <= steps.keys()
    assert all(step['status'] == 'complete' for step in steps.values())
    summary_inputs = {Path(item['path']).name for item in steps['summarize']['inputs']}
    assert {'pVal.up.vs.bg.RNAmap.txt', 'pVal.dn.vs.bg.RNAmap.txt'} <= summary_inputs
    # Refusal must not rewrite completed artifacts.
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest()
              for path in output.rglob('*') if path.is_file()}
    refusal = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=15)
    assert refusal.returncode != 0
    assert 'non-empty' in (refusal.stdout + refusal.stderr).lower()
    assert before == {path: hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in output.rglob('*') if path.is_file()}


def test_wrapper_rejects_arm_path_escape(tmp_path):
    config = {'arm': '../outside', 'inputs': {}, 'genome': {}, 'motifs': {},
              'output_root': str(tmp_path / 'output')}
    source = tmp_path / 'config.json'
    source.write_text(json.dumps(config))
    result = subprocess.run([sys.executable, str(ROOT / 'tools' / 'rmaps3_lab_run.py'),
                             '--config', str(source)], cwd=ROOT, capture_output=True,
                            text=True, timeout=15)
    assert result.returncode != 0
    assert 'arm' in (result.stdout + result.stderr).lower()
    assert not (tmp_path / 'output').exists()
