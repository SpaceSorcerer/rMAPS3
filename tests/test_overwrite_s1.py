import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rmaps_core.output_utils import ensure_output_directory


def snapshot(path):
    return {item.relative_to(path).as_posix(): item.read_bytes()
            for item in path.rglob('*') if item.is_file()}


def command(tmp_path, output):
    known = tmp_path / 'known.tsv'
    optional = tmp_path / 'optional.tsv'
    known.write_text('Protein_name\tregularExpression\nTEST\tAA\n')
    optional.write_text('Protein_name\tregularExpression\nEXTRA\tATAT\n')
    return [sys.executable, str(ROOT / 'cli.py'), 'motif-map', 'se',
            '--known-motifs', str(known), '--motifs', str(optional),
            '--fasta-root', str(tmp_path), '--genome', 'unused',
            '--output', str(output), '--overwrite']


def test_no_manifest_refuses_nonempty_overwrite_without_writes(tmp_path):
    # AUDIT S1: unknown ownership must fail before a log or archive is created.
    output = tmp_path / 'out'
    output.mkdir()
    (output / 'pVal.up.vs.bg.RNAmap.txt').write_text('USER SENTINEL')
    before = snapshot(output)
    result = subprocess.run(command(tmp_path, output), capture_output=True, text=True)
    assert result.returncode != 0
    assert 'without run_manifest.json' in result.stdout + result.stderr
    assert snapshot(output) == before
    assert not list(output.glob('_previous_*'))


@pytest.mark.parametrize('collision', [
    'pVal.up.vs.bg.RNAmap.txt', 'log.motifMap.txt',
    'exon/up.coord.txt', 'fasta/bg.TargetExon.fasta',
    'temp/sequence_metadata.npz', 'temp/sequence_region_7.npy',
    'temp/TEST.AA.pVal.up.vs.bg.txt',
    'positional/EXTRA.ATAT.TargetExon-3prime.hits.npz',
    'maps/SE.TEST-AA.pdf', 'maps/SE.EXTRA-ATAT.png', 'maps/TEST-AA.png', 'maps/jj.png',
])
def test_cli_unregistered_target_refuses_before_archival(tmp_path, collision):
    # AUDIT S1: include optional motifs and both native/fallback plot paths.
    output = tmp_path / 'out'
    output.mkdir()
    owned = output / 'old.txt'
    owned.write_text('REGISTERED')
    (output / 'run_manifest.json').write_text(json.dumps({
        'outputs': [{'path': 'old.txt', 'status': 'retained'}]}))
    target = output / collision
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text('USER SENTINEL')
    before = snapshot(output)
    result = subprocess.run(command(tmp_path, output), capture_output=True, text=True)
    assert result.returncode != 0
    assert 'unregistered output target' in result.stdout + result.stderr
    assert snapshot(output) == before
    assert not list(output.glob('_previous_*'))


def test_registered_archive_preserves_unrelated_files(tmp_path):
    # AUDIT S1: registered ownership allows archival without sweeping sentinels.
    (tmp_path / 'owned.txt').write_text('REGISTERED')
    (tmp_path / 'unrelated.txt').write_text('USER SENTINEL')
    (tmp_path / 'run_manifest.json').write_text(json.dumps({
        'outputs': [{'path': 'owned.txt', 'status': 'retained'}]}))
    ensure_output_directory(tmp_path, overwrite=True, target_paths=['owned.txt'])
    archives = list(tmp_path.glob('_previous_*'))
    assert len(archives) == 1
    assert (archives[0] / 'owned.txt').read_text() == 'REGISTERED'
    assert (tmp_path / 'unrelated.txt').read_text() == 'USER SENTINEL'


def test_deleted_registration_does_not_claim_new_user_file(tmp_path):
    # AUDIT S1: a previously deleted target can have a new owner.
    (tmp_path / 'owned.txt').write_text('USER SENTINEL')
    (tmp_path / 'run_manifest.json').write_text(json.dumps({
        'outputs': [{'path': 'owned.txt', 'status': 'deleted'}]}))
    before = snapshot(tmp_path)
    with pytest.raises(ValueError, match='unregistered output target'):
        ensure_output_directory(tmp_path, overwrite=True, target_paths=['owned.txt'])
    assert snapshot(tmp_path) == before


def test_target_parent_file_refuses_before_archival(tmp_path):
    # AUDIT S1: a target directory occupied by a file must fail before archival.
    (tmp_path / 'temp').write_text('USER SENTINEL')
    (tmp_path / 'old.txt').write_text('REGISTERED')
    (tmp_path / 'run_manifest.json').write_text(json.dumps({
        'outputs': [{'path': 'old.txt', 'status': 'retained'}]}))
    before = snapshot(tmp_path)
    with pytest.raises(ValueError, match='parent is not a directory'):
        ensure_output_directory(tmp_path, overwrite=True,
                                target_paths=['temp/sequence_metadata.npz'])
    assert snapshot(tmp_path) == before


def test_xlsx_conversion_target_is_checked_before_writes(tmp_path):
    # AUDIT S1: preflight must precede even the wrapper's XLSX conversion.
    output = tmp_path / 'out'
    target = output / 'temp' / 'input.from_xlsx.tsv'
    target.parent.mkdir(parents=True)
    target.write_text('USER SENTINEL')
    (output / 'run_manifest.json').write_text(json.dumps({'outputs': []}))
    before = snapshot(output)
    result = subprocess.run(command(tmp_path, output) + ['--rMATS', str(tmp_path / 'input.xlsx')],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert 'unregistered output target' in result.stdout + result.stderr
    assert snapshot(output) == before
