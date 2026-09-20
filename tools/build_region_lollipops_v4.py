"""rMAPS3 SE region-resolved lollipops, version 4 (rank-sum layers).

Layer 1 / MAIN  : released_ranksum_rawP  - authors' released rMAPS3 (b9a9dce) run with
                  --stat-method mannwhitney; raw regional-minimum p exactly as the tool reports it.
Layer 2 / SUPP  : calibrated_ranksum     - the same statistic with a Westfall-Young
                  label-permutation p and BH q over the 756 pooled tests.

The look is inherited verbatim from v3.1 (build_rmaps_region_lollipops.py); only the layer
inputs change. v3.1 is never modified.
"""
from __future__ import annotations
import argparse, csv, hashlib, html, itertools, json, math, re, sys
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import openpyxl
from openpyxl.styles import Alignment, Font, PatternFill
from scipy.stats import rankdata

REGIONS = ['Upstream Intron', 'Exon Body', 'Downstream Intron']
DIRECTIONS = ['INCLUDED', 'SKIPPED']
LABELS = {'MIAT_KD_A': 'MIAT knockdown (rMATS events)',
          'MIAT_KD_B': 'MIAT knockdown (rMATS ∩ VAST-tools concordant events)',
          'MIAT_KD_Beffect': 'MIAT knockdown (rMATS ∩ VAST-tools concordant, relaxed rule)',
          'MIAT_OE_A': 'MIAT overexpression (rMATS events)',
          'QKI_KO_A': 'QKI knockout (rMATS events)',
          'QKI_KO_B': 'QKI knockout (rMATS ∩ VAST-tools concordant events)',
          'QKI_KO_Beffect': 'QKI knockout (rMATS ∩ VAST-tools concordant, relaxed rule)'}
COLORS = ['#0072B2', '#56B4E9', '#999999', '#DDDDDD']
LAYERS = ['released_ranksum_rawP', 'calibrated_ranksum']
LAYER_TITLES = {'released_ranksum_rawP': "MAIN — authors' rMAPS3 rank-sum, raw p",
                'calibrated_ranksum': 'SUPPLEMENT — permutation-calibrated rank-sum'}
KINDS = ['byRBP', 'byMotif']
VARIANTS = ['main', 'noSpliceosome_noBroad']
PNG_DPI = 200
SUBREGIONS = ['upstreamExon-3prime', 'upstreamExonIntron', 'upstreamIntron', 'targetExon-5prime',
              'targetExon-3prime', 'downstreamIntron', 'downstreamExonIntron', 'downstreamExon-5prime']
POOLING = {'Upstream Intron': (1, 2), 'Exon Body': (3, 4), 'Downstream Intron': (5, 6)}
RELEASED_COMMIT = 'b9a9dce'
STAT_METHOD = 'mannwhitney'
BH_FAMILY = 756
REPO_ROOT = Path(__file__).resolve().parents[1]
NAMING_TABLE = REPO_ROOT / 'data' / 'knownMotifs.human.mouse.txt'
ESRP_TABLE = REPO_ROOT / 'data' / 'ESRP.like.motif.txt'
_DIGEST_CACHE = {}


def digest(p, algo='md5'):
    p = Path(p)
    key = (str(p.resolve()), algo, p.stat().st_size, p.stat().st_mtime_ns)
    if key not in _DIGEST_CACHE:
        h = hashlib.new(algo)
        with p.open('rb') as f:
            for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
                h.update(block)
        _DIGEST_CACHE[key] = h.hexdigest()
    return _DIGEST_CACHE[key]


def read_tsv(p):
    with open(p, encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, fields, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows({k: ('NA' if v is None else v) for k, v in r.items()} for r in rows)


def rbp_label(row):
    return row.get('_display_name', 'ESRP-like' if re.fullmatch('motif_[0-9]+', row['RBP']) else row['RBP'])


def tick_label(row, kind):
    """by-RBP: the HGNC group. by-motif: the same, except the synthetic ESRP-like hexamers,
    which carry no gene name and are shown by sequence (lab figure convention)."""
    if kind == 'byMotif' and row['_naming_action'] == 'esrp_like_group':
        return f'{row["motif_key"].split(".", 1)[1]} (ESRP-like)'
    return row['_label']


# ---------------------------------------------------------------- naming (verbatim from v3.1)
def load_gene_names(gtf):
    mapping = {}
    with Path(gtf).open(encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t', 8)
            if len(fields) != 9 or fields[2] != 'gene':
                continue
            gene_id = re.search(r'gene_id "([^"]+)"', fields[8])
            gene_name = re.search(r'gene_name "([^"]+)"', fields[8])
            if gene_id and gene_name:
                key = gene_id.group(1).split('.')[0]
                if key in mapping and mapping[key] != gene_name.group(1):
                    raise ValueError(f'Conflicting GENCODE gene_name for {key}')
                mapping[key] = gene_name.group(1)
    if not mapping:
        raise ValueError(f'No gene names found in {gtf}')
    return mapping


def load_naming(table, esrp_table, gtf, alias_table):
    exact_symbols = set(load_gene_names(gtf).values())
    aliases = {}
    for row in read_tsv(alias_table):
        if not {'table_name', 'hgnc_symbol', 'evidence', 'ambiguity_note'} <= row.keys():
            raise ValueError(f'Missing alias-table columns in {alias_table}')
        name = row['table_name'].strip()
        symbol = row['hgnc_symbol'].strip()
        if not name or not symbol or not row['evidence']:
            raise ValueError(f'Missing alias name, symbol or evidence in {alias_table}')
        if name in aliases:
            raise ValueError(f'Duplicate alias-table name: {name}')
        if symbol not in exact_symbols:
            raise ValueError(f'Alias target absent from GENCODE v49: {symbol}')
        aliases[name] = {k: (v or '') for k, v in row.items()}
    audit_groups = defaultdict(set)
    mappings = defaultdict(list)
    rows = read_tsv(table)
    if not rows:
        raise ValueError(f'Empty motif table {table}')
    for row in rows:
        name = row.get('RBP', row.get('Protein_name', '')).strip()
        motif = row.get('Binding Motif', row.get('regularExpression', '')).strip()
        if not name or not motif:
            raise ValueError(f'Missing RBP/motif columns in {table}')
        if name not in aliases and name not in exact_symbols:
            raise ValueError(f'Name absent from HGNC alias table and GENCODE v49: {name}')
        alias = aliases.get(name)
        symbol = alias['hgnc_symbol'] if alias else name
        source = 'HGNC alias table' if alias else 'exact GENCODE gene_name'
        evidence = alias['evidence'] if alias else f'Exact gene_name in {gtf}'
        ambiguity = alias['ambiguity_note'] if alias else ''
        entry = {'_display_name': symbol, '_hgnc_symbol': symbol,
                 '_naming_action': 'changed_to_hgnc' if symbol != name else 'confirmed_hgnc',
                 '_exclusion_symbol': symbol, '_exclusion_basis': source}
        mappings[name].append((motif, entry))
        audit_groups[(name, symbol, source, evidence, ambiguity)].add(motif)
    for row in read_tsv(esrp_table):
        mappings[row['name']].append((row['motif'], {'_display_name': 'ESRP-like', '_hgnc_symbol': '',
                                                     '_naming_action': 'esrp_like_group', '_exclusion_symbol': '',
                                                     '_exclusion_basis': 'motif_group'}))
    collapsed = defaultdict(set)
    for name, symbol, _, _, _ in audit_groups:
        collapsed[symbol].add(name)
    merged = {k: sorted(v) for k, v in collapsed.items() if len(v) > 1}
    audits = [dict(zip(['table_name', 'hgnc_symbol', 'source', 'evidence', 'ambiguity_note'], key),
                   n_motifs=len(motifs), merged_into=key[1] if key[1] in merged else '')
              for key, motifs in sorted(audit_groups.items())]
    unused = sorted(set(aliases) - set(mappings))
    if unused:
        raise ValueError(f'Alias-table names absent from motif table: {unused}')
    stats = {'table_names_changed': sum(a['table_name'] != a['hgnc_symbol'] for a in audits),
             'rbp_groups_merged': len(merged), 'table_names_collapsed': sum(len(v) - 1 for v in merged.values()),
             'merged_groups': merged, 'mapping_complete': True, 'alias_table': str(Path(alias_table).resolve()),
             'alias_mappings': {n: r['hgnc_symbol'] for n, r in sorted(aliases.items())},
             'ambiguity_notes': {n: r['ambiguity_note'] for n, r in sorted(aliases.items()) if r['ambiguity_note']},
             'esrp_like_exception': 'Synthetic ESRP-like motif keys retain their existing group; not HGNC gene names '
                                    'and omitted from gene naming audit.'}
    return mappings, audits, stats


def apply_naming(rows, mappings):
    for row in rows:
        candidates = mappings.get(row['RBP'], [])
        if not candidates:
            raise ValueError(f'RBP missing from motif tables: {row["RBP"]}')
        motif = row['motif_key'].split('.', 1)[1]
        exact = [entry for seq, entry in candidates if seq == motif]
        entries = exact or [entry for _, entry in candidates]
        unique = {json.dumps(e, sort_keys=True): e for e in entries}
        if len(unique) != 1:
            raise ValueError(f'Ambiguous table identity for {row["RBP"]}: {motif}')
        row.update(next(iter(unique.values())))


def load_exclusion_list(path):
    symbols = set()
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        body = line.split('#', 1)[0].strip()
        if body:
            symbols.add(body.split()[0])
    if not symbols:
        raise ValueError(f'Empty exclusion list {path}')
    return symbols


def exclusion_audit(arm, entries, lists):
    grouped = defaultdict(list)
    for r in entries:
        grouped[r['RBP']].append(r)
    audits, dropped = [], set()
    for name, rows in sorted(grouped.items()):
        identities = {(r['_display_name'], r['_exclusion_symbol'], r['_exclusion_basis']) for r in rows}
        for display, symbol, basis in sorted(identities):
            matched = [k for k, v in lists.items() if symbol in v]
            remove = bool(matched) and symbol != 'SRSF1'
            if remove:
                dropped.add(symbol)
            audits.append({'arm': arm, 'table_name': name, 'display_name': display,
                           'hgnc_symbol': rows[0]['_hgnc_symbol'], 'exclusion_symbol': symbol,
                           'symbol_basis': basis, 'matched_lists': ';'.join(matched),
                           'dropped_by': ';'.join(matched) if remove else '',
                           'action': 'dropped' if remove else ('always_keep_SRSF1' if symbol == 'SRSF1' else 'retained'),
                           'n_motifs': len({r['motif_key'] for r in rows if r['_display_name'] == display})})
    return audits, dropped


def gate_counts(event_sets_root, arm, gate):
    if gate != 'persample10_bm50_bgfdr0.5':
        raise ValueError('v4 requires persample10_bm50_bgfdr0.5')
    base, rule = arm.rsplit('_', 1)
    path = Path(event_sets_root) / base / f'{gate}_rule{rule}' / 'counts.json'
    source = json.loads(path.read_text())
    counts = {k: int(source[k]) for k in ['n_up', 'n_dn', 'n_bg', 'n_expr_unknown_in_fg', 'n_expr_unknown_in_bg']}
    rows = [r for r in read_tsv(Path(event_sets_root) / 'gate_sensitivity.tsv')
            if r['arm'] == base and r['coverage'] == 'persample10' and r['floor'] == 'bm50'
            and float(r['bg_fdr']) == .5 and r['rule'] == rule]
    if len(rows) > 1 or (not rows and rule != 'Beffect'):
        raise ValueError(f'Expected one gate_sensitivity row for {arm}; found {len(rows)}')
    if rows and any(int(rows[0][k]) != v for k, v in counts.items()):
        raise ValueError(f'Gate count sources disagree for {arm}')
    return counts, path


# ---------------------------------------------------------------- layer readers
def released_entries(arm, root, mappings):
    """Authors' released rMAPS3 --stat-method mannwhitney root tables -> pooled-region entries."""
    sources, entries = [], []
    for code, direction in [('up', 'INCLUDED'), ('dn', 'SKIPPED')]:
        path = Path(root) / arm / f'pVal.{code}.vs.bg.RNAmap.txt'
        sources.append(path)
        with path.open(newline='') as f:
            rows = list(csv.reader(f, delimiter='\t'))
        header = rows[0]
        if len(header) != 9 or header[0] != 'RBP' or any(
                header[i + 1] != f'smallest_p_in_{name}' for i, name in enumerate(SUBREGIONS)):
            raise ValueError(f'Unexpected released rank-sum schema: {path}')
        for row in rows[1:]:
            if len(row) != 9:
                raise ValueError(f'Ragged released rank-sum row in {path}')
            key = row[0]
            for pooled, (i, j) in POOLING.items():
                best = min((i, j), key=lambda k: (float(row[k + 1]), SUBREGIONS[k]))
                entries.append({'motif_key': key, 'RBP': key.split('.', 1)[0], 'direction_label': direction,
                                'pooled_region': pooled, 'region': SUBREGIONS[best],
                                '_p': float(row[best + 1]), '_native_p': float(row[best + 1]), '_q': None,
                                '_ratio': 1.0, '_calib_perms_used': None, '_calib_stage': ''})
    log = Path(root) / f'{arm}_command.log'
    sources.append(log)
    text = log.read_text()
    if RELEASED_COMMIT not in text or f'stat={STAT_METHOD}' not in text or 'exit=0' not in text:
        raise ValueError(f'Released rank-sum command log does not record commit/stat/exit: {log}')
    n_motifs = len({e['motif_key'] for e in entries})
    if len(entries) != n_motifs * 6:
        raise ValueError(f'Released rank-sum entry count {len(entries)} != {n_motifs} motifs x 6 panels')
    apply_naming(entries, mappings)
    for e in entries:
        e['_label'] = rbp_label(e)
        if not 0 < e['_p'] <= 1:
            raise ValueError(f'Invalid released rank-sum p for {e["motif_key"]}: {e["_p"]}')
    return entries, sources, n_motifs


def calibrated_entries(arm, root, mappings):
    """Westfall-Young calibrated rank-sum summary -> pooled-region entries (one per motif x panel)."""
    base = Path(root) / arm
    table = base / 'per_motif_regions.tsv'
    report = json.loads((base / 'refinement_report.json').read_text())
    sources = [table, base / 'condensed_per_rbp.tsv', base / 'refinement_report.json',
               base / 'versions.txt', base / 'command.log', base / 'input_md5s.tsv', base / 'readout.md']
    log = (base / 'command.log').read_text()
    if 'exit=0' not in log:
        raise ValueError(f'Calibration run for {arm} did not record exit=0; refusing to plot a partial summary')
    rows = [r for r in read_tsv(table) if r['plot'].strip().upper() == 'TRUE']
    grouped = defaultdict(list)
    for r in rows:
        grouped[(r['direction_label'], r['pooled_region'], r['motif_key'])].append(r)
    entries = []
    for (direction, pooled, key), group in grouped.items():
        if pooled not in REGIONS or direction not in DIRECTIONS:
            raise ValueError(f'Unexpected calibrated panel {direction} x {pooled}')
        if len({r['calibrated_q'] for r in group}) != 1 or len({r['calibrated_p_pooled'] for r in group}) != 1:
            raise ValueError(f'Calibrated pooled p/q disagree within {key} {direction} {pooled}')
        rep = min(group, key=lambda r: (float(r['native_ranksum_p']), -float(r['enrichment_ratio']), r['region']))
        if not math.isclose(min(float(r['native_ranksum_p']) for r in group),
                            float(rep['native_ranksum_p_pooled']), rel_tol=1e-12):
            raise ValueError(f'Calibrated native pooled p is not the sub-region minimum: {key}')
        entries.append({'motif_key': key, 'RBP': rep['rbp_table_name'], '_calib_rbp': rep['RBP'],
                        'direction_label': direction, 'pooled_region': pooled, 'region': rep['region'],
                        '_p': float(rep['calibrated_p_pooled']), '_native_p': float(rep['native_ranksum_p_pooled']),
                        '_q': float(rep['calibrated_q']), '_native_q': float(rep['native_q']),
                        '_ratio': float(rep['enrichment_ratio']),
                        '_calib_perms_used': int(rep['calib_perms_used_pooled']),
                        '_calib_stage': rep['calib_stage_pooled']})
    n_motifs = len({e['motif_key'] for e in entries})
    if len(entries) != n_motifs * 6 or len(entries) != report['bh_family_size']:
        raise ValueError(f'Calibrated entry count {len(entries)} disagrees with BH family {report["bh_family_size"]}')
    apply_naming(entries, mappings)
    for e in entries:
        e['_label'] = rbp_label(e)
        expected = e['_calib_rbp']
        ours = e['_display_name']
        if re.fullmatch('motif_[0-9]+', e['RBP']):
            if expected != e['RBP']:
                raise ValueError(f'Calibrated ESRP-like name changed: {e["RBP"]} vs {expected}')
        elif ours != expected:
            raise ValueError(f'HGNC naming disagrees with the calibration summary: {ours} vs {expected}')
        if not 0 < e['_p'] <= 1 or not 0 <= e['_q'] <= 1:
            raise ValueError(f'Invalid calibrated p/q for {e["motif_key"]}')
        if e['_ratio'] < 0 or not math.isfinite(e['_ratio']):
            raise ValueError(f'Invalid enrichment ratio for {e["motif_key"]}')
    return entries, sources, report


def cross_check_layers(released, calibrated, arm):
    """The calibrated layer recalibrates the released statistic: native pooled p must be identical."""
    rel = {(e['direction_label'], e['pooled_region'], e['motif_key']): e['_native_p'] for e in released}
    cal = {(e['direction_label'], e['pooled_region'], e['motif_key']): e['_native_p'] for e in calibrated}
    if set(rel) != set(cal):
        raise ValueError(f'Released and calibrated rank-sum motif universes differ for {arm}')
    bad = [k for k in rel if not math.isclose(rel[k], cal[k], rel_tol=1e-9)]
    if bad:
        raise ValueError(f'Native rank-sum p differs between released table and calibration summary: {arm} {bad[:3]}')
    return len(rel)


# ---------------------------------------------------------------- selection + drawing
def ratio_area(r):
    if not math.isfinite(r) or r < 0:
        raise ValueError(f'Invalid enrichment ratio {r}')
    return 24 + 116 * math.log2(1 + r)


def qcolor(q):
    if not math.isfinite(q):
        raise ValueError('Missing q')
    return COLORS[0 if q < .01 else 1 if q < .05 else 2 if q < .1 else 3]


def pcolor(p):
    return COLORS[0 if p < .001 else 1 if p < .01 else 2 if p < .05 else 3]


def rank_key(r):
    return (r['_p'], -r['_ratio'], r['_label'], r['motif_key'])


def select_panels(entries, layer, kind, top_n=10):
    panels = {(d, p): [] for d in DIRECTIONS for p in REGIONS}
    for e in entries:
        panels[(e['direction_label'], e['pooled_region'])].append(dict(e))
    for key, rows in panels.items():
        if kind == 'byRBP':
            groups = defaultdict(list)
            for r in rows:
                groups[r['_label']].append(r)
            rows = []
            for group in groups.values():
                r = dict(min(group, key=rank_key))
                r['_n'] = len(group)
                r['_k'] = sum((x['_p'] if layer.endswith('rawP') else x['_q']) < .05 for x in group)
                rows.append(r)
        else:
            for r in rows:
                r['_n'] = 1
                r['_k'] = int((r['_p'] if layer.endswith('rawP') else r['_q']) < .05)
        for r in rows:
            r['_tick'] = tick_label(r, kind)
        panels[key] = sorted(rows, key=rank_key)[:top_n]
    return panels


def y_scale(entries, layer):
    """One shared y-scale per arm x layer, used by byRBP/byMotif and both exclusion variants.

    The cap is data-driven: 1.25 x the largest -log10 p among RBPs other than the single top RBP of
    the arm, rounded up to a clean tick. Stems above the cap are truncated and labelled with their
    value. When the cap already exceeds the largest value nothing is truncated and the scale is the
    plain rounded maximum, so the rule is uniform but never adds a break where none is needed.
    """
    values = [(-math.log10(e['_p']), e['_label']) for e in entries]
    layer_max, top_rbp = max(values)
    others = [v for v, label in values if label != top_rbp]
    if not others:
        raise ValueError('Cannot derive a y cap: every entry belongs to one RBP')
    cap_base = max(others)
    cap = 1.25 * cap_base
    scale_value = layer_max if cap >= layer_max else cap
    if layer == 'calibrated_ranksum':
        ymax = float(max(4, math.ceil(scale_value)))
        step = 1.0
    else:
        raw_step = scale_value / 5
        magnitude = 10 ** math.floor(math.log10(raw_step)) if raw_step > 0 else 1
        step = next(v * magnitude for v in (1, 2, 2.5, 5, 10) if v * magnitude >= raw_step)
        ymax = max(step, math.ceil(scale_value / step) * step)
    return {'ymax': ymax, 'step': step, 'layer_max': layer_max, 'top_rbp': top_rbp,
            'cap_base': cap_base, 'cap': cap, 'cap_applied': bool(cap < layer_max),
            'n_entries_above_cap': sum(v > ymax + 1e-12 for v, _ in values)}


def truncation_text(value):
    return f'−log10 p = {value:.1f}'


def subtitle(layer, refinement):
    if layer == 'released_ranksum_rawP':
        return ("Authors' rMAPS3 (released code " + RELEASED_COMMIT + ", --stat-method " + STAT_METHOD + "): one-sided "
                "rank-sum on per-exon motif hit counts,\nsmallest p over the 50-nt windows of each region; raw p, "
                "no multiple-testing adjustment (none is part of the tool).\n"
                "p comes from the tie-corrected normal approximation; the RBP ORDER is what this panel claims.")
    stage1 = f"{refinement['stage1_permutations']:,}"
    stage2 = f"{refinement['stage2_permutations']:,}"
    floor = f"1/{refinement['stage2_permutations'] + 1:,}"
    return (f"Same statistic, Westfall–Young permutation-calibrated p: changed/background labels permuted, min-P over "
            f"the windows of a region.\nB = {stage1} permutations, refined to {stage2} where stage-1 p ≤ "
            f"{refinement['refine_threshold']}; floor {floor}; BH q over {refinement['bh_family_size']} tests "
            f"(126 motifs × 2 directions × 3 regions).\n"
            "Supplement to the main figure: how much of the released rank-sum signal survives a permutation null.")


def draw_legend(fig, layer, refinement, excluded=False):
    raw = layer == 'released_ranksum_rawP'
    if raw:
        lines = ['Stem height = −log10 raw rank-sum p',
                 "Layer: authors' released rMAPS3 " + RELEASED_COMMIT + ' with --stat-method ' + STAT_METHOD +
                 ' — one-sided rank-sum on per-exon motif hit counts, regional minimum over 50-nt windows, raw p',
                 'Dot size is constant: the released tool reports no enrichment ratio.',
                 'Dot colour = raw rank-sum p (no multiple-testing adjustment)', None,
                 'included = exons MORE included in treatment; skipped = exons LESS included',
                 'stems above the cap are truncated and labelled with their value']
    else:
        lines = ['Stem height = −log10 calibrated p',
                 'Layer: the same rank-sum statistic with a Westfall–Young label-permutation p (B = '
                 f'{refinement["stage1_permutations"]:,}, refined to {refinement["stage2_permutations"]:,}),',
                 f'BH q over {refinement["bh_family_size"]} tests — the supplement, not the main claim',
                 'Dot size = enrichment ratio: % of changed exons carrying the motif ÷ % of background exons carrying '
                 'it, at the most significant 50-nt window', 'SIZE_KEY',
                 'Dot colour = calibrated BH q', None,
                 'included = exons MORE included in treatment; skipped = exons LESS included',
                 'stems above the cap are truncated and labelled with their value']
    if excluded:
        lines += ['core spliceosome + broad binders excluded (lab RBP-RELI lists); SRSF1 always shown']
    height_pt = (len(lines) - 1) * 26 + 32
    ax = fig.add_axes([.055, .198 - height_pt / (fig.get_figheight() * 72), .70, height_pt / (fig.get_figheight() * 72)])
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, height_pt)
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor='#FAFAFA', edgecolor='#AAAAAA',
                           lw=.7, clip_on=False))
    for i, line in enumerate(lines):
        y = height_pt - 16 - i * 26
        if line == 'SIZE_KEY':
            for x, r, label in zip([.025, .26, .44], [1, 2, 4], ['1× (no enrichment)', '2×', '4×']):
                ax.scatter([x], [y], s=ratio_area(r), marker='o', color='#BBBBBB', edgecolor='#333333', linewidth=.7)
                ax.text(x + .022, y, label, va='center', fontsize=14)
        elif line is None:
            labels = (['p < 0.001', '0.001–0.01', '0.01–0.05', '≥0.05 (not significant)'] if raw else
                      ['q < 0.01', '0.01–0.05', '0.05–0.10', '≥0.10 (not significant)'])
            for x, c, label in zip([.025, .26, .44, .62], COLORS, labels):
                ax.scatter([x], [y], s=140, marker='o', color=c, edgecolor='#333333', linewidth=.7)
                ax.text(x + .022, y, label, fontsize=14, va='center')
        else:
            ax.text(.018, y, line, fontsize=14, va='center')


def polygons_overlap(a, b, tolerance=0.5):
    for poly in (a, b):
        for i in range(len(poly)):
            x1, y1 = poly[i]
            x2, y2 = poly[(i + 1) % len(poly)]
            nx, ny = -(y2 - y1), x2 - x1
            norm = math.hypot(nx, ny)
            if norm == 0:
                continue
            nx /= norm
            ny /= norm
            pa = [x * nx + y * ny for x, y in a]
            pb = [x * nx + y * ny for x, y in b]
            if min(max(pa), max(pb)) - max(min(pa), min(pb)) <= tolerance:
                return False
    return True


def box_polygon(b):
    return [[b.x0, b.y0], [b.x1, b.y0], [b.x1, b.y1], [b.x0, b.y1]]


def layout_audit(fig, axes, model_artists, path, ymax, max_value):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    w, h = fig.canvas.get_width_height()
    texts, objects = [], []
    for t in fig.texts:
        objects.append((t, 'figure'))
    for ax in fig.axes:
        objects.extend((t, 'annotation') for t in ax.texts)
        objects.extend((t, 'title') for t in (ax.title, ax._left_title, ax._right_title))
        if ax.axison:
            for axis in (ax.xaxis, ax.yaxis):
                objects.append((axis.label, 'axis_label'))
                objects.append((axis.offsetText, 'offset'))
                for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                    for t in (tick.label1, tick.label2):
                        objects.append((t, 'tick'))
    seen = set()
    for t, role in objects:
        if id(t) in seen or not t.get_visible() or not t.get_text():
            continue
        seen.add(id(t))
        bb = t.get_window_extent(renderer)
        if not bb.width or not bb.height:
            continue
        angle = t.get_rotation()
        t.set_rotation(0)
        flat = t.get_window_extent(renderer)
        t.set_rotation(angle)
        theta = math.radians(angle)
        c, ss = math.cos(theta), math.sin(theta)
        cx, cy = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
        poly = [[cx + x * c - y * ss, cy + x * ss + y * c] for x, y in
                [(-flat.width / 2, -flat.height / 2), (flat.width / 2, -flat.height / 2),
                 (flat.width / 2, flat.height / 2), (-flat.width / 2, flat.height / 2)]]
        texts.append({'text': t.get_text(), 'role': role, 'bbox': list(bb.extents), 'polygon': poly,
                      'inside': bool(bb.x0 >= 0 and bb.y0 >= 0 and bb.x1 <= w and bb.y1 <= h)})
    overlaps = []
    for i, t in enumerate(texts):
        for j in range(i):
            if polygons_overlap(t['polygon'], texts[j]['polygon']):
                overlaps.append({'first': j, 'second': i, 'texts': [texts[j]['text'], t['text']]})
    axes_hits, model_hits = [], []
    for i, t in enumerate(texts):
        if t['text'] == 'permutation floor' or t['text'].startswith('−log10 p = '):
            continue
        for ai, ax in enumerate(axes):
            if polygons_overlap(t['polygon'], box_polygon(ax.get_window_extent(renderer))):
                axes_hits.append({'text_index': i, 'text': t['text'], 'panel': ai})
        for mi, artist in enumerate(model_artists):
            if polygons_overlap(t['polygon'], box_polygon(artist.get_window_extent(renderer))):
                model_hits.append({'text_index': i, 'text': t['text'], 'model_part': mi})
    limits = [list(ax.get_ylim()) for ax in axes]
    shared = all(sorted(v) == [0.0, ymax] for v in limits)
    mirrored = all(v[0] == 0 for v in limits[:3]) and all(v[1] == 0 for v in limits[3:])
    return {'figure': str(path), 'canvas': [w, h], 'all_inside': all(t['inside'] for t in texts),
            'bounds': texts, 'text_overlaps': overlaps, 'text_axes_overlaps': axes_hits,
            'text_gene_model_overlaps': model_hits, 'shared_y_limit': ymax, 'plotted_max_log10p': max_value,
            'panel_y_limits': limits, 'shared_y_limit_verified': shared, 'mirrored_axes_verified': mirrored,
            'overlap_method': 'Oriented text rectangles with separating-axis intersection; 0.5 pixel tolerance',
            'intentional_exclusions': ['permutation-floor text inside its data panel',
                                       'truncated-stem value labels inside their data panel'],
            'pass': all(t['inside'] for t in texts) and not (overlaps or axes_hits or model_hits)
                    and shared and mirrored}


def draw_break(ax, x, ymax):
    """Standard axis-break glyph on a truncated stem: a white gap crossed by two short diagonals."""
    low, high = .90 * ymax, .94 * ymax
    ax.plot([x, x], [low, high], color='white', lw=3.0, zorder=2.4, clip_on=False,
            solid_capstyle='butt')
    for y in (low, high):
        ax.plot([x - .28, x + .28], [y - .014 * ymax, y + .014 * ymax], color='#333333', lw=1.1,
                zorder=2.6, clip_on=False, solid_capstyle='butt')


def draw_figure(arm, layer, kind, panels, counts, refinement, scale, path, args, excluded=False):
    fig = plt.figure(figsize=(32, 28), dpi=PNG_DPI)
    fig.text(.5, .978, LABELS[arm] + ' — skipped-exon motif map', ha='center', va='top', fontsize=27, weight='bold')
    entries = [r for rows in panels.values() for r in rows]
    fig.text(.5, .946, subtitle(layer, refinement), ha='center', va='top', fontsize=15)
    if kind == 'byMotif':
        fig.text(.5, .908, "each dot is one motif; the motif sequences are in the arm's workbook",
                 ha='center', va='top', fontsize=15)
    fig.text(.5, .884, 'included', ha='center', color='#E69F00', fontsize=23, weight='bold')
    calibrated = layer == 'calibrated_ranksum'
    floor = -math.log10(1 / (refinement['stage2_permutations'] + 1)) if calibrated else None
    show_floor = calibrated and any(math.isclose(r['_p'], 1 / (r['_calib_perms_used'] + 1), rel_tol=1e-10, abs_tol=0)
                                    and r['_calib_perms_used'] == refinement['stage2_permutations'] for r in entries)
    max_value = max(-math.log10(r['_p']) for rows in panels.values() for r in rows)
    ymax, step = scale['ymax'], scale['step']
    truncated_rows = []
    axes = []
    for ri, d in enumerate(DIRECTIONS):
        for ci, p in enumerate(REGIONS):
            ax = fig.add_axes([.065 + ci * .305, .68 if ri == 0 else .31, .265, .175])
            axes.append(ax)
            rows = panels[(d, p)]
            for i, r in enumerate(rows):
                value = -math.log10(r['_p'])
                truncated = value > ymax + 1e-12
                y = ymax if truncated else value
                ax.vlines(i, 0, y, color='#BBBBBB', lw=1.2, zorder=2)
                area = 140 if not calibrated else ratio_area(r['_ratio'])
                dot = ax.scatter(i, y, s=area, color=qcolor(r['_q']) if calibrated else pcolor(r['_p']),
                                 edgecolor='#333333', linewidth=.7, zorder=3, clip_on=False)
                dot.set_gid(f'dot_{ri}_{ci}_{i}')
                if truncated:
                    draw_break(ax, i, ymax)
                    label = ax.text(i + .38, ymax, truncation_text(value), fontsize=12, color='#333333',
                                    ha='left', va='center', clip_on=False, zorder=5)
                    label.set_gid(f'truncation_label_{ri}_{ci}_{i}')
                    truncated_rows.append({'direction': d, 'pooled_region': p, 'rank': i + 1,
                                           'label': r['_label'], 'tick': r['_tick'],
                                           'motif_key': r['motif_key'], 'log10p': value,
                                           'printed': truncation_text(value)})
            ax.set_xlim(-1.1, args.top_n + .4)
            ax.set_ylim((0, ymax) if ri == 0 else (ymax, 0))
            ax.set_xticks(range(len(rows)), [r['_tick'] for r in rows], rotation=55,
                          ha='right' if ri == 0 else 'left', rotation_mode='anchor',
                          fontsize=13 if kind == 'byMotif' else 16)
            if ri == 1:
                ax.xaxis.tick_top()
                ax.tick_params(axis='x', top=True, labeltop=True, bottom=False, labelbottom=False, pad=7)
                for t in ax.get_xticklabels():
                    t.set_va('bottom')
            ax.set_yticks([i * step for i in range(round(ymax / step) + 1)])
            ax.tick_params(axis='y', labelsize=14)
            ax.set_ylabel('−log10 p', fontsize=16)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(ri == 1)
            ax.spines['bottom'].set_visible(ri == 0)
            ax.yaxis.grid(True, color='#EEEEEE', lw=.6)
            ax.set_axisbelow(True)
            if ri == 0:
                ax.set_title(p, fontsize=21, pad=14)
            if show_floor:
                ax.axhline(floor, color='#777777', lw=.6, ls=':', zorder=1)
                ax.text(.985, floor + (.015 * ymax if ri == 0 else -.015 * ymax), 'permutation floor',
                        transform=ax.get_yaxis_transform(), ha='right', va='bottom' if ri == 0 else 'top',
                        fontsize=10, color='#666666')
    model = fig.add_axes([0, .576, 1, .018])
    model.axis('off')
    model.set_xlim(0, 1)
    model.set_ylim(0, 1)
    model_artists = []
    for ci in (0, 2):
        patch = Rectangle((.065 + ci * .305, .46), .265, .08, facecolor='black')
        model.add_patch(patch)
        model_artists.append(patch)
    exon = Rectangle((.370, .12), .265, .76, facecolor='#999999', edgecolor='black', lw=1)
    model.add_patch(exon)
    model_artists.append(exon)
    model.text(.052, .5, "5'", ha='right', va='center', weight='bold', fontsize=19)
    model.text(.953, .5, "3'", ha='left', va='center', weight='bold', fontsize=19)
    fig.text(.5, .280, 'skipped', ha='center', color='#0072B2', fontsize=23, weight='bold')
    rule = arm.rsplit('_', 1)[1]
    fig.text(.5, .249, 'Gate: per-sample junction reads ≥10 (no unknown category); gene baseMean >50 '
             f'(genes absent from the DESeq2 table are kept: {counts["n_expr_unknown_in_fg"]:,} foreground / '
             f'{counts["n_expr_unknown_in_bg"]:,} background events)', ha='center', fontsize=14, color='#555555')
    fig.text(.5, .225, f'Foreground FDR <0.05 and |ΔPSI| ≥0.10; background FDR ≥0.5; rule {rule}  |  '
             f'n included={counts["n_up"]:,} / n skipped={counts["n_dn"]:,} / n background={counts["n_bg"]:,}',
             ha='center', fontsize=14, color='#555555')
    if rule == 'Beffect':
        fig.text(.5, .208, 'rule B (relaxed): frozen ±5 bp join, VAST coverage VLOW+, VAST |dPSI| ≥ 0.10, '
                 'same direction; no MV requirement', ha='center', fontsize=14, color='#555555')
    draw_legend(fig, layer, refinement, excluded=excluded)
    tail = ('Released engine commit ' + RELEASED_COMMIT + ', --stat-method ' + STAT_METHOD +
            '; input event sets are md5-identical to every other engine run of this arm.' if not calibrated else
            f'Calibration seed {refinement["seed"]} on the same released engine run ({RELEASED_COMMIT}, '
            f'--stat-method {STAT_METHOD}); no p-value recomputed for the main figure.')
    fig.text(.5, .055, tail, ha='center', fontsize=14, color='#555555')
    report = layout_audit(fig, axes, model_artists, path, ymax, max_value)
    report.update({'layer': layer, 'kind': kind, 'arm': arm, 'excluded': excluded,
                   'refined_entries': sum(r['_calib_perms_used'] > refinement['stage1_permutations']
                                          for r in entries) if calibrated else None,
                   'refined_floor_drawn': show_floor, 'y_scale': scale,
                   'n_stems_truncated': len(truncated_rows), 'truncated_stems': truncated_rows})
    fig.savefig(str(path) + '.svg', metadata={'Date': None})
    fig.savefig(str(path) + '.png', dpi=PNG_DPI)
    plt.close(fig)
    if not report['pass']:
        Path(str(path) + '_layout_failure.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        raise ValueError(f'Layout audit failed: {path}; see _layout_failure.json')
    return report


# ---------------------------------------------------------------- rank comparison workbook
V31_LAYERS = ['released_rawP', 'audited_rawP', 'audited_BH', 'audited_calibrated']


def _sheet(book, name, fields, rows, freeze='A2'):
    sheet = book.create_sheet(name)
    sheet.append(list(fields))
    for row in rows:
        sheet.append([row.get(f) if row.get(f) is not None else 'NA' for f in fields])
    for cell in sheet[1]:
        cell.font = Font(name='Arial', bold=True, color='FFFFFF')
        cell.fill = PatternFill('solid', fgColor='0072B2')
        cell.alignment = Alignment(wrap_text=True, vertical='top')
    for cells in sheet.iter_rows(min_row=2):
        for cell in cells:
            cell.font = Font(name='Arial', size=11)
            cell.alignment = Alignment(vertical='top')
            if isinstance(cell.value, float):
                cell.number_format = '0.00000E+00' if 0 < abs(cell.value) < .001 else '0.0000'
    sheet.freeze_panes = freeze
    sheet.auto_filter.ref = sheet.dimensions
    sheet.row_dimensions[1].height = 42
    for index, field in enumerate(fields, 1):
        sheet.column_dimensions[openpyxl.utils.get_column_letter(index)].width = min(42, max(17, len(field) + 2))
    return sheet


def build_rank_comparison(arm, layer_panels, v31_tsv, method_tsv, refinement, dest):
    """Ranks for every available v4 layer, joined to the archived v3.1 layer ranks."""
    v4_layers = [l for l in LAYERS if l in layer_panels]
    panels = sorted({key for l in v4_layers for key in layer_panels[l]})
    prior = {}
    prior_layers = []
    if v31_tsv and Path(v31_tsv).is_file():
        for row in read_tsv(v31_tsv):
            prior[(row['direction'], row['pooled_region'], row['RBP'])] = row
        prior_layers = [l for l in V31_LAYERS if f'{l}_rank' in next(iter(prior.values()))]
    all_layers = v4_layers + prior_layers
    fields = ['arm', 'direction', 'pooled_region', 'RBP']
    for l in v4_layers:
        fields += [f'{l}_p', f'{l}_rank', f'{l}_motif', f'{l}_top10']
        if l == 'calibrated_ranksum':
            fields += ['calibrated_ranksum_q', 'calibrated_ranksum_native_p', 'calib_perms_used', 'calib_stage']
    for l in prior_layers:
        fields += [f'v3.1_{l}_p', f'v3.1_{l}_rank']
    ranks, comparisons, universe_notes = [], [], []
    for direction, region in panels:
        key = (direction, region)
        rank_maps, top_sets, selected = {}, {}, {}
        for l in v4_layers:
            rows = sorted(layer_panels[l][key], key=rank_key)
            selected[l] = {r['_label']: r for r in rows}
            if len(selected[l]) != len(rows):
                raise ValueError(f'Duplicate HGNC group {arm} {l} {key}')
            rank_maps[l] = dict(zip([r['_label'] for r in rows], rankdata([r['_p'] for r in rows], method='average')))
            top_sets[l] = {r['_label'] for r in rows[:10]}
        names = set(selected[v4_layers[0]])
        if any(set(selected[l]) != names for l in v4_layers):
            raise ValueError(f'RBP universes differ between v4 layers for {arm} {key}')
        prior_names = {n for n in names if (direction, region, n) in prior}
        if prior_layers and prior_names != names:
            universe_notes.append({'arm': arm, 'direction': direction, 'pooled_region': region,
                                   'v4_only': ';'.join(sorted(names - prior_names)) or 'none',
                                   'n_shared': len(prior_names)})
        for l in prior_layers:
            rank_maps[l] = {n: float(prior[(direction, region, n)][f'{l}_rank']) for n in prior_names}
            top_sets[l] = {n for n in prior_names if prior[(direction, region, n)][f'{l}_top10'] == '1'}
        for name in sorted(names):
            record = dict(arm=arm, direction=direction, pooled_region=region, RBP=name)
            for l in v4_layers:
                r = selected[l][name]
                record.update({f'{l}_p': float(r['_p']), f'{l}_rank': float(rank_maps[l][name]),
                               f'{l}_motif': r['motif_key'], f'{l}_top10': int(name in top_sets[l])})
                if l == 'calibrated_ranksum':
                    record.update({'calibrated_ranksum_q': float(r['_q']),
                                   'calibrated_ranksum_native_p': float(r['_native_p']),
                                   'calib_perms_used': r['_calib_perms_used'], 'calib_stage': r['_calib_stage']})
            for l in prior_layers:
                row = prior.get((direction, region, name))
                record[f'v3.1_{l}_p'] = float(row[f'{l}_p']) if row else None
                record[f'v3.1_{l}_rank'] = float(row[f'{l}_rank']) if row else None
            ranks.append(record)
        for first, second in itertools.combinations(all_layers, 2):
            shared = sorted(set(rank_maps[first]) & set(rank_maps[second]))
            x = np.array([rank_maps[first][n] for n in shared])
            y = np.array([rank_maps[second][n] for n in shared])
            defined = len(x) >= 2 and np.ptp(x) > 0 and np.ptp(y) > 0
            comparisons.append(dict(arm=arm, direction=direction, pooled_region=region, layer1=first, layer2=second,
                                    n_rbps=len(shared), rho=float(np.corrcoef(x, y)[0, 1]) if defined else None,
                                    rho_status='defined' if defined else 'undefined_constant_or_single_rank',
                                    top10_overlap=len(top_sets[first] & top_sets[second]),
                                    top10_n_layer1=len(top_sets[first]), top10_n_layer2=len(top_sets[second])))
    summaries = []
    for first, second in itertools.combinations(all_layers, 2):
        rows = [r for r in comparisons if r['layer1'] == first and r['layer2'] == second]
        valid = [r['rho'] for r in rows if r['rho'] is not None]
        summaries.append(dict(arm=arm, layer1=first, layer2=second,
                              mean_rho=float(np.mean(valid)) if valid else None,
                              mean_top10_overlap=float(np.mean([r['top10_overlap'] for r in rows])),
                              n_panels=len(rows), n_defined_rho_panels=len(valid)))
    write_tsv(dest / f'{arm}_rank_comparison_v4.tsv', fields, ranks)
    write_tsv(dest / f'{arm}_rank_panel_comparisons_v4.tsv', list(comparisons[0]), comparisons)
    method_rows = []
    if method_tsv and Path(method_tsv).is_file():
        wide = defaultdict(dict)
        for row in read_tsv(method_tsv):
            if row['arm'] != arm:
                continue
            k = (row['pooled_region'], row['direction_label'], row['rbp'])
            wide[k][row['method'] + '_rank'] = float(row['rank'])
            wide[k][row['method'] + '_p'] = float(row['p_value'])
            wide[k][row['method'] + '_motif'] = row['best_motif']
        for (region, direction, rbp), values in sorted(wide.items()):
            method_rows.append(dict(arm=arm, pooled_region=region, direction_label=direction, rbp=rbp, **values))
    readme = [
        ('Purpose', f'{arm}: RBP ranks per pooled region × direction under every rank-sum layer drawn in figure '
                    'version 4, joined to the four layers of figure version 3.1 that were archived under '
                    '_v3.1_archive_2026-09-17. No p-value or q-value is recomputed here.'),
        ('Sheet: Ranks', 'One row per panel × HGNC-grouped RBP. <layer>_p is the plotted statistic, <layer>_rank is '
                         'the average rank of that statistic across the whole panel, <layer>_motif is the motif that '
                         'represented the RBP, <layer>_top10 flags membership of the deterministic top 10.'),
        ('Sheet: Panel comparisons', 'Spearman rho of the paired average-rank vectors and top-10 overlap for every '
                                     'layer pair inside one panel; rho is NA when a rank vector is constant.'),
        ('Sheet: Summary', 'Arithmetic mean of rho and of top-10 overlap over the six panels; '
                           'n_defined_rho_panels is the denominator for the mean rho.'),
        ('Sheet: Method comparison', 'Ranks and p-values for the six statistical methods computed independently at '
                                     'the table given by --method-comparison, for this arm only. '
                                     'Empty when the arm was not part of that comparison run.'),
        ('Sheet: Methods', 'Definitions, tie policy and the reason the rank-sum p-values run far below the Fisher '
                           'p-values on the same data.'),
        ('MAIN vs SUPPLEMENT', 'released_ranksum_rawP is the main figure layer. calibrated_ranksum is the supplement. '
                               'The v3.1_* columns are the Fisher-based layers, kept as comparison columns only.'),
        ('Scope', 'Homo sapiens / GRCh38 (hg38) / GENCODE v49; skipped-exon (SE) events only; main variant '
                  '(no exclusion list applied); all HGNC-grouped RBPs.'),
    ]
    methods = [
        ('Layer 1 statistic', f"Released rMAPS3 at commit {RELEASED_COMMIT} run with --stat-method {STAT_METHOD}: a "
                              'one-sided Mann-Whitney rank-sum test on per-exon motif hit counts, taking the smallest '
                              'p over the 50-nt windows of a region. Raw p; the tool applies no adjustment.'),
        ('Layer 2 statistic', f'The same statistic with a Westfall-Young min-P label-permutation p '
                              f'(B = {refinement["stage1_permutations"]}, refined to '
                              f'{refinement["stage2_permutations"]} where stage-1 p <= {refinement["refine_threshold"]}'
                              f', seed {refinement["seed"]}) and BH q over {refinement["bh_family_size"]} tests.'
                              if refinement else 'Not built for this arm; the calibration summary was absent.'),
        ('Why the rank-sum p is so small', 'The ranked unit is one observation per eligible exon per window. With more '
                                           'than 99% of exons carrying no motif hit, the tie correction collapses the '
                                           'variance of the normal approximation, so the approximate p runs many '
                                           'orders of magnitude below an exact test on the same 2x2 table. The '
                                           'p-values are anti-conservative in the sparse tail; the RBP ORDER is the '
                                           'claim; see the statistical-method comparison report.'),
        ('Pooled regions', 'Upstream Intron = min(upstreamExonIntron, upstreamIntron); Exon Body = '
                           'min(targetExon-5prime, targetExon-3prime); Downstream Intron = min(downstreamIntron, '
                           'downstreamExonIntron). The two flanking-exon sub-regions are not plotted.'),
        ('Selection', 'One best motif per RBP per layer, panel and direction: smallest layer p, then enrichment ratio '
                      'descending, then HGNC display label, then motif key. Released ratios do not exist and are held '
                      'constant for tie ordering.'),
        ('Tied ranks', 'scipy.stats.rankdata(method="average").'),
        ('Grouping', 'Table names are mapped to HGNC symbols with the supplied alias table and verified against '
                     'GENCODE v49 gene_name. The twelve synthetic ESRP-like hexamers form one ESRP-like group.'),
        ('v3.1 join', 'The v3.1 ranks are read verbatim from the archived per-arm rank table and joined on '
                      '(direction, pooled region, RBP). Nothing in v3.1 is recomputed.'),
        ('Export', 'Immutable numerical export, not a calculation model.'),
    ]
    book = openpyxl.Workbook()
    book.remove(book.active)
    sheet = _sheet(book, 'README', ('topic', 'description'), [dict(topic=k, description=v) for k, v in readme])
    sheet.column_dimensions['A'].width = 27
    sheet.column_dimensions['B'].width = 115
    for row in sheet.iter_rows(min_row=2):
        row[1].alignment = Alignment(wrap_text=True, vertical='top')
        sheet.row_dimensions[row[0].row].height = 58
    _sheet(book, 'Summary', ('arm', 'layer1', 'layer2', 'mean_rho', 'mean_top10_overlap', 'n_panels',
                             'n_defined_rho_panels'), summaries)
    _sheet(book, 'Ranks', fields, ranks, freeze='D2')
    _sheet(book, 'Panel comparisons', list(comparisons[0]), comparisons)
    if method_rows:
        _sheet(book, 'Method comparison', list(method_rows[0]), method_rows, freeze='E2')
    sheet = _sheet(book, 'Methods', ('topic', 'definition'), [dict(topic=k, definition=v) for k, v in methods])
    sheet.column_dimensions['A'].width = 27
    sheet.column_dimensions['B'].width = 115
    for row in sheet.iter_rows(min_row=2):
        row[1].alignment = Alignment(wrap_text=True, vertical='top')
        sheet.row_dimensions[row[0].row].height = 58
    book.save(dest / f'{arm}_rank_comparison.xlsx')
    return summaries, universe_notes, bool(method_rows)


# ---------------------------------------------------------------- index
def write_index(out, records, skipped, archive_name):
    parts = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>rMAPS3 rank-sum motif maps (v4)</title>'
             '<style>body{font:16px Arial;margin:30px;line-height:1.5;max-width:1500px}'
             '.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}img{width:100%}'
             'section{margin:50px 0}article{border:1px solid #ddd;padding:8px}'
             'table{border-collapse:collapse;margin:12px 0}td,th{border:1px solid #ccc;padding:5px 10px;'
             'text-align:left}.main{background:#EAF3FA}.supp{background:#F6F6F6}'
             '@media(max-width:900px){.grid{grid-template-columns:1fr}}</style>',
             '<h1>Skipped-exon motif maps — rank-sum layers (version 4)</h1>',
             '<p>Human / GRCh38 (hg38) / GENCODE v49. Skipped-exon (SE) events only. '
             'Every figure is a 600-equivalent raster at ' + str(PNG_DPI) + ' dpi on a 32 × 28 in canvas plus an '
             'editable-text Arial SVG.</p>',
             '<h2>How to read</h2>',
             '<p class="main"><b>MAIN — Authors\u2019 rMAPS3 rank-sum, raw p.</b> The collaborators\u2019 released '
             'rMAPS3 code at commit ' + RELEASED_COMMIT + ' run with <code>--stat-method ' + STAT_METHOD + '</code>: a '
             'one-sided rank-sum test on per-exon motif hit counts, reduced to the smallest p over the 50-nt windows '
             'of each region. Stems and ranking are the raw p exactly as the tool reports it, colours are raw-p bins '
             '(&lt;0.001, &lt;0.01, &lt;0.05, \u22650.05) and dots have constant size because the released tool '
             'reports no enrichment ratio. No multiple-testing adjustment is applied, because none is part of the '
             'tool. The p-values come from the tie-corrected normal approximation and are anti-conservative in the '
             'sparse tail; the RBP order is the claim.</p>',
             '<p class="supp"><b>SUPPLEMENT — permutation-calibrated rank-sum.</b> The same statistic with its p '
             'recalibrated by Westfall\u2013Young label permutation (min-P over the windows of a region), and a BH q '
             'over the ' + str(BH_FAMILY) + ' pooled tests. Stems and ranking use the calibrated p, colour uses the '
             'calibrated q, and dot size shows the enrichment ratio the calibration summary provides. This says how '
             'much of the main figure\u2019s signal survives a permutation null of the same statistic.</p>',
             '<p>Included stems rise, skipped stems hang down, and the six panels of one figure share a y-scale. '
             'Best motifs are chosen after HGNC grouping; the twelve synthetic ESRP-like hexamers form one group. '
             'The <code>_noSpliceosome_noBroad</code> variants drop the lab core-spliceosome and broad-binder lists '
             'and always retain SRSF1. The authors\u2019 Fisher default and our audited Fisher-binary layer are not '
             'plotted here; they remain as comparison columns in each arm\u2019s rank workbook and in the archived '
             'version 3.1 figures.</p>',
             f'<p><a href="{archive_name}/index.html">Archived version 3.1 index (Fisher layers)</a> · '
             '<a href="naming_audit_v4.tsv">Naming audit</a> · <a href="selection_audit_v4.tsv">Selections</a> · '
             '<a href="exclusion_audit_v4.tsv">Exclusions</a> · '
             '<a href="positive_control_audit.tsv">QKI positive controls</a> · '
             '<a href="rank_agreement_summary_v4.tsv">Rank agreement</a> · '
             '<a href="y_scale_audit_v4.tsv">Y-scale caps</a> · '
             '<a href="truncation_audit_v4.tsv">Truncated stems</a> · '
             '<a href="figures_manifest_v4.tsv">Manifest (includes skipped layers)</a></p>']
    if skipped:
        parts.append('<h2>Layers not built</h2><table><tr><th>Arm</th><th>Layer</th><th>Missing input</th></tr>' +
                     ''.join(f'<tr><td>{html.escape(s["arm"])}</td><td>{html.escape(s["layer"])}</td>'
                             f'<td><code>{html.escape(s["file"])}</code></td></tr>' for s in skipped) + '</table>')

    def card(arm, kind, layer, variant):
        stem = f'{arm}_SE_{kind}_{layer}' + ('' if variant == 'main' else '_' + variant)
        return (f'<article><h4>{kind} · {variant}</h4><a href="{arm}/{stem}.png">'
                f'<img loading="lazy" src="{arm}/{stem}.png" alt="{stem}"></a>'
                f'<a href="{arm}/{stem}.svg">Editable SVG</a></article>')

    for rec in records:
        arm = rec['arm']
        c = rec['counts']
        parts.append(f'<section><h2>{html.escape(LABELS[arm])}</h2>'
                     f'<p>Included {c["n_up"]:,}; skipped {c["n_dn"]:,}; background {c["n_bg"]:,}. '
                     'Foreground FDR &lt;0.05 and |ΔPSI| ≥0.10; background FDR ≥0.5; per-sample junction reads ≥10; '
                     'gene baseMean &gt;50.</p>')
        for layer in LAYERS:
            if layer not in rec['layers']:
                continue
            s = rec['scales'][layer]
            axis = (f'Shared y-axis 0 to {s["ymax"]:g}. Cap = 1.25 × {s["cap_base"]:.1f}, the largest −log10 p '
                    f'outside the arm’s top RBP ({html.escape(s["top_rbp"])}); {s["n_entries_above_cap"]} '
                    f'of the layer’s tests exceed it (largest {s["layer_max"]:.1f}) and every plotted stem '
                    'above the cap is truncated with its value printed beside the dot.' if s['cap_applied'] else
                    f'Shared y-axis 0 to {s["ymax"]:g}. The data-driven cap (1.25 × {s["cap_base"]:.1f}) already '
                    f'exceeds the largest value ({s["layer_max"]:.1f}), so no stem is truncated.')
            parts.append(f'<h3 class="{"main" if layer == LAYERS[0] else "supp"}">{LAYER_TITLES[layer]}</h3>'
                         f'<p>{axis}</p><div class="grid">' +
                         ''.join(card(arm, kind, layer, variant) for kind in KINDS for variant in VARIANTS) +
                         '</div>')
        rows = [r for r in rec['rank_summary'] if r['mean_rho'] is not None]
        if rows:
            parts.append('<h3>Rank agreement between layers</h3><table>'
                         '<tr><th>Layer 1</th><th>Layer 2</th><th>Mean Spearman rho</th>'
                         '<th>Mean top-10 overlap</th></tr>' +
                         ''.join(f'<tr><td>{html.escape(r["layer1"])}</td><td>{html.escape(r["layer2"])}</td>'
                                 f'<td>{r["mean_rho"]:.3f}</td><td>{r["mean_top10_overlap"]:.1f}/10</td></tr>'
                                 for r in rows) + '</table>')
        parts.append(f'<p><a href="{arm}/{arm}_rank_comparison.xlsx">Rank workbook (.xlsx)</a> · '
                     f'<a href="{arm}/{arm}_rank_comparison_v4.tsv">Ranks TSV</a> · '
                     f'<a href="{arm}/{arm}_figure_provenance_v4.md">Provenance sidecar</a></p></section>')
    (out / 'index.html').write_text('\n'.join(parts + ['</html>']), encoding='utf-8')


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description='rMAPS3 SE rank-sum lollipops, version 4')
    ap.add_argument('--arms', nargs='+', required=True)
    ap.add_argument('--out-root', required=True)
    ap.add_argument('--released-root', required=True,
                    help='root holding <arm>/pVal.{up,dn}.vs.bg.RNAmap.txt from the released engine')
    ap.add_argument('--calibrated-root', required=True,
                    help='root holding <arm>/per_motif_regions.tsv from tools/calibrate_ranksum.py')
    ap.add_argument('--method-comparison', required=True,
                    help='method_rank_comparison.tsv from tools/compare_stat_methods.py')
    ap.add_argument('--archive-name', default='_v3.1_archive_2026-09-17')
    ap.add_argument('--event-sets-root', required=True)
    ap.add_argument('--alias-table', required=True,
                    help='HGNC alias table: table_name<TAB>hgnc_symbol')
    ap.add_argument('--naming-table', default=str(NAMING_TABLE))
    ap.add_argument('--esrp-table', default=str(ESRP_TABLE))
    ap.add_argument('--gtf', required=True,
                    help='GENCODE GTF used only for gene-name resolution')
    ap.add_argument('--spliceosome-list', required=True)
    ap.add_argument('--broad-binders-list', required=True)
    ap.add_argument('--provenance-source', action='append', default=[],
                    help='repeatable; extra file hashed into every figure manifest')
    ap.add_argument('--gate', default='persample10_bm50_bgfdr0.5')
    ap.add_argument('--top-n', type=int, default=10)
    args = ap.parse_args()
    out = Path(args.out_root)
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'Arial', 'svg.fonttype': 'none', 'svg.hashsalt': 'rmaps-v4',
                         'pdf.fonttype': 42})
    unknown = set(args.arms) - set(LABELS)
    if unknown:
        raise ValueError(f'Missing arm labels: {sorted(unknown)}')
    mappings, naming_rows, naming_stats = load_naming(args.naming_table, args.esrp_table, args.gtf, args.alias_table)
    write_tsv(out / 'naming_audit_v4.tsv',
              ['table_name', 'hgnc_symbol', 'source', 'evidence', 'ambiguity_note', 'n_motifs', 'merged_into'],
              naming_rows)
    (out / 'naming_report_v4.json').write_text(json.dumps(naming_stats, indent=2) + '\n', encoding='utf-8')
    lists = {'spliceosome_census': load_exclusion_list(args.spliceosome_list),
             'broad_binders': load_exclusion_list(args.broad_binders_list)}
    aliases = naming_stats['alias_mappings']
    lists = {k: {aliases.get(s, s) for s in v} for k, v in lists.items()}
    records, layouts, audits, exclusion_rows, controls, rank_summary = [], [], [], [], [], []
    skipped, manifest_rows, universe_rows, truncation_rows, scale_rows = [], [], [], [], []
    archive = out / args.archive_name
    for arm in args.arms:
        dest = out / arm
        dest.mkdir(exist_ok=True)
        counts, countpath = gate_counts(args.event_sets_root, arm, args.gate)
        base_sources = [countpath, Path(args.event_sets_root) / 'gate_sensitivity.tsv', Path(args.naming_table),
                        Path(args.alias_table), Path(args.esrp_table), Path(args.spliceosome_list),
                        Path(args.broad_binders_list), Path(__file__).resolve()]
        base_sources += [Path(p) for p in args.provenance_source]
        layers, sources, refinement = {}, list(base_sources), None
        released_path = Path(args.released_root) / arm / 'pVal.up.vs.bg.RNAmap.txt'
        if released_path.is_file() and (released_path.with_name('pVal.dn.vs.bg.RNAmap.txt')).is_file():
            entries, released_sources, n_motifs = released_entries(arm, args.released_root, mappings)
            layers['released_ranksum_rawP'] = entries
            sources += released_sources
        else:
            skipped.append({'arm': arm, 'layer': 'released_ranksum_rawP', 'file': str(released_path),
                            'reason': 'released rank-sum root table absent'})
        calib_path = Path(args.calibrated_root) / arm / 'per_motif_regions.tsv'
        calib_log = calib_path.with_name('command.log')
        if (calib_path.is_file() and calib_path.with_name('refinement_report.json').is_file()
                and calib_log.is_file() and 'exit=0' in calib_log.read_text()):
            entries, calib_sources, report = calibrated_entries(arm, args.calibrated_root, mappings)
            layers['calibrated_ranksum'] = entries
            sources += calib_sources
            refinement = dict(report)
            refinement['seed'] = report['seed']
        else:
            skipped.append({'arm': arm, 'layer': 'calibrated_ranksum', 'file': str(calib_path),
                            'reason': 'calibration summary absent or incomplete'})
        if not layers:
            print(f'SKIP {arm}: no layer input present', flush=True)
            continue
        crosschecked = (cross_check_layers(layers['released_ranksum_rawP'], layers['calibrated_ranksum'], arm)
                        if len(layers) == 2 else None)
        arm_exclusion, dropped = exclusion_audit(arm, layers[next(iter(layers))], lists)
        exclusion_rows.extend(arm_exclusion)
        scales = {layer: y_scale(entries, layer) for layer, entries in layers.items()}
        selections = {}
        for layer, entries in layers.items():
            for variant in VARIANTS:
                chosen = entries if variant == 'main' else [e for e in entries if e['_exclusion_symbol'] not in dropped]
                for kind in KINDS:
                    panels = select_panels(chosen, layer, kind, args.top_n)
                    selections[(variant, kind, layer)] = panels
                    if arm.startswith('QKI_KO'):
                        for panel in [('INCLUDED', 'Upstream Intron'), ('SKIPPED', 'Downstream Intron')]:
                            first = panels[panel][0]
                            controls.append({'arm': arm, 'layer': layer, 'variant': variant, 'kind': kind,
                                             'panel': ' × '.join(panel), 'first': first['_label'],
                                             'first_p': first['_p'], 'first_q': first['_q'],
                                             'pass': first['_label'] == 'QKI'})
                    for (d, region), rows in panels.items():
                        for i, r in enumerate(rows, 1):
                            audits.append({'arm': arm, 'layer': layer, 'variant': variant, 'type': kind,
                                           'direction': d, 'pooled_region': region, 'rank': i, 'label': r['_label'],
                                           'figure_tick_label': r['_tick'],
                                           'source_rbp': r['RBP'], 'motif_key': r['motif_key'],
                                           'selected_sub_region': r['region'], 'p': r['_p'], 'q': r['_q'],
                                           'native_p': r['_native_p'],
                                           'enrichment_ratio': None if layer == LAYERS[0] else r['_ratio'],
                                           'k_sig_motifs': r['_k'], 'n_motifs_in_group': r['_n'],
                                           'hgnc_symbol': r['_hgnc_symbol'], 'naming_action': r['_naming_action'],
                                           'calib_perms_used': r['_calib_perms_used'], 'calib_stage': r['_calib_stage']})
        for (variant, kind, layer), panels in sorted(selections.items()):
            stem = dest / (f'{arm}_SE_{kind}_{layer}' + ('' if variant == 'main' else '_' + variant))
            report = draw_figure(arm, layer, kind, panels, counts, refinement, scales[layer], stem, args,
                                 excluded=variant != 'main')
            layouts.append(report)
            for row in report['truncated_stems']:
                truncation_rows.append({'arm': arm, 'layer': layer, 'variant': variant, 'kind': kind,
                                        'shared_ymax': scales[layer]['ymax'], **row})
        v31_tsv = archive / arm / f'{arm}_rank_comparison.tsv'
        summaries, notes, has_methods = build_rank_comparison(
            arm, {l: select_panels(layers[l], l, 'byRBP', 10000) for l in layers},
            v31_tsv if v31_tsv.is_file() else None, args.method_comparison, refinement, dest)
        rank_summary.extend(summaries)
        universe_rows.extend(notes)
        if v31_tsv.is_file():
            sources.append(v31_tsv)
        if has_methods:
            sources.append(Path(args.method_comparison))
        prov = [f'# {arm} figure provenance — version 4 (rank-sum layers)',
                '',
                '- Scope: Homo sapiens / GRCh38 (hg38) / GENCODE v49; skipped-exon (SE) events only. '
                'GTF gene_name values verify exact names and HGNC alias targets. No mouse arm exists.',
                f'- Layers built: {", ".join(sorted(layers))}.',
                f'- Layers skipped for this arm: '
                f'{", ".join(s["layer"] for s in skipped if s["arm"] == arm) or "none"}.',
                f'- Gate {args.gate}, rule {arm.rsplit("_", 1)[1]}: included {counts["n_up"]}, '
                f'skipped {counts["n_dn"]}, background {counts["n_bg"]}; genes absent from the DESeq2 table retained '
                f'({counts["n_expr_unknown_in_fg"]} foreground / {counts["n_expr_unknown_in_bg"]} background events).',
                f"- MAIN layer released_ranksum_rawP: authors' released rMAPS3 at commit {RELEASED_COMMIT}, "
                f'--stat-method {STAT_METHOD}; raw regional-minimum p read verbatim from '
                'pVal.{up,dn}.vs.bg.RNAmap.txt. No adjustment, no ratio, constant dot size. The event-set md5s in the '
                'run log are identical to the audited and released-Fisher runs of the same arm.',
                '- Pooled regions: Upstream Intron = min(upstreamExonIntron, upstreamIntron); Exon Body = '
                'min(targetExon-5prime, targetExon-3prime); Downstream Intron = min(downstreamIntron, '
                'downstreamExonIntron). The two flanking-exon sub-regions are not plotted.']
        if refinement:
            prov += [f'- SUPPLEMENT layer calibrated_ranksum: same statistic; Westfall-Young min-P label permutation, '
                     f'seed {refinement["seed"]}, stage 1 B = {refinement["stage1_permutations"]}, stage 2 B = '
                     f'{refinement["stage2_permutations"]} for stage-1 p <= {refinement["refine_threshold"]} '
                     f'({refinement["pairs_refined"]} of {refinement["pairs_total"]} motif-direction pairs refined, '
                     f'{refinement["pairs_deferred_by_cap"]} deferred by cap); BH q over '
                     f'{refinement["bh_family_size"]} tests; permutation floor '
                     f'1/{refinement["stage2_permutations"] + 1}.',
                     '- Stem and rank on the supplement use the calibrated p; colour uses the calibrated BH q; dot '
                     'size uses the enrichment ratio supplied by the calibration summary. This follows the v3.1 '
                     'calibrated-layer convention.']
        if crosschecked:
            prov.append(f'- Cross-check: the native rank-sum p of all {crosschecked} pooled tests in the calibration '
                        'summary equals the released root table value to a relative tolerance of 1e-9, so both layers '
                        'describe the same statistic.')
        prov += ['- No p-value or q-value is recomputed by this script. Colour bins: raw p <0.001, <0.01, <0.05, '
                 '>=0.05; calibrated q <0.01, <0.05, <0.10, >=0.10.',
                 '- Representative sub-region for a pooled region: smallest native rank-sum p, then enrichment ratio '
                 'descending, then sub-region name.',
                 '- byRBP: best motif per HGNC group by layer p, ties by enrichment ratio descending then display '
                 'label then motif key. No k/n numerals are drawn on any figure; the per-group significant-motif '
                 'counts are in selection_audit_v4.tsv.',
                 '- Tick labels: HGNC group everywhere, except that byMotif panels label the twelve synthetic '
                 'ESRP-like hexamers by sequence, because they carry no gene name and would otherwise repeat the '
                 'same word. Ranking is unaffected: it uses the group label, not the tick label. Every drawn tick '
                 'label is recorded in selection_audit_v4.tsv.',
                 '- Naming: display and grouping use the supplied HGNC alias table, otherwise the exact GENCODE v49 '
                 'gene_name. Mappings: ' + '; '.join(f'{n} → {s}' for n, s in naming_stats['alias_mappings'].items()) +
                 '. HGNC-ambiguous aliases: 9G8 (SRSF7/SLU7; supplied resolution SRSF7) and SRp40 (SRSF5/NOLC1; '
                 'supplied resolution SRSF5). The twelve synthetic ESRP-like hexamers keep their group and are '
                 'omitted from the gene naming audit.',
                 f'- Exclusions in the _noSpliceosome_noBroad variants: {len(dropped)} distinct confirmed symbols '
                 f'removed: {", ".join(sorted(dropped)) or "none"}. SRSF1 is always retained.',
                 '- One y limit per arm x layer, shared by all six panels and by byRBP/byMotif and both exclusion '
                 'variants. The limit is a data-driven cap: 1.25 x the largest -log10 p among RBPs other than the '
                 'single top RBP of the arm, rounded up to a clean tick (calibrated: rounded up to an integer, '
                 'minimum 4; raw p: a 1/2/2.5/5/10 scale). When that cap already exceeds the largest value nothing '
                 'is truncated and the limit is the plain rounded maximum. Stems above the cap are drawn to the cap, '
                 'carry a two-diagonal axis-break glyph and have their exact value printed beside the dot; the '
                 'legend box says so. Every truncated stem is listed in truncation_audit_v4.tsv. Included up, '
                 'skipped down; gene model between the rows.']
        for layer in sorted(scales):
            s = scales[layer]
            prov.append(f'  - {layer}: top RBP {s["top_rbp"]}, largest -log10 p {s["layer_max"]:.4f}, largest '
                        f'outside that RBP {s["cap_base"]:.4f}, cap {s["cap"]:.4f}, shared y limit {s["ymax"]:g}, '
                        f'tick step {s["step"]:g}, cap applied {s["cap_applied"]}, entries above the cap '
                        f'{s["n_entries_above_cap"]}.')
        prov += [
                 f'- Rendering: Arial, editable-text SVG (svg.fonttype none, hashsalt rmaps-v4) plus a {PNG_DPI}-dpi '
                 'PNG of the same 32 × 28 in canvas. Every figure passed the automated text-overlap and '
                 'out-of-bounds audit.',
                 '- Source SHA256:']
        prov += [f'  - {str(p)}: {digest(p, "sha256")}' for p in sources]
        (dest / f'{arm}_figure_provenance_v4.md').write_text('\n'.join(prov) + '\n', encoding='utf-8')
        (dest / 'command_v4.log').write_text(
            ' '.join([sys.executable, str(Path(__file__).resolve())] + sys.argv[1:]) + '\n', encoding='utf-8')
        (dest / 'versions_v4.txt').write_text(
            f'Python\t{sys.version.split()[0]}\nmatplotlib\t{matplotlib.__version__}\n'
            f'openpyxl\t{openpyxl.__version__}\nnumpy\t{np.__version__}\nPNG_DPI\t{PNG_DPI}\n', encoding='utf-8')
        for layer, s in sorted(scales.items()):
            scale_rows.append({'arm': arm, 'layer': layer, **s})
        records.append({'arm': arm, 'counts': counts, 'layers': sorted(layers), 'rank_summary': summaries,
                        'scales': scales})
        print(f'OK {arm}: layers={sorted(layers)} figures={len(selections)}', flush=True)
    write_index(out, records, skipped, args.archive_name)
    if controls:
        write_tsv(out / 'positive_control_audit.tsv', list(controls[0]), controls)
    write_tsv(out / 'selection_audit_v4.tsv', list(audits[0]), audits)
    write_tsv(out / 'exclusion_audit_v4.tsv', list(exclusion_rows[0]), exclusion_rows)
    write_tsv(out / 'rank_agreement_summary_v4.tsv', list(rank_summary[0]), rank_summary)
    write_tsv(out / 'y_scale_audit_v4.tsv', list(scale_rows[0]), scale_rows)
    write_tsv(out / 'truncation_audit_v4.tsv',
              ['arm', 'layer', 'variant', 'kind', 'shared_ymax', 'direction', 'pooled_region', 'rank', 'label',
               'tick', 'motif_key', 'log10p', 'printed'],
              truncation_rows or [{'arm': 'NONE', 'layer': 'NONE', 'variant': 'NONE', 'kind': 'NONE',
                                   'shared_ymax': 'NA', 'direction': 'NA', 'pooled_region': 'NA', 'rank': 'NA',
                                   'label': 'NA', 'tick': 'NA', 'motif_key': 'NA', 'log10p': 'NA',
                                   'printed': 'no stem exceeded its cap in this build'}])
    if universe_rows:
        write_tsv(out / 'rank_universe_notes_v4.tsv', list(universe_rows[0]), universe_rows)
    (out / 'layout_report_v4.json').write_text(json.dumps(layouts, indent=2), encoding='utf-8')
    for s in skipped:
        manifest_rows.append({'arm': s['arm'], 'layer': s['layer'], 'status': 'skipped_input_missing',
                              'file': s['file'], 'bytes': 'NA', 'md5': 'NA'})
    built = {p for rec in records for p in (out / rec['arm']).rglob('*')}
    built |= {p for p in out.glob('*') if p.is_file()}
    for p in sorted(built):
        if not p.is_file() or p.name == 'figures_manifest_v4.tsv':
            continue
        layer = next((l for l in LAYERS if l in p.name), 'ALL')
        manifest_rows.append({'arm': p.parent.name if p.parent != out else 'ALL', 'layer': layer,
                              'status': 'built', 'file': str(p.resolve()), 'bytes': p.stat().st_size,
                              'md5': digest(p)})
    write_tsv(out / 'figures_manifest_v4.tsv', ['arm', 'layer', 'status', 'file', 'bytes', 'md5'], manifest_rows)
    print(f'OK v4 index, audits, layout report and manifest: {out}', flush=True)
    if skipped:
        print('SKIPPED LAYERS: ' + '; '.join(f'{s["arm"]}/{s["layer"]}' for s in skipped), flush=True)


if __name__ == '__main__':
    main()
