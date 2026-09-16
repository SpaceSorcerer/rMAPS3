from pyx import *
import re, os, sys, logging, time, argparse, numpy, subprocess, shutil
import glob
import operator
from scipy import stats
import warnings
import timeit
import multiprocessing
import pickle
import types
from pathlib import Path

from rmaps_core import drawutils
from rmaps_core.genome_access import load_genome, fetch_seq as fetch_seq_from_fasta
from rmaps_core.stat_utils import normalize_stat_method, pvalue_header_label
from rmaps_core.se_windows import (REGION_NAMES, read_event_sets, event_regions,
                                   overlapping_hits, binary_windows, binary_table,
                                   binary_density, WindowCounts, max_motif_width)
from rmaps_core.output_utils import ensure_output_directory, write_run_manifest

DRAW_MOTIF_PLOTS = (
    os.environ.get("RMAPS_FORCE_MOTIF_FALLBACK") != "1"
    and drawutils.use_unicode_text_engine()
)

def run_command(cmd):
    completed = subprocess.run(cmd, capture_output=True, text=True)
    status = completed.returncode
    output = (completed.stdout or "") + (completed.stderr or "")
    return status, output


def copy_file(src, dst):
    try:
        shutil.copy2(src, dst)
        return 0, ""
    except Exception as exc:
        return 1, str(exc)


def _safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _mean_psi_field(field):
    values = []
    for token in field.replace('"', '').split(','):
        val = _safe_float(token)
        if val is not None:
            values.append(val)
    if not values:
        return None
    return sum(values) / float(len(values))


def setup_runtime():
    parser = argparse.ArgumentParser(
        description=
        'Making motif map from a list of motifs and three lists of exon coordinates'
    )
    parser.add_argument('-m',
                        '--motif',
                        dest='motif',
                        required=True,
                        help='List of motifs. A file with only one header')
    parser.add_argument(
        '-k',
        '--knownMotifs',
        dest='knownMotifs',
        required=True,
        help=
        'A tab delimited list of motifs with the following headers. MotifName regularExpression'
    )

    parser.add_argument('-f',
                        '--fasta-root',
                        '--fastaRoot',
                        dest='fastaRoot',
                        required=True,
                        help='FASTA root path. e.g., /path/to/genomedata')
    parser.add_argument(
        '-g',
        '--genome',
        dest='genome',
        required=True,
        help='The UCSC genome build name, hg38, hg19, mm10, or dm3')
    parser.add_argument('-o',
                        '--output',
                        dest='output',
                        required=True,
                        help='output directory')
    parser.add_argument('-r',
                        '--rMATS',
                        dest='rMATS',
                        required=True,
                        help='an rMATS output file from SE event')
    parser.add_argument('-mi',
                        '--miso',
                        dest='miso',
                        required=True,
                        help='an miso output file from SE event')
    start = timeit.default_timer()

    # AUDIT F13: advertise exactly the schema enforced by the coordinate validator.
    parser.add_argument(
        '-u',
        '--up',
        dest='up',
        required=True,
        help=
        'Eight tab-separated columns with header: chr strand exonStart exonEnd firstExonStart firstExonEnd secondExonStart secondExonEnd'
    )
    parser.add_argument(
        '-d',
        '--down',
        dest='dn',
        required=True,
        help=
        'Eight tab-separated columns with header: chr strand exonStart exonEnd firstExonStart firstExonEnd secondExonStart secondExonEnd'
    )
    parser.add_argument(
        '-b',
        '--background',
        dest='bg',
        required=True,
        help=
        'Eight tab-separated columns with header: chr strand exonStart exonEnd firstExonStart firstExonEnd secondExonStart secondExonEnd'
    )
    parser.add_argument('--label',
                        type=str,
                        dest='label',
                        default="RBP",
                        help='Label of motif. e.g., TG-rich')
    parser.add_argument(
        '--intron',
        type=int,
        dest='intron',
        default=250,
        help='grab sequence up to a this number of NTs away from the exon junction'
    )
    parser.add_argument(
        '--exon',
        type=int,
        dest='exon',
        default=50,
        help=
        'grab sequence up to this number of NTs from both ends of the exon body')
    parser.add_argument(
        '--window',
        type=int,
        dest='window',
        default=50,
        help='number of NTs examined for the given motif at a time')
    parser.add_argument('--step',
                        type=int,
                        dest='step',
                        default=1,
                        help='slide window by this number of NTs at a time')
    # AUDIT R1: exon windows inherit --window unless explicitly overridden.
    parser.add_argument('--exon-window', type=int, default=None)
    parser.add_argument('--sigFDR',
                        type=float,
                        dest='sigFDR',
                        default=0.05,
                        help='FDR cutoff for significant events.')
    parser.add_argument(
        '--sigDeltaPSI',
        type=float,
        dest='sigDeltaPSI',
        default=0.05,
        help='inclusion level difference cutoff for significant events.')
    parser.add_argument('--separate',
                        dest='separate',
                        default=False,
                        action='store_const',
                        const=True)
    # AUDIT F13: cross-set overlap requires explicit acknowledgement.
    parser.add_argument('--allow-overlap', action='store_true')
    # AUDIT F1: retain directional Fisher by default; expose two-sided testing.
    parser.add_argument('--fisher-alternative', choices=['greater', 'two-sided'], default='greater')
    # AUDIT F20: bound automatic concurrency and permit explicit allocation.
    parser.add_argument('--workers', type=int, default=max(1, min(4, multiprocessing.cpu_count() - 1)))
    # AUDIT F7: positional scientific outputs are preserved by default.
    parser.add_argument('--delete-temp', action='store_true')
    # AUDIT F8: existing outputs require explicit reuse authorization.
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if args.exon_window is None:
        args.exon_window = args.window
    if any(value <= 0 for value in (args.intron, args.exon, args.window, args.exon_window, args.step, args.workers)):
        parser.error('intron, exon, window, step and workers must be positive')
    stat_method = normalize_stat_method(os.environ.get('RMAPS_STAT_METHOD', 'fisher'))
    try:
        stat_permutations = max(50, int(os.environ.get('RMAPS_STAT_PERMUTATIONS', '500')))
    except Exception:
        stat_permutations = 500
    try:
        stat_seed = max(0, int(os.environ.get('RMAPS_STAT_SEED', '1337')))
    except Exception:
        stat_seed = 1337


    def listToString(x):  ## log command
        rVal = ''
        for a in x:
            rVal += a + ' '
        return rVal


    outDir = args.output
    ensure_output_directory(outDir, overwrite=args.overwrite)
    run_outputs = []
    run_motif_ids = []
    outPath = os.path.abspath(outDir)
    exonDir = args.output + '/exon'
    os.makedirs(exonDir, exist_ok=True)
    exonPath = os.path.abspath(exonDir)
    fastaDir = args.output + '/fasta'
    os.makedirs(fastaDir, exist_ok=True)
    fastaPath = os.path.abspath(fastaDir)
    mapsDir = args.output + '/maps'
    os.makedirs(mapsDir, exist_ok=True)
    mapsPath = os.path.abspath(mapsDir)
    tempDir = args.output + '/temp'
    os.makedirs(tempDir, exist_ok=True)
    tempPath = os.path.abspath(tempDir)
    positionalPath = os.path.join(os.path.abspath(outDir), 'positional')
    os.makedirs(positionalPath, exist_ok=True)
    scriptPath = os.path.abspath(os.path.dirname(__file__))
    binPath = os.path.abspath(os.path.join(scriptPath, '..', 'bin'))

    VER = "2.0.1"
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(message)s',
        filename=outPath + '/log.motifMap' + '.txt',
        filemode='w')

    logging.debug('Start the program with [%s]\n', listToString(sys.argv))
    run_outputs.append(os.path.join(outPath, 'log.motifMap.txt'))
    startTime = time.time()
    logging.debug('motifTools version: %s', VER)

    fasta_root_PATH = args.fastaRoot

    rMATS = args.rMATS
    miso = args.miso
    up = args.up
    dn = args.dn
    bg = args.bg
    cdist_up = {}
    cdist_dn = {}
    cdist_bg = {}
    global mFile
    if rMATS == "NA" and miso == "NA" and (up == "NA" or dn == "NA"
                                           or bg == "NA"):  ## not going to work
        print("Incorrect Input! Need to have rMATS input or miso input or coordinates for all three types of exons")
        sys.exit(-99)

    iLen = args.intron  # intron length to examine
    eLen = args.exon
    wLen = args.window
    sLen = args.step

    sigFDR = args.sigFDR
    sigDeltaPSI = args.sigDeltaPSI
    bgFDR = 0.5
    minPSI = 0.85
    maxPSI = 0.15
    u = {}
    d = {}
    b = {}

    if args.motif == 'NA':  ## user didn't give motif
        pass
    else:  ## optional motif was given
        mFile = open(args.motif)
    kFile = open(args.knownMotifs)
    region = [
        "UpstreamExon", "UpstreamExonIntron", "UpstreamIntron", "TargetExon",
        "DownstreamIntron", "DownstreamExonIntron", "DownstreamExon"
    ]
    rName = [
        "UpstreamExon", "UpstreamExonIntron", "UpstreamIntron", "TargetExon",
        "TargetExon", "DownstreamIntron", "DownstreamExonIntron", "DownstreamExon"
    ]
    pRegion = ["UpstreamExonIntron", "DownstreamIntron", "DownstreamExon"]
    nRegion = ["UpstreamExon", "UpstreamIntron", "DownstreamExonIntron"]

    totalExonCount = {}
    motifLabel = args.label
    boxHeight = 1.0
    uNum = 0
    dNum = 0
    bNum = 0
    if miso != "NA":  ## got miso input here
        rMATS = tempPath + '/converted.rMATS.se.txt'
        convcmd = ['perl', binPath + '/miso2rMATS.SE.pl', '1', '100', miso, rMATS]
        status, output = run_command(convcmd)
        logging.debug("Conversion from miso to rMATS is done. Status: %s" % status)
        if (int(status) != 0):  ## it did not go well
            logging.debug("error in converting miso file")
            logging.debug("error detail: %s" % output)
            sys.exit(-101)
        logging.debug(output)
        run_outputs.append(rMATS)



    globals().update(locals())

def fetch_seq(genome, strand, chr, start, end):
    """
    Fetch sequence from genome using the shared pyfaidx-based implementation.
    """
    return fetch_seq_from_fasta(genome, strand, chr, start, end)


def makeInputFiles(
    rmats
):  ## make input coordinate files. Don't need to this if rMATS was "NA"
    global up, dn, bg
    global totalExonCount
    global nu, nd, nb
    logging.debug("Making input files from rMATS output")
    if rmats == "NA":  ## do not need to make new up,dn,bg. Just copy it over
        status, output = copy_file(up, exonPath + '/up.coord.txt')
        logging.debug("Copying up file is done. Status: %s" % status)
        if (int(status) != 0):  ## it did not go well
            logging.debug("error in copying up file")
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)

        status, output = copy_file(dn, exonPath + '/dn.coord.txt')
        logging.debug("Copying dn file is done. Status: %s" % status)
        if (int(status) != 0):  ## it did not go well
            logging.debug("error in copying dn file")
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)

        status, output = copy_file(bg, exonPath + '/bg.coord.txt')
        logging.debug("Copying bg file is done. Status: %s" % status)
        if (int(status) != 0):  ## it did not go well
            logging.debug("error in copying bg file")
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)

        nb = wccount(exonPath + '/bg.coord.txt') - 1
        nd = wccount(exonPath + '/dn.coord.txt') - 1
        nu = wccount(exonPath + '/up.coord.txt') - 1
        totalExonCount = {'up': nu, 'dn': nd, 'bg': nb}

        return

    rFile = open(rmats, 'r')

    uFile = open(exonPath + '/up.coord.txt', 'w')
    dFile = open(exonPath + '/dn.coord.txt', 'w')
    bFile = open(exonPath + '/bg.coord.txt', 'w')

    header = 'chr\tstrand\texonStart\texonEnd\tfirstExonStart\tfirstExonEnd\tsecondExonStart\tsecondExonEnd'
    uFile.write(header + '\n')
    dFile.write(header + '\n')
    bFile.write(header + '\n')

    line = rFile.readline()
    fallback_bg = []
    for line in rFile:  ## process each line

        ele = line.strip().split('\t')
        if len(ele) < 11:
            continue
        key = ':'.join(ele[3:7])
        value = [1, '\t'.join(ele[3:11])]
        fdr = _safe_float(ele[-4])
        deltaPSI = _safe_float(ele[-1])
        if fdr is None:
            continue
        if fdr < sigFDR and deltaPSI is not None:  ## it could be significant
            if deltaPSI >= sigDeltaPSI:  ## it's upregulated. high in sample 1
                u[key] = value
            elif deltaPSI <= -sigDeltaPSI:  ## it's downregulated. high in sample 2
                d[key] = value
        elif fdr > bgFDR:  ## it could be background
            PSI_1 = _mean_psi_field(ele[-3])
            PSI_2 = _mean_psi_field(ele[-2])

            if PSI_1 is not None and PSI_2 is not None and min(PSI_1, PSI_2) < minPSI and max(
                    PSI_1, PSI_2) > maxPSI:  ## it is a background event
                b[key] = value

        # Keep non-significant candidates as fallback background when strict background is empty.
        if deltaPSI is not None and fdr >= sigFDR:
            fallback_bg.append((abs(deltaPSI), -fdr, key, value))

    logging.debug(
        "Done populating initial dictionaries with possible duplicates")
    logging.debug("Number of up, down, and background exons are: %d, %d, %d" %
                  (len(u), len(d), len(b)))

    logging.debug("Removing exons included in more than one dictionaries..")

    for key in u:  ## going through u
        if key in d:  ## same exon in d
            d[key][0] += 1
            u[key][0] += 1
        if key in b:  ## same exon in b
            b[key][0] += 1
            u[key][0] += 1

    for key in d:  ## going through d
        if key in b:  ## same exon in b
            b[key][0] += 1
            d[key][0] += 1

    nu = 0
    nd = 0
    nb = 0
    for key in u:
        if u[key][0] == 1:  ## it is unique
            nu += 1
            uFile.write(u[key][1] + '\n')
    for key in d:
        if d[key][0] == 1:  ## it is unique
            nd += 1
            dFile.write(d[key][1] + '\n')

    if len(b) == 0 and fallback_bg:
        logging.debug("No strict background events found; building fallback background from non-significant events")
        fallback_bg.sort()
        target_bg = max(200, int((max(nu, 1) + max(nd, 1)) / 4))
        added = 0
        for _, _, key, value in fallback_bg:
            if key in u or key in d or key in b:
                continue
            b[key] = [1, value[1]]
            added += 1
            if added >= target_bg:
                break

    for key in b:
        if b[key][0] == 1:  ## it is unique
            nb += 1
            bFile.write(b[key][1] + '\n')

    logging.debug("Number of upregulated (high in sample_1) exons: %d" % nu)
    logging.debug("Number of downregulated (high in sample_2) exons: %d" % nd)
    logging.debug("Number of background exons: %d" % nb)

    rFile.close()
    uFile.close()
    dFile.close()
    bFile.close()

    totalExonCount = {'up': nu, 'dn': nd, 'bg': nb}

    logging.debug("Done making input file from rMATS")


def getFasta(t):
    # AUDIT F6: use true features, clipping chromosome edges before orientation.
    region_data[t] = [event_regions(genome, event, iLen, eLen, motif_padding) for event in event_sets[t]]
    feature_indices = (0, 1, 2, 3, 5, 6, 7)
    for rg, index in zip(region, feature_indices):
        filename = os.path.join(fastaPath, t + '.' + rg + '.fasta')
        with open(filename, 'w') as handle:
            for event, regions in zip(event_sets[t], region_data[t]):
                handle.write('>' + event.event_id + '\n' + regions[index].sequence + '\n')
        run_outputs.append(filename)


def findAll(re_motif, s_seq):
    # AUDIT F9: preserve overlapping matches and each individual span.
    hits = overlapping_hits(re_motif, s_seq)
    return [list(hit) for hit in hits], [start for start, _ in hits]

def initMotif(mF, mc):  ## initialize motif counts

    motifs = []
    logging.debug("Initializing motif counts")
    header = mF.readline()
    for line in mF:  ## process each motif
        motifs.append(line.strip().split('\t')[1])
    logging.debug("The number of motifs: %d" % len(motifs))
    return motifs


def prepare_sequence_store():
    # AUDIT R4: workers reopen immutable mmap arrays, never pickle region_data.
    global sequence_paths, sequence_arrays, sequence_origins, sequence_lengths, sequence_labels
    rows = [regions for label in ('up', 'dn', 'bg') for regions in region_data[label]]
    sequence_paths = []
    sequence_origins = numpy.array([[r.origin for r in row] for row in rows], dtype=numpy.int32)
    sequence_lengths = numpy.array([[len(r.sequence) for r in row] for row in rows], dtype=numpy.int32)
    sequence_labels = numpy.concatenate([numpy.full(totalExonCount[label], code, dtype=numpy.uint8)
                                         for code, label in enumerate(('up', 'dn', 'bg'))])
    for index in range(8):
        filename = os.path.join(tempPath, f'sequence_region_{index}.npy')
        numpy.save(filename, numpy.array([row[index].sequence for row in rows]))
        sequence_paths.append(filename)
    metadata = os.path.join(tempPath, 'sequence_metadata.npz')
    numpy.savez(metadata, origins=sequence_origins, lengths=sequence_lengths, labels=sequence_labels)
    run_outputs.extend(sequence_paths + [metadata])
    sequence_arrays = [numpy.load(filename, mmap_mode='r', allow_pickle=False) for filename in sequence_paths]
    region_data.clear()


def countMotif(mc, ttName, ttMotif, motifs):
    # AUDIT R4: sparse intervals plus range differences make work proportional to hits.
    cdist = {label: {} for label in ('up', 'dn', 'bg')}
    padding = max(max_motif_width(pattern) for pattern in motifs) - 1
    for index, name in enumerate(REGION_NAMES):
        is_exon = index in (0, 3, 4, 7)
        length = eLen if is_exon else iLen
        window = args.exon_window if is_exon else wLen
        n_windows = max(0, (length - window) // sLen + 1)
        origins = sequence_origins[:, index]
        sizes = sequence_lengths[:, index]
        first = numpy.maximum(0, (origins + sLen - 1) // sLen)
        last = numpy.minimum(n_windows, (origins + sizes - window) // sLen + 1)
        valid = first < last
        elig_lo = numpy.where(valid, first * sLen, -1).astype(numpy.int32)
        elig_hi = numpy.where(valid, (last - 1) * sLen + 1, -1).astype(numpy.int32)
        eligible_diff = numpy.zeros((3, n_windows + 1), dtype=numpy.int32)
        numpy.add.at(eligible_diff, (sequence_labels[valid], first[valid]), 1)
        numpy.add.at(eligible_diff, (sequence_labels[valid], last[valid]), -1)
        hit_diff = numpy.zeros_like(eligible_diff)
        hit_events, hit_starts, hit_ends = [], [], []
        for row, value in enumerate(sequence_arrays[index]):
            origin = int(origins[row])
            crop_lo = max(0, -padding - origin)
            crop_hi = min(len(value), length + padding - origin)
            sequence = str(value[crop_lo:crop_hi])
            offset = origin + crop_lo
            spans = sorted((a + offset, b + offset) for pattern in motifs
                           for a, b in overlapping_hits(pattern, sequence)
                           if a + offset < length and b + offset > 0)
            merged_lo = merged_hi = -1
            for start, end in spans:
                hit_events.append(row)
                hit_starts.append(start)
                hit_ends.append(end)
                lo = max(int(first[row]), (start - window) // sLen + 1)
                hi = min(int(last[row]), (end - 1) // sLen + 1)
                if lo >= hi:
                    continue
                if lo <= merged_hi:
                    merged_hi = max(merged_hi, hi)
                else:
                    if merged_lo >= 0:
                        hit_diff[sequence_labels[row], merged_lo] += 1
                        hit_diff[sequence_labels[row], merged_hi] -= 1
                    merged_lo, merged_hi = lo, hi
            if merged_lo >= 0:
                hit_diff[sequence_labels[row], merged_lo] += 1
                hit_diff[sequence_labels[row], merged_hi] -= 1
        positives = numpy.cumsum(hit_diff[:, :-1], axis=1)
        eligible = numpy.cumsum(eligible_diff[:, :-1], axis=1)
        for code, label in enumerate(('up', 'dn', 'bg')):
            cdist[label][index] = {j: WindowCounts(int(positives[code, j]), int(eligible[code, j]), totalExonCount[label])
                                   for j in range(n_windows)}
        filename = os.path.join(positionalPath, ttName + '.' + ttMotif + '.' + name + '.hits.npz')
        numpy.savez_compressed(filename, schema_version=numpy.int16(2),
                               event_index=numpy.asarray(hit_events, dtype=numpy.int32),
                               hit_start=numpy.asarray(hit_starts, dtype=numpy.int16),
                               hit_end=numpy.asarray(hit_ends, dtype=numpy.int16),
                               elig_lo=elig_lo, elig_hi=elig_hi, set_label=sequence_labels,
                               window=numpy.int32(window), step=numpy.int32(sLen), region_length=numpy.int32(length))
        run_outputs.append(filename)
    return cdist

def drawNode(c, eH, eW, iW, indent, gap, sGap, scale):  ## draw node

    mls = style.linewidth.THIck

    upstreamExon = path.path(path.moveto(0, 0), path.lineto(eW * scale, 0),
                             path.lineto(scale * eW, scale * eH),
                             path.lineto(0, scale * eH),
                             path.lineto(scale * indent, scale * eH / 2),
                             path.lineto(0, 0))
    intron_1 = path.line(scale * eW, scale * eH / 2, scale * (eW + iW),
                         scale * eH / 2)
    divider_1 = [
        path.line((eW + iW - sGap) * scale, (eH / 2 - gap) * scale,
                  (eW + iW + sGap) * scale, (eH / 2 + gap) * scale),
        path.line((eW + iW - sGap + gap) * scale, (eH / 2 - gap) * scale,
                  (eW + iW + sGap + gap) * scale, (eH / 2 + gap) * scale)
    ]
    intron_2 = path.line((eW + iW + gap) * scale, scale * eH / 2,
                         (eW + iW + gap + iW) * scale, scale * eH / 2)

    targetExon = path.rect((eW + iW + gap + iW) * scale, 0,
                           (eW + gap + gap + eW) * scale, eH * scale)

    newX = eW + iW + gap + iW + eW + gap + gap + eW
    intron_3 = path.line(newX * scale, (eH / 2) * scale, (newX + iW) * scale,
                         scale * eH / 2)
    divider_2 = [
        path.line((newX + iW - sGap) * scale, (eH / 2 - gap) * scale,
                  (newX + iW + sGap) * scale, (eH / 2 + gap) * scale),
        path.line((newX + iW - sGap + gap) * scale, (eH / 2 - gap) * scale,
                  (newX + iW + sGap + gap) * scale, (eH / 2 + gap) * scale)
    ]
    intron_4 = path.line((newX + iW + gap) * scale, scale * eH / 2,
                         (newX + iW + gap + iW) * scale, (eH / 2) * scale)
    newX = newX + iW + gap + iW
    downstreamExon = path.path(
        path.moveto(newX * scale, 0), path.lineto((newX + eW) * scale, 0),
        path.lineto((newX + eW - indent) * scale, scale * eH / 2),
        path.lineto((newX + eW) * scale, eH * scale),
        path.lineto(newX * scale, eH * scale), path.lineto(newX * scale, 0))

    c.stroke(upstreamExon, [mls, deco.filled([color.gray(0.8)])])
    c.stroke(intron_1, [mls])
    c.stroke(divider_1[0], [mls])
    c.stroke(divider_1[1], [mls])
    c.stroke(intron_2, [mls])
    c.stroke(targetExon, [mls, deco.filled([color.rgb.green])])
    c.stroke(intron_3, [mls])
    c.stroke(divider_2[0], [mls])
    c.stroke(divider_2[1], [mls])
    c.stroke(intron_4, [mls])
    c.stroke(downstreamExon, [mls, deco.filled([color.gray(0.8)])])

    logging.debug("Done drawing node structure")


def drawBox(c,
            max_point_value,
            exon_height,
            exon_width,
            intron_width,
            indent,
            divider_gap,
            slope_gap,
            scale,
            exon_length,
            intron_length,
            map_name,
            maxNegP,
            separate=False):  ## draw plot box
    global boxHeight

    mls = style.linewidth.THIck
    ldash = style.linestyle.dashed
    largeLab = text.size.huge
    yLabAtt = [largeLab, text.halign.boxright, text.valign.middle]
    ypLabAtt = [largeLab, text.halign.boxleft, text.valign.middle]

    boxY = exon_height + exon_height
    boxHeight = intron_width * 1.2
    bW = exon_width + intron_width

    second_box_offset = -(boxHeight + 2 * boxY)

    myMax = max_point_value * 1.01
    yLab = drawutils.make_label(myMax)
    ypLab = drawutils.make_label(maxNegP)

    units_per_bp = float(exon_width) / exon_length

    xLab_1_3 = [-exon_length, 0, intron_length / 2, intron_length]
    xLab_2_4 = [-intron_length, -intron_length / 2, 0, exon_length]

    xS = 0
    xE = xS + bW
    rect1, r1_line, r1_ss = drawutils.boxes(xS, bW, scale, boxY, boxHeight,
                                            exon_width)


    drawutils.draw_y_axis(c, scale, yLab, xS, boxY, boxHeight)

    x_axis_canvas = canvas.canvas()

    drawutils.draw_x_axis_segment(x_axis_canvas, scale, xS, boxY, xLab_1_3, units_per_bp)

    xS = xE + divider_gap
    xE = xS + bW
    rect2, r2_line, r2_ss = drawutils.boxes(xS, bW, scale, boxY, boxHeight,
                                            intron_width)

    drawutils.draw_x_axis_segment(x_axis_canvas, scale, xS, boxY, xLab_2_4, units_per_bp)

    drawutils.title_and_legend(c, scale, xE, divider_gap, boxY, boxHeight,
                               exon_width, intron_width, map_name, nu, nd, nb)

    xS = xE + divider_gap + divider_gap
    xE = xS + bW
    rect3, r3_line, r3_ss = drawutils.boxes(xS, bW, scale, boxY, boxHeight,
                                            exon_width)

    drawutils.draw_x_axis_segment(x_axis_canvas, scale, xS, boxY, xLab_1_3, units_per_bp)

    xS = xE + divider_gap
    xE = xS + bW
    rect4, r4_line, r4_ss = drawutils.boxes(xS, bW, scale, boxY, boxHeight,
                                            intron_width)

    drawutils.draw_x_axis_segment(x_axis_canvas, scale, xS, boxY, xLab_2_4, units_per_bp)

    c.insert(x_axis_canvas)
    if separate:
        c.insert(x_axis_canvas, [trafo.translate(0, -boxY * scale)])

    box_canvas = canvas.canvas()

    box_canvas.stroke(rect1, [mls])
    box_canvas.stroke(rect2, [mls])
    box_canvas.stroke(rect3, [mls])
    box_canvas.stroke(rect4, [mls])

    for i in range(3):
        box_canvas.stroke(r1_line[i], [ldash])
        box_canvas.stroke(r2_line[i], [ldash])
        box_canvas.stroke(r3_line[i], [ldash])
        box_canvas.stroke(r4_line[i], [ldash])

    box_canvas.stroke(r1_ss)
    box_canvas.stroke(r2_ss)
    box_canvas.stroke(r3_ss)
    box_canvas.stroke(r4_ss)

    c.insert(box_canvas)
    if separate:
        c.insert(box_canvas, [trafo.translate(0, second_box_offset * scale)])

    yp_y_offset = second_box_offset if separate else 0
    xE = xE + exon_width
    drawutils.draw_yp_axis(c, scale, ypLab, yp_y_offset, xE, boxY, boxHeight,
                           divider_gap)

    logging.debug("Done drawing plot box")


def plotRegions(
    tRegion_1,
    tRegion_2,
    xS,
    boxY,
    bH,
    fW,
    sW,
    mpv,
    scale,
    howmany=3
):  ## plot given regions, fW, sW tells the width of the first and the second region
    pup = []
    pdn = []
    pbg = []
    pvalup = []
    pvaldn = []
    y1 = []
    y2 = []

    for jj in range(len(tRegion_1) -
                    1):  ## for each point in the first region (upstream exon)
        x1 = xS + jj * float(fW) / len(tRegion_1) + (float(fW) /
                                                     len(tRegion_1)) / 2.0
        x2 = xS + (jj + 1) * float(fW) / len(tRegion_1) + (
            float(fW) / len(tRegion_1)) / 2.0
        if howmany == 3:  ## actual points
            y1 = [
                boxY + tRegion_1[jj][0] * bH / mpv,
                boxY + tRegion_1[jj][1] * bH / mpv,
                boxY + tRegion_1[jj][2] * bH / mpv
            ]
            y2 = [
                boxY + tRegion_1[jj + 1][0] * bH / mpv,
                boxY + tRegion_1[jj + 1][1] * bH / mpv,
                boxY + tRegion_1[jj + 1][2] * bH / mpv
            ]

            pup.append([x1 * scale, y1[0] * scale])
            pdn.append([x1 * scale, y1[1] * scale])
            pbg.append([x1 * scale, y1[2] * scale])

            if jj == len(
                    tRegion_1
            ) - 2:  ## second to the last element. add coords for the last element
                pup.append([x2 * scale, y2[0] * scale])
                pdn.append([x2 * scale, y2[1] * scale])
                pbg.append([x2 * scale, y2[2] * scale])

        elif howmany == 2:  ## pvalue points
            y1 = [
                boxY + tRegion_1[jj][0] * bH / mpv,
                boxY + tRegion_1[jj][1] * bH / mpv
            ]
            y2 = [
                boxY + tRegion_1[jj + 1][0] * bH / mpv,
                boxY + tRegion_1[jj + 1][1] * bH / mpv
            ]

            pvalup.append([x1 * scale, y1[0] * scale])
            pvaldn.append([x1 * scale, y1[1] * scale])

            if jj == len(
                    tRegion_1
            ) - 2:  ## second to the last element. add coords for the last element
                pvalup.append([x2 * scale, y2[0] * scale])
                pvaldn.append([x2 * scale, y2[1] * scale])

    xS = xS + fW
    for jj in range(len(tRegion_2) -
                    1):  ## for each point in the second region (intron)
        x1 = xS + jj * float(sW) / len(tRegion_2) + (float(sW) /
                                                     len(tRegion_2)) / 2.0
        x2 = xS + (jj + 1) * float(sW) / len(tRegion_2) + (
            float(sW) / len(tRegion_2)) / 2.0
        if howmany == 3:  ## actual points
            y1 = [
                boxY + tRegion_2[jj][0] * bH / mpv,
                boxY + tRegion_2[jj][1] * bH / mpv,
                boxY + tRegion_2[jj][2] * bH / mpv
            ]
            y2 = [
                boxY + tRegion_2[jj + 1][0] * bH / mpv,
                boxY + tRegion_2[jj + 1][1] * bH / mpv,
                boxY + tRegion_2[jj + 1][2] * bH / mpv
            ]

            pup.append([x1 * scale, y1[0] * scale])
            pdn.append([x1 * scale, y1[1] * scale])
            pbg.append([x1 * scale, y1[2] * scale])

            if jj == len(
                    tRegion_2
            ) - 2:  ## second to the last element. add coords for the last element
                pup.append([x2 * scale, y2[0] * scale])
                pdn.append([x2 * scale, y2[1] * scale])
                pbg.append([x2 * scale, y2[2] * scale])
        elif howmany == 2:  ## pvalue points
            y1 = [
                boxY + tRegion_2[jj][0] * bH / mpv,
                boxY + tRegion_2[jj][1] * bH / mpv
            ]
            y2 = [
                boxY + tRegion_2[jj + 1][0] * bH / mpv,
                boxY + tRegion_2[jj + 1][1] * bH / mpv
            ]

            pvalup.append([x1 * scale, y1[0] * scale])
            pvaldn.append([x1 * scale, y1[1] * scale])

            if jj == len(
                    tRegion_2
            ) - 2:  ## second to the last element. add coords for the last element
                pvalup.append([x2 * scale, y2[0] * scale])
                pvaldn.append([x2 * scale, y2[1] * scale])

    logging.debug("Done plotRegions function")
    if howmany == 3:
        return pup, pdn, pbg
    elif howmany == 2:
        return pvalup, pvaldn




def fillUpPath(tPoints):
    # AUDIT F15: unavailable windows form gaps, not fabricated density values.
    rPath = path.path()
    connected = False
    for x, y in tPoints:
        if not numpy.isfinite(x) or not numpy.isfinite(y):
            connected = False
            continue
        rPath.append(path.lineto(x, y) if connected else path.moveto(x, y))
        connected = True
    return rPath


def drawAcutalPlot(
    c,
    dp,
    mpv,
    eH,
    eW,
    iW,
    gap,
    sGap,
    scale,
    wl,
    sl,
    nPP,
    mNP,
    separate=False
):  ## draw plots now.. nPP is negativePvaluePoins, mNP is maximumNegPval
    global boxHeight

    mls = style.linewidth.THIck
    bgmls = style.linewidth.Thick
    ldash = style.linestyle.dashed
    largeLab = text.size.huge

    boxY = eH + eH
    bH = boxHeight
    bW = eW + iW
    second_box_offset = -(bH + 2 * boxY)
    p_transform = trafo.translate(0,
                                  second_box_offset * scale if separate else 0)

    upColor = color.rgb.red
    dnColor = color.rgb.blue
    bgColor = color.rgb.black

    xS = 0
    points_up, points_dn, points_bg = plotRegions(dp[0], dp[1], xS, boxY, bH,
                                                  eW, iW, mpv, scale, 3)
    path_up = fillUpPath(points_up)
    path_dn = fillUpPath(points_dn)
    path_bg = fillUpPath(points_bg)
    c.stroke(path_up, [upColor, mls])
    c.stroke(path_dn, [dnColor, mls])
    c.stroke(path_bg, [bgColor, bgmls])
    points_pup, points_pdn = plotRegions(nPP[0], nPP[1], xS, boxY, bH, eW, iW,
                                         mNP, scale, 2)
    path_pup = fillUpPath(points_pup)
    path_pdn = fillUpPath(points_pdn)
    c.stroke(path_pup, [upColor, ldash, p_transform])
    c.stroke(path_pdn, [dnColor, ldash, p_transform])

    xS = xS + eW + iW + gap
    points_up, points_dn, points_bg = plotRegions(dp[2], dp[3], xS, boxY, bH,
                                                  iW, eW, mpv, scale)
    path_up = fillUpPath(points_up)
    path_dn = fillUpPath(points_dn)
    path_bg = fillUpPath(points_bg)
    c.stroke(path_up, [upColor, mls])
    c.stroke(path_dn, [dnColor, mls])
    c.stroke(path_bg, [bgColor, bgmls])
    points_pup, points_pdn = plotRegions(nPP[2], nPP[3], xS, boxY, bH, iW, eW,
                                         mNP, scale, 2)
    path_pup = fillUpPath(points_pup)
    path_pdn = fillUpPath(points_pdn)
    c.stroke(path_pup, [upColor, ldash, p_transform])
    c.stroke(path_pdn, [dnColor, ldash, p_transform])

    xS = xS + iW + eW + gap + gap
    points_up, points_dn, points_bg = plotRegions(dp[4], dp[5], xS, boxY, bH,
                                                  eW, iW, mpv, scale)
    path_up = fillUpPath(points_up)
    path_dn = fillUpPath(points_dn)
    path_bg = fillUpPath(points_bg)
    c.stroke(path_up, [upColor, mls])
    c.stroke(path_dn, [dnColor, mls])
    c.stroke(path_bg, [bgColor, bgmls])
    points_pup, points_pdn = plotRegions(nPP[4], nPP[5], xS, boxY, bH, eW, iW,
                                         mNP, scale, 2)
    path_pup = fillUpPath(points_pup)
    path_pdn = fillUpPath(points_pdn)
    c.stroke(path_pup, [upColor, ldash, p_transform])
    c.stroke(path_pdn, [dnColor, ldash, p_transform])

    xS = xS + eW + iW + gap
    points_up, points_dn, points_bg = plotRegions(dp[6], dp[7], xS, boxY, bH,
                                                  iW, eW, mpv, scale)
    path_up = fillUpPath(points_up)
    path_dn = fillUpPath(points_dn)
    path_bg = fillUpPath(points_bg)
    c.stroke(path_up, [upColor, mls])
    c.stroke(path_dn, [dnColor, mls])
    c.stroke(path_bg, [bgColor, bgmls])
    points_pup, points_pdn = plotRegions(nPP[6], nPP[7], xS, boxY, bH, iW, eW,
                                         mNP, scale, 2)
    path_pup = fillUpPath(points_pup)
    path_pdn = fillUpPath(points_pdn)
    c.stroke(path_pup, [upColor, ldash, p_transform])
    c.stroke(path_pdn, [dnColor, ldash, p_transform])

    logging.debug("Done drawing acutal plots")


def countPerWindow(cpw, ic, sign, rLen):
    global wLen, sLen
    sInd = 0
    cInd = 0
    if sign == -1:  ## need to to this in reverse way
        sInd = -rLen
        cInd = -1
    for k in range(len(cpw)):  ## do it this many times
        for mm in motifs:  ## examine all motifs
            for n in range(sInd, sInd + sign * (wLen), sign):
                cpw[k] += ic[mm][n][cInd]
        sInd += sLen


def ccc(ic):
    # AUDIT F15: density is the proportion of eligible events, without a cap.
    return {index: [binary_density(ic[index][locus]) for locus in range(len(ic[index]))]
            for index in range(8)}

def _plot_snapshot(paths):
    return {path: (os.stat(path).st_mtime_ns, os.stat(path).st_size)
            if os.path.isfile(path) else None for path in paths}


def _record_plot_writes(before):
    # AUDIT F8: stale plots from overwritten runs must not enter this manifest.
    after = _plot_snapshot(before)
    run_outputs.extend(path for path, status in after.items()
                       if status is not None and status != before[path])


def plotMotifs_finale(mapName,
                      maxPointValue,
                      maxNegPval,
                      drawPoints,
                      negPvalPoints,
                      separate=False):

    exonHeight = 30
    exonWidth = eLen
    intronWidth = iLen
    indentation = 5
    dividerGap = 12
    slopeGap = 20
    Scale = 0.03
    canv = canvas.canvas()

    drawNode(canv, exonHeight, exonWidth, intronWidth, indentation, dividerGap,
             slopeGap, Scale)
    drawBox(canv, maxPointValue, exonHeight, exonWidth, intronWidth,
            indentation, dividerGap, slopeGap, Scale, eLen, iLen, mapName,
            maxNegPval, separate)
    drawAcutalPlot(canv, drawPoints, maxPointValue, exonHeight, exonWidth,
                   intronWidth, dividerGap, slopeGap, Scale, wLen, sLen,
                   negPvalPoints, maxNegPval, separate)

    safe_map_name = re.sub(r'[^A-Za-z0-9._-]+', '_', mapName)
    pdf_path = mapsPath + '/' + 'SE.' + safe_map_name + '.pdf'
    png_path = mapsPath + '/' + safe_map_name + '.png'
    before = _plot_snapshot((pdf_path, png_path))
    pdf_ok, png_ok = drawutils.export_canvas_outputs(
        canv, pdf_path, png_path, png_resolution=100, logger=logging)
    _record_plot_writes(before)
    if not pdf_ok:
        logging.debug("PDF export failed for %s", mapName)
    if not png_ok:
        logging.debug("PNG export failed for %s", mapName)


def plotMotifs(
    d, mc, mapName, pdic_up, pdic_dn, motifs,
    cdist):  ## plotting motifs (count, name, pvalue.up.vs.bg, pvalue.dn.vs.bg)

    global uNum, dNum, bNum
    upc = {}
    dnc = {}
    bgc = {}

    upc = ccc(cdist['up'])
    dnc = ccc(cdist['dn'])
    bgc = ccc(cdist['bg'])

    drawPoints = []
    negPvalPoints = []
    maxPointValue = 0.0000000000000000000000000000000000000000000000000001
    maxNegPval = 0.00000000000000000000000000000000000000000000000000000000000000001
    for zz in range(8):  ## for 8 regions
        drawPoints.append([])
        negPvalPoints.append([])
        for ind in range(
                len(upc[zz])
        ):  ## everything has the same items. It's okay to use this.
            drawPoints[zz].append([
                upc[zz][ind],
                dnc[zz][ind],
                bgc[zz][ind]
            ])
            negPvalPoints[zz].append([
                # AUDIT F14: retain underflowed p=0 in plots without masking it as unavailable.
                -numpy.log10(max(pdic_up[zz][ind], numpy.nextafter(0.0, 1.0))),
                -numpy.log10(max(pdic_dn[zz][ind], numpy.nextafter(0.0, 1.0)))
            ])
            for yy in drawPoints[zz][ind]:  ## examine three values
                if yy > maxPointValue:  ## new max here
                    maxPointValue = yy
            for yy in negPvalPoints[zz][ind]:  ## examine two values
                if yy > maxNegPval:  ## new max here
                    maxNegPval = yy
    logging.debug("Max average height is: %f" % maxPointValue)
    logging.debug("Max negPval is: %f" % maxNegPval)
    val = [mapName, maxPointValue, maxNegPval, drawPoints, negPvalPoints]
    d.append(val)


def printCountDist(cdist, tName, tMotif):
    # AUDIT F6: expose per-window eligible denominators alongside binary values.
    for label in ('up', 'dn', 'bg'):
        filename = os.path.join(tempPath, tName + '.' + tMotif + '.countDist.' + label + '.txt')
        with open(filename, 'w') as handle:
            handle.write('Region\tposition\tsum\teligible\tdensity\tvalues\n')
            for index, name in enumerate(REGION_NAMES):
                for locus, values in cdist[label][index].items():
                    # AUDIT R4: binary histogram replaces redundant dense text values.
                    density = binary_density(values)
                    handle.write(f'{name}\t{locus * sLen}\t{values.hits}\t{values.eligible}\t' +
                                 ('NA' if not numpy.isfinite(density) else str(density)) + '\t' +
                                 f'0:{values.eligible-values.hits},1:{values.hits},NA:{values.total-values.eligible}\n')
        run_outputs.append(filename)

def computePValues(cdist_one, cdist_two,
                     test_p):  ## count p value for one vs. two
    rName = {
        0: 'UpstreamExon_3prime',
        1: 'UpstreamExonIntron',
        2: 'UpstreamIntron',
        3: 'TargetExon_5prime',
        4: 'TargetExon-3prime',
        5: 'DownstreamIntron',
        6: 'DownstreamExonIntron',
        7: 'DownstreamExon_5prime'
    }
    rng = numpy.random.default_rng(stat_seed)
    fisher_cache = {}
    reasons = {}
    globals().setdefault("stat_reasons", {})[id(test_p)] = reasons

    for zz in range(8):  ## for 8 regions
        test_p[zz] = {}
        for locus in range(len(cdist_one[zz])):
            # AUDIT F6: ineligible events leave both numerator and denominator.
            one = cdist_one[zz][locus]
            two = cdist_two[zz][locus]
            if not one.eligible or not two.eligible:
                test_p[zz][locus] = float('nan')
                reasons[(zz, locus)] = 'no eligible events in one or both sets'
                continue
            try:
                if stat_method != 'fisher':
                    first, second = one.observations(), two.observations()
                if stat_method == 'mannwhitney_greater':
                    pvalue = float(stats.mannwhitneyu(first, second, alternative='greater')[1])
                elif stat_method == 'brunnermunzel_greater':
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore', category=RuntimeWarning, module=r'scipy\\.stats.*')
                        pvalue = float(stats.brunnermunzel(first, second, alternative='greater', distribution='normal')[1])
                elif stat_method == 'permutation_one_sided':
                    first_arr = numpy.asarray(first, dtype=float)
                    second_arr = numpy.asarray(second, dtype=float)
                    n_first = first_arr.size
                    observed = float(numpy.mean(first_arr) - numpy.mean(second_arr))
                    pooled = numpy.concatenate([first_arr, second_arr])
                    ge_count = 0
                    for _ in range(stat_permutations):
                        perm = rng.permutation(pooled)
                        stat = float(numpy.mean(perm[:n_first]) - numpy.mean(perm[n_first:]))
                        if stat >= observed:
                            ge_count += 1
                    pvalue = (ge_count + 1.0) / (stat_permutations + 1.0)
                else:
                    # AUDIT F1: Fisher uses eligible binary events, not motif occurrences.
                    # AUDIT R4: repeated binary tables share an exact Fisher result.
                    key = (one.hits, one.eligible, two.hits, two.eligible)
                    if key not in fisher_cache:
                        fisher_cache[key] = float(stats.fisher_exact(
                            [[one.hits, one.eligible-one.hits], [two.hits, two.eligible-two.hits]],
                            alternative=args.fisher_alternative)[1])
                    pvalue = fisher_cache[key]
            except Exception as exc:
                # AUDIT F14: statistical exceptions must retain diagnostic context.
                raise RuntimeError(f'{stat_method}: {REGION_NAMES[zz]} window {locus * sLen}: {exc}') from exc

            if not numpy.isfinite(pvalue):
                reasons[(zz, locus)] = 'nonfinite statistical result'
                pvalue = float('nan')
            test_p[zz][locus] = pvalue
    return test_p


def printPval(pdic, dFile, eNum):
    # AUDIT F14: unavailable tests are NA with an explicit reason, never p=1.
    header = pvalue_header_label(stat_method)
    dFile.write('Region\tposition\t' + header + '\treason\n')
    reasons = globals().get('stat_reasons', {}).get(id(pdic), {})
    for index, name in enumerate(REGION_NAMES):
        for locus, pvalue in pdic[index].items():
            value = str(pvalue) if numpy.isfinite(pvalue) else 'NA'
            dFile.write(f'{name}\t{locus * sLen}\t{value}\t{reasons.get((index, locus), "")}\n')

def makeIndividualMaps(d, line):
    global mCount, totalExonCount

    test_up1 = {}
    test_dn1 = {}
    ele = line.strip().split('\t')
    tName = ele[0]
    tMotif = ele[1]

    tmpFile = open(tempPath + '/' + tName + '.' + tMotif + '.txt', 'w')
    tmpFile.write(tHeader + '\n')
    tmpFile.write(tName + '\t' + tMotif + '\n')
    tmpFile.close()
    run_outputs.append(tempPath + '/' + tName + '.' + tMotif + '.txt')
    tmpFile = open(tempPath + '/' + tName + '.' + tMotif + '.txt')
    mCount = {}
    motifs = initMotif(tmpFile, mCount)
    cdist = countMotif(mCount, tName, tMotif, motifs)
    printCountDist(cdist, tName, tMotif)
    test_up = computePValues(cdist['up'], cdist['bg'], test_up1)
    test_dn = computePValues(cdist['dn'], cdist['bg'], test_dn1)
    upbgPFile = open(
        tempPath + '/' + tName + '.' + tMotif + '.pVal.up.vs.bg.txt', 'w')
    dnbgPFile = open(
        tempPath + '/' + tName + '.' + tMotif + '.pVal.dn.vs.bg.txt', 'w')
    printPval(test_up, upbgPFile, totalExonCount['up'])
    printPval(test_dn, dnbgPFile, totalExonCount['dn'])
    plotMotifs(d, mCount, tName + '-' + tMotif, test_up, test_dn, motifs,
               cdist)
    upbgPFile.close()
    dnbgPFile.close()
    run_outputs.extend([upbgPFile.name, dnbgPFile.name])
    tmpFile.close()



def wccount(filename):
    count = 0
    with open(filename, "rb") as f:
        for _ in f:
            count += 1
    return count




def minPvalueOut(exonType):
    # AUDIT F8: aggregate only positional tables produced for current-run motifs.
    suffix = 'up' if exonType == 'up' else 'dn'
    rows = []
    for motif_id in run_motif_ids:
        filename = os.path.join(tempPath, motif_id + '.pVal.' + suffix + '.vs.bg.txt')
        buckets = {name: [] for name in REGION_NAMES}
        with open(filename) as handle:
            next(handle)
            for line in handle:
                fields = line.rstrip('\n').split('\t')
                if fields[2] != 'NA':
                    buckets[fields[0]].append(float(fields[2]))
        # AUDIT F14: a region with no usable test remains explicitly unavailable.
        rows.append([motif_id] + [min(buckets[name]) if buckets[name] else None for name in REGION_NAMES])
    rows.sort(key=lambda row: (row[2] is None, row[2] if row[2] is not None else 0))
    filename = os.path.join(outPath, 'pVal.' + suffix + '.vs.bg.RNAmap.txt')
    with open(filename, 'w') as handle:
        handle.write('RBP\tsmallest_p_in_upstreamExon-3prime\tsmallest_p_in_upstreamExonIntron\tsmallest_p_in_upstreamIntron\tsmallest_p_in_targetExon-5prime\tsmallest_p_in_targetExon-3prime\tsmallest_p_in_downstreamIntron\tsmallest_p_in_downstreamExonIntron\tsmallest_p_in_downstreamExon-5prime\n')
        for row in rows:
            handle.write('\t'.join('NA' if value is None else str(value) for value in row) + '\n')
    run_outputs.append(filename)

def _build_worker_state():
    state = {}
    for key, value in globals().items():
        # AUDIT R4: shared sequence arrays reopen read-only instead of being pickled.
        if key in {'region_data', 'event_sets', 'sequence_arrays', 'sequence_origins',
                   'sequence_lengths', 'sequence_labels', 'genome'}:
            continue
        if key.startswith('__'):
            continue
        if callable(value):
            continue
        if isinstance(value, types.ModuleType):
            continue
        try:
            pickle.dumps(value)
        except Exception:
            continue
        state[key] = value
    return state


def _init_worker(state):
    globals().update(state)
    global sequence_arrays, sequence_origins, sequence_lengths, sequence_labels
    sequence_arrays = [numpy.load(filename, mmap_mode='r', allow_pickle=False) for filename in sequence_paths]
    with numpy.load(os.path.join(tempPath, 'sequence_metadata.npz'), allow_pickle=False) as metadata:
        sequence_origins = metadata['origins']
        sequence_lengths = metadata['lengths']
        sequence_labels = metadata['labels']
    for array in (sequence_origins, sequence_lengths, sequence_labels):
        array.flags.writeable = False


def _make_individual_map_worker(line):
    global run_outputs
    run_outputs = []
    wall_start, cpu_start = timeit.default_timer(), time.process_time()
    d = []
    makeIndividualMaps(d, line)
    return (d[0] if d else None, run_outputs,
            timeit.default_timer() - wall_start, time.process_time() - cpu_start)


def maybe_plot_motif_map(mapname, maxPointValue, maxNegPval, drawPoints,
                         negPvalPoints, separate=False, out_dir=None, event_type=None):
    def write_fallback(reason):
        logging.warning("%s; attempting motif fallback plot for %s", reason,
                        mapname)
        if out_dir and event_type:
            safe_map_name = re.sub(r'[^A-Za-z0-9._-]+', '_', mapname)
            pdf_path = os.path.join(out_dir, 'maps',
                                    f'{event_type}.{safe_map_name}.pdf')
            png_path = os.path.join(out_dir, 'maps',
                                    f'{event_type}.{safe_map_name}.png')
            before = _plot_snapshot((pdf_path, png_path))
            try:
                pdf_ok, png_ok = drawutils.export_motif_map_fallback(
                    event_type, mapname, drawPoints, negPvalPoints,
                    maxPointValue, maxNegPval, pdf_path, png_path,
                    logger=logging)
                _record_plot_writes(before)
                if pdf_ok or png_ok:
                    logging.info(
                        "Motif fallback plot generated for %s (PDF: %s, PNG: %s)",
                        mapname, pdf_ok, png_ok)
                else:
                    logging.warning(
                        "Motif fallback plot failed for %s; text outputs only",
                        mapname)
            except Exception as exc:
                logging.warning(
                    "Exception in motif fallback plot for %s: %s", mapname, exc)

    if not DRAW_MOTIF_PLOTS:
        write_fallback("PyX text fonts unavailable")
        return
    try:
        plotMotifs_finale(mapname, maxPointValue, maxNegPval, drawPoints,
                          negPvalPoints, separate)
    except Exception as exc:
        logging.warning("Could not draw motif plot for %s: %s", mapname, exc)
        drawutils.suppress_pyx_text_cleanup(logging)
        write_fallback("PyX motif plot failed")


def run_individual_map_workers(lines):
    global run_outputs
    lines = [line for line in lines if line.strip()]
    if not lines:
        return []
    # AUDIT F8: summaries use only motif identifiers produced by this run.
    motif_ids = ['.'.join(line.strip().split('\t')[:2]) for line in lines]
    if len(set(motif_ids)) != len(motif_ids) or set(motif_ids) & set(run_motif_ids):
        raise ValueError('Duplicate RBP/motif output identifier in motif inputs')
    run_motif_ids.extend(motif_ids)
    saved_outputs = run_outputs
    worker_state = _build_worker_state()
    # AUDIT F20: use the requested allocation, bounded by the actual motif count.
    worker_count = min(args.workers, len(lines))
    if worker_count == 1:
        results = [_make_individual_map_worker(line) for line in lines]
    else:
        try:
            with multiprocessing.get_context("spawn").Pool(
                    processes=worker_count,
                    initializer=_init_worker,
                    initargs=(worker_state, )) as pool:
                results = pool.map(_make_individual_map_worker, lines)
        except PermissionError:
            logging.warning(
                "Multiprocessing pool unavailable; falling back to serial motif map generation"
            )
            _init_worker(worker_state)
            results = [_make_individual_map_worker(line) for line in lines]
    # AUDIT R4: separate measured motif task cost from spawn/rendering overhead.
    args.motif_task_wall_seconds = getattr(args, 'motif_task_wall_seconds', 0) + sum(item[2] for item in results)
    args.motif_task_cpu_seconds = getattr(args, 'motif_task_cpu_seconds', 0) + sum(item[3] for item in results)
    run_outputs = saved_outputs + [path for item in results for path in item[1]]
    return [item[0] for item in results if item[0] is not None]


def run_pipeline():
    global genome, start, uNum, dNum, bNum, tHeader, event_sets, region_data, totalExonCount, motif_padding
    # AUDIT R4: validate maximum match width before any sequence extraction.
    patterns = []
    for motif_path in (args.knownMotifs, args.motif):
        if motif_path != 'NA':
            with open(motif_path) as handle:
                next(handle)
                patterns.extend(line.strip().split('\t')[1] for line in handle if line.strip())
    motif_padding = max((max_motif_width(pattern) - 1 for pattern in patterns), default=0)
    if max(iLen, eLen) + motif_padding > numpy.iinfo(numpy.int16).max:
        raise ValueError('Region length plus motif overhang exceeds int16 positional coordinates')
    extraction_start = timeit.default_timer()
    logging.debug("================================")
    logging.debug("GETTING GENOME FASTA OBJECT")
    genome = load_genome(args.genome, fasta_root_PATH)
    logging.debug("DONE GETTING GENOME FASTA OBJECT")
    logging.debug("================================")

    logging.debug("================================")
    logging.debug("MAKING INPUT FILES FROM rMATS FILE")
    try:
        makeInputFiles(rMATS)
    except:
        logging.debug("There is an exception in making input from rMATS output")
        logging.debug("Exception: %s" % sys.exc_info()[0])
        logging.debug("Detail: %s" % sys.exc_info()[1])
        sys.exit(-1)
    logging.debug("DONE MAKING INPUT FILES FROM rMATS FILE")
    logging.debug("================================")

    # AUDIT F13: all input paths converge on one validated eight-column schema.
    coord_paths = {label: os.path.join(exonPath, label + '.coord.txt') for label in ('up', 'dn', 'bg')}
    event_sets = read_event_sets(coord_paths, genome, allow_overlap=args.allow_overlap)
    totalExonCount = {label: len(events) for label, events in event_sets.items()}
    run_outputs.extend(coord_paths.values())
    region_data = {}

    logging.debug("================================")
    logging.debug("MAKING FASTA FILES")
    try:
        getFasta('up')
        getFasta('dn')
        getFasta('bg')
        prepare_sequence_store()
        for jj in ['up', 'dn', 'bg']:
            rg = region[0]
            fFile = open(fastaPath + '/' + jj + '.' + rg + '.fasta')
            c = 0
            for dummy in fFile:
                noUse = next(fFile).strip()
                c += 1

            if jj == 'up':
                uNum = c
            elif jj == 'dn':
                dNum = c
            elif jj == 'bg':
                bNum = c
                fFile.close()

        logging.debug(
            "Number of events for upregulated, downregulated, and background: %d, %d, %d"
            % (uNum, dNum, bNum))
        pass

    except:
        logging.debug("There is an exception in making fasta files")
        logging.debug("Exception: %s" % sys.exc_info()[0])
        logging.debug("Detail: %s" % sys.exc_info()[1])
        sys.exit(-2)
    logging.debug("DONE MAKING FASTA FILES")
    args.extraction_seconds = timeit.default_timer() - extraction_start
    logging.debug("================================")




    logging.debug("================================")
    logging.debug("MAKING INDIVIDUAL MAPS")
    try:
        kFile.seek(0)
        tHeader = kFile.readline().strip()

        known_lines = [line for line in kFile]
        motif_start = timeit.default_timer()
        d = run_individual_map_workers(known_lines)
        args.motif_seconds = timeit.default_timer() - motif_start

        for i in range(len(d)):

            [mapname, maxPointValue, maxNegPval, drawPoints, negPvalPoints] = d[i]
            maybe_plot_motif_map(mapname, maxPointValue, maxNegPval, drawPoints,
                                 negPvalPoints, args.separate, outPath, 'SE')

        stop = timeit.default_timer()
        print(stop - start)
        if args.motif != 'NA':
            mFile.seek(0)
            tHeader = mFile.readline().strip()

            motif_lines = [line2 for line2 in mFile]
            motif_start = timeit.default_timer()
            d2 = run_individual_map_workers(motif_lines)
            args.motif_seconds += timeit.default_timer() - motif_start

            for i in range(len(d2)):

                [mapname, maxPointValue, maxNegPval, drawPoints,
                 negPvalPoints] = d2[i]
                maybe_plot_motif_map(mapname, maxPointValue, maxNegPval,
                                     drawPoints, negPvalPoints, False, outPath, 'SE')

        pass

    except:
        logging.debug("There is an exception in making individual maps")
        logging.exception("Exception while making individual maps")
        sys.exit(-3)
    logging.debug("DONE MAKING INDIVIDUAL MAPS")
    logging.debug("================================")

    logging.debug("================================")
    logging.debug("SORTING BBPs by minimum P-values")
    logging.debug("processing up regulated exons..")
    minPvalueOut("up")
    logging.debug("Done processing up regulated exons..")
    logging.debug("processing dn regulated exons..")
    minPvalueOut("down")
    logging.debug("Done processing dn regulated exons..")
    logging.debug("DONE SORTING BBPs by minimum P-values")
    logging.debug("================================")

    if args.motif != "NA":
        mFile.close()

    # AUDIT R4: release Windows mmap handles before registered temporary cleanup.
    for array in sequence_arrays:
        array._mmap.close()

    logging.debug("Program ended")
    currentTime = time.time()
    runningTime = currentTime - startTime
    logging.debug("Program ran %.2d:%.2d:%.2d" %
                  (runningTime / 3600,
                   (runningTime % 3600) / 60, runningTime % 60))

    # AUDIT R2/R3: hash every registered summary input before deleting current temp files.
    logging.shutdown()
    parameters = dict(vars(args), stat_method=stat_method,
                      stat_permutations=stat_permutations, stat_seed=stat_seed)
    completed = [Path(filename) for filename in run_outputs]
    temporary = [path for path in completed if path.resolve().parent == Path(tempPath).resolve()]
    write_run_manifest(outPath, completed, parameters,
                       delete_files=temporary if args.delete_temp else [])
    if args.delete_temp and not any(Path(tempPath).iterdir()):
        Path(tempPath).rmdir()

    sys.exit(0)




def main():
    setup_runtime()
    run_pipeline()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()


