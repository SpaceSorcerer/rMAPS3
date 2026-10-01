"""Validated SE events and binary, transcript-oriented window measurements."""

from dataclasses import dataclass
from pathlib import Path
import re
import warnings
from functools import lru_cache

from rmaps_core.genome_access import fetch_seq


HEADER = ("chr", "strand", "exonStart", "exonEnd", "firstExonStart",
          "firstExonEnd", "secondExonStart", "secondExonEnd")
REGION_NAMES = ("UpstreamExon_3prime", "UpstreamExonIntron", "UpstreamIntron",
                "TargetExon_5prime", "TargetExon-3prime", "DownstreamIntron",
                "DownstreamExonIntron", "DownstreamExon_5prime")


@dataclass(frozen=True)
class Event:
    chrom: str
    strand: str
    exon_start: int
    exon_end: int
    first_start: int
    first_end: int
    second_start: int
    second_end: int

    @property
    def event_id(self):
        return ":".join(str(value) for value in self.__dict__.values())


@dataclass(frozen=True)
class RegionSequence:
    sequence: str
    origin: int
    length: int


def read_event_sets(paths, genome, allow_overlap=False):
    # AUDIT F13: validate the complete eight-column identity before extraction.
    groups = {}
    for label, path in paths.items():
        events = []
        seen = set()
        with Path(path).open(encoding="utf-8-sig") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            if tuple(header) != HEADER:
                raise ValueError(f"{label}: {path}: required header is {' '.join(HEADER)}")
            for line_no, line in enumerate(handle, 2):
                fields = line.rstrip("\r\n").split("\t")
                context = f"{label}: {path}:{line_no}: event {':'.join(fields)}"
                if len(fields) != 8:
                    raise ValueError(f"{context}: expected 8 columns")
                try:
                    event = Event(fields[0], fields[1], *(int(v) for v in fields[2:]))
                except ValueError as exc:
                    raise ValueError(f"{context}: coordinates must be integers") from exc
                if event.strand not in {"+", "-"}:
                    raise ValueError(f"{context}: strand must be + or -")
                if event.chrom not in genome:
                    raise ValueError(f"{context}: chromosome absent from FASTA")
                pairs = ((event.first_start, event.first_end),
                         (event.exon_start, event.exon_end),
                         (event.second_start, event.second_end))
                if any(start < 0 or start >= end for start, end in pairs):
                    raise ValueError(f"{context}: require nonnegative start < end")
                if not (event.first_end <= event.exon_start and
                        event.exon_end <= event.second_start):
                    raise ValueError(f"{context}: require genomic-left, target, genomic-right exon order")
                if event in seen:
                    raise ValueError(f"{context}: duplicate event within {label}")
                seen.add(event)
                events.append(event)
        if not events:
            raise ValueError(f"{label}: {path}: empty event set")
        print(f"{label}: {len(events)} events")
        groups[label] = events
    labels = list(groups)
    for i, label in enumerate(labels):
        for other in labels[i + 1:]:
            common = set(groups[label]) & set(groups[other])
            if common:
                message = f"{label}/{other}: {len(common)} overlapping events; example {sorted(e.event_id for e in common)[0]}"
                if not allow_overlap:
                    raise ValueError(message + "; use --allow-overlap to permit")
                warnings.warn(message, stacklevel=2)
    return groups


@lru_cache(maxsize=512)
def _compiled_motif(pattern):
    pattern = pattern.replace("U", "T").replace("u", "t")
    return re.compile(pattern, re.IGNORECASE)


def max_motif_width(pattern):
    # AUDIT R4: bounded extraction requires a finite, context-free match width.
    regex = _compiled_motif(pattern)
    parsed = re._parser.parse(regex.pattern, regex.flags)
    def check(items):
        for opcode, value in items:
            if str(opcode) in {"ASSERT", "ASSERT_NOT", "AT"}:
                raise ValueError("Bounded motif scans do not support anchors or lookarounds")
            if str(opcode) == "SUBPATTERN":
                check(value[-1])
            elif str(opcode) == "BRANCH":
                for branch in value[1]:
                    check(branch)
            elif str(opcode) in {"MAX_REPEAT", "MIN_REPEAT", "POSSESSIVE_REPEAT"}:
                check(value[-1])
            elif str(opcode) == "ATOMIC_GROUP":
                check(value)
            elif str(opcode) == "GROUPREF_EXISTS":
                check(value[1])
                if value[2] is not None:
                    check(value[2])
    check(parsed)
    minimum, maximum = parsed.getwidth()
    if minimum == 0 or maximum >= re._constants.MAXREPEAT:
        raise ValueError(f"Motif requires a finite, positive match width: {pattern!r}")
    return maximum


def overlapping_hits(pattern, sequence):
    # AUDIT R4: native search skips nonmatching starts; next start preserves overlaps.
    regex = _compiled_motif(pattern)
    result = []
    position = 0
    while position <= len(sequence):
        match = regex.search(sequence, position)
        if match is None:
            break
        start, end = match.span()
        if start == end:
            raise ValueError(f"Motif regex must not match an empty sequence: {pattern!r}")
        result.append((start, end))
        position = start + 1
    return result


def event_regions(genome, event, intron, exon, padding=0):
    # AUDIT F6: true exon/intron boundaries and chromosome clipping define eligibility.
    if intron <= 0 or exon <= 0:
        raise ValueError("intron and exon lengths must be positive")
    if event.chrom not in genome:
        raise ValueError(f"event {event.event_id}: chromosome absent from FASTA")
    chromosome_length = len(genome[event.chrom])
    features = [(event.first_start, event.first_end),
                (event.first_end, event.exon_start),
                (event.exon_start, event.exon_end),
                (event.exon_end, event.second_start),
                (event.second_start, event.second_end)]
    if event.strand == "-":
        features.reverse()
    specs = ((0, exon, True), (1, intron, False), (1, intron, True),
             (2, exon, False), (2, exon, True), (3, intron, False),
             (3, intron, True), (4, exon, False))
    result = []
    for feature, length, end_anchored in specs:
        start, end = features[feature]
        feature_length = end - start
        anchor = feature_length - length if end_anchored else 0
        # AUDIT R4: fetch only the region plus finite motif overhang within its feature.
        lo, hi = max(0, anchor - padding), min(feature_length, anchor + length + padding)
        left, right = ((start + lo, start + hi) if event.strand == "+"
                       else (end - hi, end - lo))
        left = min(chromosome_length, max(0, left))
        right = min(chromosome_length, max(0, right))
        try:
            sequence = fetch_seq(genome, event.strand, event.chrom, left, right)
        except Exception as exc:
            raise ValueError(f"event {event.event_id}: fetch failed for [{left}, {right}): {exc}") from exc
        if len(sequence) != right - left:
            raise ValueError(f"event {event.event_id}: incomplete FASTA fetch [{left}, {right})")
        offset = left - start if event.strand == "+" else end - right
        result.append(RegionSequence(sequence, offset - anchor, length))
    return result


def binary_windows(region, patterns, window, step=1):
    """Window j covers [j, j+window) in transcript orientation in every region.

    Start positions are range(0, region.length-window+1, step). A hit overlaps when any
    nucleotide intersects the window. A window is eligible only when all its
    nucleotides belong to the event's feature and are present in the chromosome.
    End-anchored regions use j=0 at feature_end-region.length; other regions use
    j=0 at feature_start. None denotes ineligible, never a motif-negative event.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window and step must be positive")
    # AUDIT F2: end-anchored and start-anchored regions are counted independently.
    hits = [(start + region.origin, end + region.origin)
            for pattern in patterns for start, end in overlapping_hits(pattern, region.sequence)]
    # AUDIT F3: one half-open convention, independent of strand and motif length.
    result = []
    # AUDIT R1: complete windows only, including one window when length == window.
    for j in range(0, region.length - window + 1, step):
        if j < region.origin or j + window > region.origin + len(region.sequence):
            result.append(None)
        else:
            # AUDIT F1: the observational unit is one binary event per window.
            result.append(int(any(start < j + window and end > j for start, end in hits)))
    return result


def binary_table(first, second):
    rows = []
    for values in (first, second):
        eligible = [value for value in values if value is not None]
        if any(value not in (0, 1) for value in eligible):
            raise ValueError("Event/window counts must be binary 0/1 or None")
        rows.append([sum(eligible), len(eligible) - sum(eligible)])
    return rows


def binary_density(values):
    # AUDIT F15: event proportion uses the same eligible binary observations as Fisher.
    if isinstance(values, WindowCounts):
        return values.hits / values.eligible if values.eligible else float("nan")
    eligible = [value for value in values if value is not None]
    return sum(eligible) / len(eligible) if eligible else float("nan")


@dataclass(frozen=True)
class WindowCounts:
    # AUDIT R4: sufficient binary counts avoid n_events x n_windows Python objects.
    hits: int
    eligible: int
    total: int

    def observations(self):
        return [1] * self.hits + [0] * (self.eligible - self.hits)
