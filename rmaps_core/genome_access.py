from pathlib import Path
from functools import lru_cache

from pyfaidx import Fasta


_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


@lru_cache(maxsize=None)
def load_genome(build: str, base_dir: str) -> Fasta:
    """
    Load a genome FASTA for the given build from a genomedata-style layout:
    base_dir/<build>/<build>.fa
    """
    base = Path(base_dir)
    fasta_path = base / build / f"{build}.fa"
    return Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True)


def fetch_seq(fasta: Fasta, strand: str, chrom: str, start: int, end: int) -> str:
    """
    Fetch sequence from a pyfaidx Fasta:
    - 0-based, end-exclusive slicing
    - reverse-complement on '-' strand
    - clip chromosome edges and pad in genomic orientation before reverse complement
    - missing chromosomes and fetch failures raise with coordinate context
    """
    length = end - start
    if length <= 0:
        return ""
    # AUDIT F5: fail explicitly and pad before orientation conversion.
    try:
        chromosome_length = len(fasta[chrom])
        lo = min(chromosome_length, max(0, start))
        hi = min(chromosome_length, max(0, end))
        seq = str(fasta[chrom][lo:hi])
        if len(seq) != hi - lo:
            raise ValueError("FASTA returned a truncated sequence")
        left_pad = min(length, max(0, -start))
        right_pad = length - left_pad - len(seq)
        seq = "N" * left_pad + seq + "N" * right_pad
        if strand == "-":
            seq = revcomp(seq)
        return seq.upper()
    except Exception as exc:
        raise ValueError(f"Cannot fetch {chrom}:{start}-{end} ({strand}): {exc}") from exc

