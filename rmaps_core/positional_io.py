"""Sparse positional schema 2; coordinates are transcript-oriented, half-open."""

import numpy as np
from contextlib import ExitStack
from zipfile import ZipFile


# AUDIT S2: reject incompatible sparse archives before interpreting their fields.
def _validate_schema(data):
    if "schema_version" not in data:
        raise ValueError("Unsupported sparse schema_version: missing; expected 2")
    version = np.asarray(data["schema_version"])
    if version.shape != () or version.dtype.kind not in "iu" or version.item() != 2:
        raise ValueError(f"Unsupported sparse schema_version: {version!r}; expected 2")


def iter_hits(path):
    # AUDIT R4: stream hit tuples without allocating an event-by-window matrix.
    with np.load(path, allow_pickle=False) as data:
        _validate_schema(data)
    with ExitStack() as stack:
        archive = stack.enter_context(ZipFile(path))
        streams, dtypes, sizes = [], [], []
        for key in ("event_index", "hit_start", "hit_end"):
            stream = stack.enter_context(archive.open(key + ".npy"))
            version = np.lib.format.read_magic(stream)
            reader = (np.lib.format.read_array_header_1_0 if version == (1, 0)
                      else np.lib.format.read_array_header_2_0)
            shape, _, dtype = reader(stream)
            if len(shape) != 1 or dtype.hasobject:
                raise ValueError("Sparse hit columns must be one-dimensional numeric arrays")
            streams.append(stream)
            dtypes.append(dtype)
            sizes.append(shape[0])
        if len(set(sizes)) != 1:
            raise ValueError("Sparse hit columns have unequal lengths")
        for offset in range(0, sizes[0], 65536):
            count = min(65536, sizes[0] - offset)
            chunks = [np.frombuffer(stream.read(count * dtype.itemsize), dtype=dtype, count=count)
                      for stream, dtype in zip(streams, dtypes)]
            for event, start, end in zip(*chunks):
                yield int(event), int(start), int(end)


def load_hits(path):
    # AUDIT R1/R4: reconstruct only complete windows from interval overlaps.
    with np.load(path, allow_pickle=False) as data:
        _validate_schema(data)
        window, step, length = (int(data[key]) for key in ("window", "step", "region_length"))
        positions = np.arange(0, length - window + 1, step)
        labels = data["set_label"].copy()
        lo, hi = data["elig_lo"], data["elig_hi"]
        eligible = (positions[None, :] >= lo[:, None]) & (positions[None, :] < hi[:, None]) & (lo[:, None] >= 0)
        dense = np.zeros(eligible.shape, dtype=bool)
        for event, start, end in zip(data["event_index"], data["hit_start"], data["hit_end"]):
            first = max(0, (int(start) - window) // step + 1)
            last = min(len(positions), (int(end) - 1) // step + 1)
            dense[int(event), first:last] = True
        dense &= eligible
        return dense, eligible, labels
