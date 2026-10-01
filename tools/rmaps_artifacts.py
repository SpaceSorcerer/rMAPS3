"""Registry of the files a run's writers produce, so the completeness inventory is checked against the writers.

Every output-writing function of tools/rmaps3_skill_run.py and tools/build_region_lollipops_v4.py records each path it
writes here; a subprocess step is recorded by a before/after listing of the run directory (`capture`). The run manifest
hashes every recorded path, and `tests/test_artifact_inventory.py` asserts that the recorded set equals the wrapper's
REQUIRED_ARTIFACTS for each mode.
"""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path


class ArtifactRegistry:
    """Paths written under one root, each with the writer that wrote it (posix paths relative to the root)."""

    def __init__(self, root):
        self.root = Path(root).resolve()
        self.entries: dict[str, str] = {}

    def rel(self, path) -> str:
        path = Path(path)
        path = path if path.is_absolute() else self.root / path
        return path.resolve().relative_to(self.root).as_posix()

    def add(self, path, writer: str):
        """Record a written file; returns the path so a writer can register inline."""
        self.entries[self.rel(path)] = writer
        return path

    def discard(self, path) -> None:
        self.entries.pop(self.rel(path), None)

    def files_on_disk(self) -> set:
        return {p.relative_to(self.root).as_posix() for p in self.root.rglob("*") if p.is_file()}

    @contextmanager
    def capture(self, writer: str):
        """Record every file a step creates or rewrites under the root, and forget files it deletes."""
        before = {rel: (self.root / rel).stat().st_mtime_ns for rel in self.files_on_disk()}
        try:
            yield self
        finally:
            after = self.files_on_disk()
            for rel in sorted(after):
                path = self.root / rel
                if rel not in before or path.stat().st_mtime_ns != before[rel]:
                    self.entries.setdefault(rel, writer)
            for rel in set(before) - after:
                self.entries.pop(rel, None)

    def paths(self) -> list:
        return sorted(self.entries)

    def by_writer(self) -> dict:
        out: dict[str, list] = {}
        for rel, writer in sorted(self.entries.items()):
            out.setdefault(writer, []).append(rel)
        return out
