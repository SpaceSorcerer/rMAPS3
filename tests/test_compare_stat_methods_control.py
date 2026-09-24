"""compare_stat_methods.py checks the --positive-control symbol on every arm, and no symbol by default; the arm name is
never read as a control (review 5, arm semantics)."""
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COLUMNS = ["smallest_p_in_upstreamExonIntron", "smallest_p_in_upstreamIntron", "smallest_p_in_targetExon-5prime",
           "smallest_p_in_targetExon-3prime", "smallest_p_in_downstreamIntron", "smallest_p_in_downstreamExonIntron"]
ARMS = ("QKI_KO_T", "ARBITRARY")


def inputs(tmp_path):
    runs = tmp_path / "runs"
    for arm in ARMS:
        root = runs / arm / "released_b9a9dce"
        root.mkdir(parents=True)
        for direction in ("up", "dn"):
            (root / f"pVal.{direction}.vs.bg.RNAmap.txt").write_text(
                "RBP\t" + "\t".join(COLUMNS) + "\n"
                + "".join(f"{motif}\t" + "\t".join([p] * len(COLUMNS)) + "\n"
                          for motif, p in (("QKI.ACUAAY", "0.001"), ("PTBP1.UCUU", "0.01"), ("RBFOX2.UGCAUG", "0.5"))),
                encoding="utf-8")
    alias = tmp_path / "alias.tsv"
    alias.write_text("table_name\thgnc_symbol\n", encoding="utf-8")
    return runs, alias


def run(tmp_path, name, *extra):
    runs, alias = inputs(tmp_path) if not (tmp_path / "runs").exists() else (tmp_path / "runs", tmp_path / "alias.tsv")
    out = tmp_path / name
    proc = subprocess.run([sys.executable, str(ROOT / "tools" / "compare_stat_methods.py"), "--arms", *ARMS,
                           "--runs-root", str(runs), "--stat-root", str(tmp_path / "stat"), "--alias", str(alias),
                           "--out", str(out), *extra], capture_output=True, text=True, cwd=ROOT)
    assert proc.returncode == 0, proc.stderr[-2000:]
    with open(out / "positive_control.tsv", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t")), (out / "REPORT.md").read_text(encoding="utf-8")


def test_no_positive_control_by_default_even_for_a_qki_named_arm(tmp_path):
    rows, report = run(tmp_path, "default")
    assert rows == [] and "no --positive-control given" in report


def test_the_flag_checks_its_symbol_on_every_arm(tmp_path):
    rows, report = run(tmp_path, "ptbp1", "--positive-control", "PTBP1")
    assert "## Positive control PTBP1" in report and "QKI positive control" not in report
    assert {r["arm"] for r in rows} == set(ARMS) and {r["symbol"] for r in rows} == {"PTBP1"}
    released = [r for r in rows if r["method"] == "released_fisher"]
    assert len(released) == 2 * len(ARMS) and {r["symbol_rank"] for r in released} == {"2"}
    assert {(r["pooled_region"], r["direction"]) for r in released} == {("Upstream Intron", "included"),
                                                                        ("Downstream Intron", "skipped")}
