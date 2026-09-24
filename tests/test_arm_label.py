"""No arm-name semantics (review 4 must-fix 2): the figure title comes from --arm-label, the gate record from an
explicit path, the positive control from --positive-control; the arm name is a label and nothing is read from it."""
import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_region_lollipops_v4 as lol  # noqa: E402
import rmaps3_skill_run as skill  # noqa: E402
from test_skill_run import fork_data_engine_run, parse  # noqa: E402

ARM = "QKI_KO_B"                               # a name that used to select the "concordant" title
LABEL = "QKI knockout (rMATS-strict)"
LABEL_TABLE = ROOT / "data" / "arm_labels_dissertation.tsv"


def rule_a_arm(tmp_path):
    out, engine, counts_root, roots, motifs, alias, gtf = fork_data_engine_run(tmp_path, arm=ARM)
    return out, engine, counts_root, roots, motifs, alias, gtf


def quick_figure_args(out, tmp_path, gtf, *extra):
    return parse("--mode", "quick", "--arm", ARM, "--out", str(out), "--up", "u", "--dn", "d", "--bg", "b",
                 "--gate-counts", "c.json", "--gate-rule", "A", "--engine-root", str(ROOT), "--gtf", str(gtf),
                 "--spliceosome-list", str(tmp_path / "splice.txt"),
                 "--broad-binders-list", str(tmp_path / "broad.txt"), *extra)


def svg_texts(directory):
    return {p.name: p.read_text(encoding="utf-8") for p in Path(directory).rglob("*.svg")}


def test_the_builder_has_no_label_table_keyed_by_arm_name():
    assert not hasattr(lol, "LABELS")
    source = (ROOT / "tools" / "build_region_lollipops_v4.py").read_text(encoding="utf-8")
    assert "startswith('QKI_KO')" not in source and "QKI_MAP_NAME" not in source
    assert "arm[:-len(" not in source, "no storage path may be derived from an arm suffix"


def test_dissertation_labels_are_a_data_file_with_the_seven_arms():
    with open(LABEL_TABLE, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["arm"] for r in rows] == ["MIAT_KD_A", "MIAT_KD_B", "MIAT_KD_Beffect", "MIAT_OE_A", "QKI_KO_A",
                                        "QKI_KO_B", "QKI_KO_Beffect"]
    assert all(r["arm_label"].strip() for r in rows)


def test_quick_figures_of_a_rule_a_arm_named_qki_ko_b_carry_the_given_title(tmp_path):
    out, engine, counts_root, roots, motifs, alias, gtf = rule_a_arm(tmp_path)
    args = quick_figure_args(out, tmp_path, gtf, "--arm-label", LABEL)
    scores = {m: skill.motif_scores(counts_root, ARM, m) for m in motifs}
    skill.build_quick_figures(args, out, engine, roots, motifs, scores, alias,
                              ROOT / "data" / "knownMotifs.human.mouse.txt", ROOT / "data" / "ESRP.like.motif.txt",
                              lambda m: None)
    svgs = svg_texts(out / "figures")
    assert len(svgs) == 4
    for name, text in svgs.items():
        assert LABEL in text, name
        assert "concordant" not in text.lower(), name


def test_builder_main_draws_the_given_title_and_refuses_an_arm_without_one(tmp_path):
    out, engine, counts_root, roots, motifs, alias, gtf = rule_a_arm(tmp_path)
    scores_only = tmp_path / "motif_scores"
    scores = {m: skill.motif_scores(counts_root, ARM, m) for m in motifs}
    skill.write_score_table(scores_only / ARM / "per_motif_regions.tsv", ARM, roots, motifs, scores, alias)
    (scores_only / ARM / "command.log").write_text("scores\nexit=0\n", encoding="utf-8")
    base = ["--arms", ARM, "--released-root", str(engine.parent), "--calibrated-root", str(scores_only),
            "--counts-json", f"{ARM}={out / 'event_counts.json'}", "--gate-rule", f"{ARM}=A",
            "--gtf", str(gtf), "--spliceosome-list", str(tmp_path / "splice.txt"),
            "--broad-binders-list", str(tmp_path / "broad.txt"), "--released-commit", skill.git_revision(ROOT)[:7]]
    with pytest.raises(ValueError, match="no figure title for arm QKI_KO_B"):
        lol.main(base + ["--out-root", str(tmp_path / "refused")])
    assert not list((tmp_path / "refused").rglob("*.svg")) if (tmp_path / "refused").exists() else True
    figures = tmp_path / "figs"
    lol.main(base + ["--out-root", str(figures), "--arm-label", f"{ARM}={LABEL}"])
    svgs = svg_texts(figures)
    assert len(svgs) == 4
    for name, text in svgs.items():
        assert LABEL in text and "concordant" not in text.lower(), name
    index = (figures / "index.html").read_text(encoding="utf-8")
    assert LABEL in index and "concordant" not in index.lower()
    assert not (figures / "positive_control_audit.tsv").exists(), "no control is implied by the arm name"


def test_builder_positive_control_comes_from_the_flag_not_the_arm_name(tmp_path):
    out, engine, counts_root, roots, motifs, alias, gtf = fork_data_engine_run(tmp_path, arm="ARBITRARY")
    scores_only = tmp_path / "motif_scores"
    scores = {m: skill.motif_scores(counts_root, "ARBITRARY", m) for m in motifs}
    skill.write_score_table(scores_only / "ARBITRARY" / "per_motif_regions.tsv", "ARBITRARY", roots, motifs,
                            scores, alias)
    (scores_only / "ARBITRARY" / "command.log").write_text("scores\nexit=0\n", encoding="utf-8")
    figures = tmp_path / "figs"
    lol.main(["--arms", "ARBITRARY", "--out-root", str(figures), "--released-root", str(engine.parent),
              "--calibrated-root", str(scores_only), "--counts-json", f"ARBITRARY={out / 'event_counts.json'}",
              "--gate-rule", "ARBITRARY=A", "--arm-label", "ARBITRARY=Arbitrary arm", "--gtf", str(gtf),
              "--spliceosome-list", str(tmp_path / "splice.txt"), "--broad-binders-list", str(tmp_path / "broad.txt"),
              "--released-commit", skill.git_revision(ROOT)[:7], "--positive-control", "QKI"])
    with open(figures / "positive_control_audit.tsv", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 8 and {r["first"] for r in rows} == {"QKI"} and {r["pass"] for r in rows} == {"True"}


def test_gate_counts_take_an_explicit_record_path_only(tmp_path):
    with pytest.raises(ValueError, match="--counts-json QKI_KO_B=<path>"):
        lol.gate_counts("QKI_KO_B", None)
    with pytest.raises(SystemExit):
        lol.main(["--arms", "X", "--out-root", "o", "--released-root", "r", "--calibrated-root", "c", "--gtf", "g",
                  "--spliceosome-list", "s", "--broad-binders-list", "b", "--event-sets-root", "tree"])


def test_wrapper_refuses_figures_without_an_arm_label_and_full_mode_forwards_it():
    presplit = ("--up", "u", "--dn", "d", "--bg", "b", "--gate-counts", "c.json", "--gate-rule", "A")
    with pytest.raises(ValueError, match="--arm-label"):
        skill.validate(parse("--mode", "quick", "--arm", ARM, "--out", "o", *presplit))
    with pytest.raises(ValueError, match="--arm-label"):
        skill.validate(parse("--mode", "full", "--arm", ARM, "--out", "o", *presplit))
    skill.validate(parse("--mode", "quick", "--arm", ARM, "--out", "o", *presplit, "--no-figures"))
    args = parse("--mode", "full", "--arm", ARM, "--out", "o", *presplit, "--arm-label", LABEL,
                 "--engine-root", str(ROOT), "--gtf", "g", "--spliceosome-list", "s", "--broad-binders-list", "b",
                 "--positive-control", "QKI")
    command = skill.figure_command(args, Path("o"), Path("o/engine") / ARM, Path("o/summary"))
    assert command[command.index("--arm-label") + 1] == f"{ARM}={LABEL}"
    assert command[command.index("--positive-control") + 1] == "QKI"
