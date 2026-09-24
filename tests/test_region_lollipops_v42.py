"""Figure builder v4.2: parameterised footer/legend text and a single-layer main() run."""
import argparse
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_region_lollipops_v4 as lol  # noqa: E402
import rmaps3_skill_run as skill  # noqa: E402
from test_skill_run import fork_data_engine_run  # noqa: E402

COUNTS = {"n_up": 20, "n_dn": 18, "n_bg": 200, "n_expr_unknown_in_fg": 0, "n_expr_unknown_in_bg": 2}


def text_args(**overrides):
    base = dict(gate_text=None, direction_text=None, tail_text=None, tail_text_supplement=None,
                gate_record=None, released_commit="abc1234", released_stat_method="mannwhitney")
    base.update(overrides)
    return argparse.Namespace(**base)


def test_version_is_4_3_3():
    assert lol.FIG_VERSION == "4.3.3"


def test_text_precedence_is_command_line_then_record_then_default(tmp_path):
    record = tmp_path / "gate.json"
    record.write_text(json.dumps({"gate_text": "RECORD GATE", "tail_text": "RECORD TAIL"}), encoding="utf-8")
    texts = lol.resolve_texts(text_args(gate_record=str(record), gate_text="CLI GATE"))
    assert texts["gate_text"] == "CLI GATE" and texts["tail_text"] == "RECORD TAIL"
    assert texts["direction_text"] == lol.DEFAULT_DIRECTION_TEXT
    assert texts["source"] == "command line"
    assert lol.resolve_texts(text_args(gate_record=str(record)))["source"].startswith("gate record")
    assert lol.resolve_texts(text_args())["source"] == "dissertation default"
    record.write_text(json.dumps({"gate_txt": "typo"}), encoding="utf-8")
    with pytest.raises(ValueError, match="Unknown gate-record keys"):
        lol.resolve_texts(text_args(gate_record=str(record)))


def test_resolve_texts_accepts_a_namespace_without_text_flags():
    """The skill's quick-mode figures call draw_figure with the skill's own namespace."""
    texts = lol.resolve_texts(argparse.Namespace(stat_method="mannwhitney"))
    assert texts["source"] == "dissertation default" and texts["commit"] == lol.RELEASED_COMMIT


def test_gate_lines_fill_placeholders_and_cap_at_three_lines():
    texts = lol.resolve_texts(text_args(gate_text="n={n_up}/{n_dn}/{n_bg} rule {rule}\\nsecond {unknown}"))
    texts["rule"] = "B"                                          # the stated rule, never the arm-name suffix
    assert lol.gate_lines(texts, "ARM_A", COUNTS) == ["n=20/18/200 rule B", "second {unknown}"]
    default = lol.gate_lines(dict(lol.resolve_texts(text_args()), rule="Beffect"), "ARM_A", COUNTS)
    assert len(default) == 3 and default[-1] == lol.DEFAULT_BEFFECT_LINE
    with pytest.raises(ValueError, match="at most three lines"):
        lol.gate_lines(lol.resolve_texts(text_args(gate_text="a\\nb\\nc\\nd")), "ARM_A", COUNTS)


def test_tail_line_uses_the_given_commit_and_the_supplement_override():
    texts = lol.resolve_texts(text_args(tail_text="T {commit} {stat_method}",
                                        tail_text_supplement="S seed {seed}"))
    assert lol.tail_line(texts, "released_ranksum_rawP", None) == "T abc1234 mannwhitney"
    assert lol.tail_line(texts, "calibrated_ranksum", {"seed": 149}) == "S seed 149"
    default = lol.tail_line(lol.resolve_texts(text_args()), "released_ranksum_rawP", None)
    assert "abc1234" in default and default.startswith("Released engine commit")


def test_length_clauses_are_neutral_without_a_length_root():
    assert lol.length_clause("ARM_A", None) == ("changed-vs-background exon length: see methods", None,
                                                "not_measured")
    assert lol.length_sensitivity_clause("ARM_A", None) == ("", None, "not_measured")


def test_single_layer_main_draws_the_released_layer_with_overridden_text(tmp_path):
    out, engine, counts_root, roots, motifs, alias, gtf = fork_data_engine_run(tmp_path)
    figures = tmp_path / "figs"
    scores_only = tmp_path / "motif_scores"  # released-layer scores only: no refinement_report.json
    scores = {m: skill.motif_scores(counts_root, "QKI_KO_T", m) for m in motifs}
    skill.write_score_table(scores_only / "QKI_KO_T" / "per_motif_regions.tsv", "QKI_KO_T", roots, motifs,
                            scores, alias)
    (scores_only / "QKI_KO_T" / "command.log").write_text("scores\nexit=0\n", encoding="utf-8")
    lol.main(["--arms", "QKI_KO_T", "--out-root", str(figures), "--released-root", str(engine.parent),
              "--calibrated-root", str(scores_only),
              "--counts-json", f"QKI_KO_T={out / 'event_counts.json'}",
              "--gtf", str(gtf), "--spliceosome-list", str(tmp_path / "splice.txt"),
              "--broad-binders-list", str(tmp_path / "broad.txt"),
              "--released-commit", skill.git_revision(ROOT)[:7],
              "--gate-text", "TEST GATE n={n_up}", "--direction-text", "TEST DIRECTION",
              "--tail-text", "TEST TAIL {commit}"])
    drawn = sorted(p.name for p in (figures / "QKI_KO_T").glob("*.svg"))
    assert drawn == sorted(f"QKI_KO_T_SE_{k}_released_ranksum_rawP{v}.svg"
                           for k in ("byRBP", "byMotif") for v in ("", "_noSpliceosome_noBroad"))
    svg = (figures / "QKI_KO_T" / "QKI_KO_T_SE_byRBP_released_ranksum_rawP.svg").read_text(encoding="utf-8")
    assert "TEST GATE n=20" in svg and "TEST DIRECTION" in svg
    assert "TEST TAIL " + skill.git_revision(ROOT)[:7] in svg
    assert "use for RBP ORDER only" in svg
    sidecar = (figures / "QKI_KO_T" / "QKI_KO_T_figure_provenance_v43.md").read_text(encoding="utf-8")
    assert "TEST GATE n=20" in sidecar and "released-Fisher" not in sidecar and "treatment C" not in sidecar
    manifest = (figures / "figures_manifest_v43.tsv").read_text(encoding="utf-8")
    assert "calibrated_ranksum" in manifest, "the skipped supplement must be recorded, not silently dropped"
    assert (figures / "index.html").is_file()


def test_sidecar_statements_follow_the_text_flags_and_default_to_the_dissertation():
    default = lol.sidecar_statements(dict(lol.resolve_texts(text_args()), rule="B"), "QKI_KO_B", COUNTS,
                                     "persample10_bm50_bgfdr0.5", "not_measured")
    assert default[0].startswith("- Gate persample10_bm50_bgfdr0.5, rule B: included 20")
    assert default[1] == lol.DEFAULT_SIDECAR_MD5 and default[2] == lol.DEFAULT_SIDECAR_LENGTH
    texts = lol.resolve_texts(text_args(gate_text="jc10 gate n={n_up}", tail_text="TAIL {commit}"))
    gate, md5, length = lol.sidecar_statements(texts, "osd258_flight_vs_ground", COUNTS, "g", "not_measured")
    joined = gate + md5 + length
    assert "jc10 gate n=20" in gate and "rule" not in gate and "DESeq2" not in gate
    assert "released-Fisher" not in joined and "treatment C" not in joined
    assert md5 == " Tail statement (as printed): TAIL abc1234" and length == ""
    assert "treatment C" not in lol.sidecar_statements(texts, "x_A", COUNTS, "g", "measured_order_unchanged")[2]
