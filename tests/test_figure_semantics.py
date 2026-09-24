"""Main-layer colour semantics and engine provenance in every rendered text (review findings 2 and 9)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_region_lollipops_v4 as lol  # noqa: E402
import rmaps3_skill_run as skill  # noqa: E402
import synthetic_calibrated_arm as syn  # noqa: E402

SENTENCE = "raw rank-sum p (released engine): ORDER ONLY — grey = p ≥ 0.05, not a significance claim"


def svg(figures, kind, layer):
    return (figures / syn.ARM / f"{syn.ARM}_SE_{kind}_{layer}.svg").read_text(encoding="utf-8")


def test_the_sentence_is_one_constant():
    assert lol.RAW_COLOUR_SENTENCE == SENTENCE


def test_main_layer_key_and_legend_say_order_only_and_never_not_significant(calibrated_arm):
    for kind in ("byRBP", "byMotif"):
        text = svg(calibrated_arm["figures"], kind, "released_ranksum_rawP")
        assert SENTENCE in text, kind                                   # legend line
        assert lol.RAW_KEY_TITLE in text and lol.RAW_KEY_GREY in text   # colour key, two lines
        assert "not significant" not in text
        assert "use for RBP ORDER only" in text
    supplement = svg(calibrated_arm["figures"], "byRBP", "calibrated_ranksum")
    assert "q ≥ 0.05 (not significant)" in supplement and SENTENCE not in supplement


def test_engine_commit_and_stat_method_propagate_into_every_text(calibrated_arm):
    figures = calibrated_arm["figures"]
    main = svg(figures, "byRBP", "released_ranksum_rawP")
    sidecar = (figures / syn.ARM / f"{syn.ARM}_figure_provenance_v43.md").read_text(encoding="utf-8")
    index = (figures / "index.html").read_text(encoding="utf-8")
    import openpyxl
    book = openpyxl.load_workbook(figures / syn.ARM / f"{syn.ARM}_rank_comparison.xlsx")
    methods = " ".join(str(c.value) for row in book["Methods"].iter_rows() for c in row)
    for name, text in (("main svg", main), ("sidecar", sidecar), ("index", index), ("methods", methods)):
        assert syn.COMMIT in text and "mannwhitney" in text, name
        assert lol.RELEASED_COMMIT not in text, name
    for name, text in (("sidecar", sidecar), ("index", index)):
        assert SENTENCE in text, name
        assert "constant size" not in text and "constant dot size" not in text, name
        assert "count_ratio" in text, name


def test_quick_mode_figures_and_index_carry_the_same_semantics(tmp_path):
    from test_skill_run import build_fork_figures, fake_counts
    args, out, engine, result = build_fork_figures(tmp_path)
    commit = skill.git_revision(ROOT)[:7]
    text = (out / "figures" / "QKI_KO_T_SE_byRBP_released_ranksum_rawP.svg").read_text(encoding="utf-8")
    assert SENTENCE in text and commit in text and "not significant" not in text
    html = skill.write_index(args, out, engine, fake_counts(), [], {"verified": True}, [], result).read_text(
        encoding="utf-8")
    assert SENTENCE in html and commit in html and "bins the same raw p" not in html
