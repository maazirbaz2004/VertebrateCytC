"""Phase 8 pytest suite: validate final report, asset manifest, checklist, and submission package."""
from __future__ import annotations

import pandas as pd

from cytc.utils import REPO_ROOT


def test_report_exists_and_has_sections():
    txt = (REPO_ROOT / "report/final_report.md").read_text(encoding="utf-8")
    for h in [
        "Abstract",
        "Introduction",
        "Methodology",
        "Results",
        "Discussion",
        "Limitations",
        "Conclusion",
        "References",
    ]:
        assert h in txt, h


def test_all_figures_in_manifest_exist():
    man = pd.read_csv(REPO_ROOT / "report/asset_manifest.csv")
    assert len(man) >= 24
    for p in man["report_path"]:
        assert (REPO_ROOT / p).exists(), p


def test_deliverables_checklist_complete():
    txt = (REPO_ROOT / "report/deliverables_checklist.md").read_text(encoding="utf-8")
    assert txt.count("[x]") >= 8


def test_no_placeholder_text_left():
    txt = (REPO_ROOT / "report/final_report.md").read_text(encoding="utf-8")
    for bad in ["TODO", "TBD", "XXX", "CHANGE_ME", "lorem"]:
        assert bad not in txt, bad


def test_audits_and_reproducibility_files_exist():
    for rel in [
        "report/citation_audit.md",
        "report/software_versions.md",
        "results/tables/reproducibility_check.csv",
        "environment.lock.yml",
        "dist/cytc_project_submission.zip",
    ]:
        p = REPO_ROOT / rel
        assert p.exists() and p.stat().st_size > 0, rel

