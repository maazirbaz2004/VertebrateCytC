# PHASE 0: Environment Setup and Repository Scaffolding

**Read `CLAUDE.md` first.**

## 0. Objective

Create a reproducible workspace: folder structure, conda environment with every bioinformatics tool needed later (Clustal Omega, MUSCLE v5, IQ-TREE, FastTree, trimAl), a central configuration file, shared utilities (logging, config loading), a Makefile, a progress tracker, and a smoke test that proves the toolchain works.

**Serves:** all research questions (infrastructure).

## 1. Prerequisites

- A machine with internet access (needed to download packages and, in Phase 1, sequences).
- `conda` or `mamba` installed (Miniforge is recommended: https://github.com/conda-forge/miniforge).
- **Windows users:** run everything inside **WSL2 (Ubuntu)**. Clustal Omega, MUSCLE, and IQ-TREE are Linux/macOS-native tools, and native Windows is painful.
- **Apple Silicon (M1/M2/M3) users:** if a bioconda package has no `osx-arm64` build, create the environment with `CONDA_SUBDIR=osx-64 mamba env create -f environment.yml` (runs under Rosetta).

## 2. Inputs

None (this phase creates the foundation).

## 3. Tasks

### Task 0.1: Create the directory tree

From the intended project parent directory:

```bash
mkdir -p cytc-phylo && cd cytc-phylo
mkdir -p config phases \
         data/raw data/interim data/processed \
         results/alignments results/conservation results/similarity \
         results/trees results/interpretation results/figures results/tables \
         src/cytc scripts tests logs report
touch src/cytc/__init__.py
```

Copy `CLAUDE.md`, `AGENTS.md` and the `phases/` folder into the repo root if they are not already there.

### Task 0.2: Write `environment.yml`

```yaml
name: cytc
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - biopython>=1.83
  - pandas
  - numpy
  - scipy
  - matplotlib
  - seaborn
  - pyyaml
  - requests
  - dendropy
  - jupyterlab
  - pytest
  - tqdm
  # command-line bioinformatics tools
  - clustalo          # Clustal Omega
  - muscle>=5         # MUSCLE v5 (NOTE: CLI differs from v3, see Phase 3)
  - iqtree            # IQ-TREE (binary may be `iqtree2` or `iqtree3`; detect at runtime)
  - fasttree
  - trimal
  - pip
  - pip:
      - logomaker     # sequence logos
```

Create and activate:

```bash
mamba env create -f environment.yml     # or: conda env create -f environment.yml
conda activate cytc
```

**Pitfalls**
- Channel order matters: `conda-forge` before `bioconda`. Do not add `defaults`.
- If `muscle>=5` resolves to v3.8 (older builds), pin explicitly: `muscle=5.1`. Verify with `muscle -version`.
- IQ-TREE's executable name varies by version (`iqtree`, `iqtree2`, `iqtree3`). The smoke test below auto-detects it. Store the detected name in `config.yaml` (key `tools.iqtree`).
- FastTree's executable may be `FastTree` or `fasttree`. Detect likewise.

### Task 0.3: Write `config/config.yaml`

This file centralizes every tunable parameter. Later phases **read** it and must not duplicate its values. Create it exactly as below, then adjust tool names after Task 0.6.

```yaml
project:
  name: cytc-phylo
  random_seed: 42
  n_threads: 4            # adjust to the machine

paths:
  raw: data/raw
  interim: data/interim
  processed: data/processed
  alignments: results/alignments
  conservation: results/conservation
  similarity: results/similarity
  trees: results/trees
  interpretation: results/interpretation
  figures: results/figures
  tables: results/tables
  logs: logs
  species_table: config/species.csv

tools:
  clustalo: clustalo
  muscle: muscle
  iqtree: iqtree2          # set to the detected binary in Task 0.6
  fasttree: FastTree       # set to the detected binary in Task 0.6
  trimal: trimal

data_collection:           # Phase 1
  uniprot_rest: "https://rest.uniprot.org/uniprotkb/search"
  ncbi_email: "CHANGE_ME@example.com"   # NCBI requires an email for Entrez
  prefer_reviewed: true
  vertebrate_len_min: 100
  vertebrate_len_max: 112
  outgroup_len_min: 100
  outgroup_len_max: 115
  heme_motif_regex: "C..CH"    # CxxCH heme-binding motif
  request_pause_sec: 0.5

sequence_prep:             # Phase 2
  allowed_alphabet: "ACDEFGHIKLMNPQRSTVWY"
  max_ambiguous_fraction: 0.0   # no X/B/Z/U/J allowed in the final dataset

alignment:                 # Phase 3
  primary: clustalo        # clustalo | muscle
  clustalo_extra: []
  muscle_extra: []
  trimal_mode: "-gt 0.8"   # only for the sensitivity analysis, not the primary result

conservation:              # Phase 4
  vertebrate_only_for_conservation: true
  invariant_threshold: 1.0
  high_threshold: 0.90
  moderate_threshold: 0.70
  smoothing_window: 5
  min_region_length: 3
  max_gap_fraction_per_column: 0.5

similarity:                # Phase 5
  gap_open: -10.0
  gap_extend: -0.5
  substitution_matrix: BLOSUM62

phylogenetics:             # Phase 6
  outgroups: [Saccharomyces_cerevisiae, Drosophila_melanogaster]
  primary_root_outgroup: Saccharomyces_cerevisiae
  bootstrap_replicates_nj: 1000
  ufboot_replicates: 1000
  alrt_replicates: 1000
  model_selection: MFP       # ModelFinder Plus
  fasttree_model_flags: "-lg -gamma"
  support_threshold: 70      # report clades with >= this UFBoot support as "supported"

reporting:                   # Phase 8
  dpi: 300
  figure_formats: [png, svg]
  font_size: 10
```

### Task 0.4: Write `src/cytc/utils.py`

```python
"""Shared helpers: repo root, config loading, logging, path resolution."""
from __future__ import annotations

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load config/config.yaml (or a given path) as a dict."""
    cfg_path = Path(path) if path else REPO_ROOT / "config" / "config.yaml"
    with open(cfg_path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def resolve(cfg: dict[str, Any], key: str) -> Path:
    """Return an absolute Path for cfg['paths'][key], creating the directory
    if it looks like a directory (no file suffix)."""
    p = REPO_ROOT / cfg["paths"][key]
    if not p.suffix:
        p.mkdir(parents=True, exist_ok=True)
    return p


def get_logger(script_name: str, cfg: dict[str, Any] | None = None) -> logging.Logger:
    """Create a logger writing to stdout and logs/<script>_<timestamp>.log."""
    cfg = cfg or load_config()
    log_dir = REPO_ROOT / cfg["paths"]["logs"]
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{script_name}_{stamp}.log"

    logger = logging.getLogger(script_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    for handler in (logging.StreamHandler(sys.stdout), logging.FileHandler(log_file)):
        handler.setFormatter(fmt)
        logger.addHandler(handler)
    logger.info("Log file: %s", log_file)
    return logger
```

### Task 0.5: Write the smoke-test script `scripts/00_check_env.py`

```python
"""Phase 0 smoke test: verify Python packages and CLI tools, detect binary names."""
import importlib
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cytc.utils import get_logger, load_config  # noqa: E402

PY_PACKAGES = ["Bio", "pandas", "numpy", "scipy", "matplotlib", "seaborn",
               "yaml", "requests", "dendropy", "logomaker", "tqdm"]

CANDIDATES = {
    "clustalo": ["clustalo"],
    "muscle": ["muscle"],
    "iqtree": ["iqtree3", "iqtree2", "iqtree"],
    "fasttree": ["FastTree", "fasttree"],
    "trimal": ["trimal"],
}
VERSION_FLAGS = {"clustalo": ["--version"], "muscle": ["-version"],
                 "iqtree": ["--version"], "fasttree": ["-expert"],  # FastTree prints version on bad usage
                 "trimal": ["--version"]}


def main() -> None:
    cfg = load_config()
    log = get_logger("00_check_env", cfg)
    failures = []

    for pkg in PY_PACKAGES:
        try:
            mod = importlib.import_module(pkg)
            log.info("OK   python:%s %s", pkg, getattr(mod, "__version__", ""))
        except ImportError:
            log.error("MISSING python package: %s", pkg)
            failures.append(pkg)

    detected = {}
    for tool, names in CANDIDATES.items():
        found = next((n for n in names if shutil.which(n)), None)
        if not found:
            log.error("MISSING tool: %s (tried %s)", tool, names)
            failures.append(tool)
            continue
        detected[tool] = found
        try:
            out = subprocess.run([found] + VERSION_FLAGS[tool], capture_output=True,
                                 text=True, timeout=20)
            first = (out.stdout or out.stderr).strip().splitlines()[:1]
            log.info("OK   tool:%s -> %s | %s", tool, found, first)
        except Exception as exc:  # noqa: BLE001
            log.warning("Tool %s found at %s but version check failed: %s", tool, found, exc)

    log.info("DETECTED BINARY NAMES (copy into config.yaml tools:): %s", detected)
    if failures:
        log.error("SUMMARY: FAILED, missing: %s", failures)
        sys.exit(1)
    log.info("SUMMARY: environment OK; detected=%s", detected)


if __name__ == "__main__":
    main()
```

Run it: `python scripts/00_check_env.py`. Copy the detected binary names into `config.yaml` -> `tools:`.

### Task 0.6: Tiny functional test of the aligners and tree builder

Prove the tools actually run (not just that they exist), using a toy FASTA:

```bash
cat > /tmp/toy.fasta << 'EOF'
>a
MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKG
>b
MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKG
>c
MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNRG
>d
GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKG
EOF
clustalo -i /tmp/toy.fasta -o /tmp/toy_clustalo.fasta --force -v
muscle -align /tmp/toy.fasta -output /tmp/toy_muscle.fasta      # MUSCLE v5 syntax
iqtree2 -s /tmp/toy_clustalo.fasta -m LG --prefix /tmp/toy_iq -redo   # use the detected binary name
FastTree -lg /tmp/toy_clustalo.fasta > /tmp/toy_ft.tree
```

All four commands must succeed and produce output files. (This toy FASTA is only a tool test, **not** project data.)

### Task 0.7: Makefile

```makefile
PY=python
.PHONY: env check phase1 phase2 phase3 phase4 phase5 phase6 phase7 phase8 all test clean-results

check:
	$(PY) scripts/00_check_env.py

phase1:
	$(PY) scripts/10_select_and_fetch.py
phase2:
	$(PY) scripts/20_prepare_sequences.py
phase3:
	$(PY) scripts/30_align.py
phase4:
	$(PY) scripts/40_conservation.py
phase5:
	$(PY) scripts/50_similarity.py
phase6:
	$(PY) scripts/60_phylogenetics.py
phase7:
	$(PY) scripts/70_interpretation.py
phase8:
	$(PY) scripts/80_make_report_assets.py

all: phase1 phase2 phase3 phase4 phase5 phase6 phase7 phase8

test:
	pytest -q

clean-results:
	rm -rf results/alignments/* results/conservation/* results/similarity/* \
	       results/trees/* results/interpretation/* results/figures/* results/tables/*
```

(Later phases may add helper scripts. The Makefile targets above are the single entry points. Keep them in sync.)

### Task 0.8: Other scaffolding files

`PROGRESS.md`:

```markdown
# PROGRESS

| Phase | Status | Date | Notes |
|-------|--------|------|-------|
| 0 Setup | [ ] | | |
| 1 Data collection | [ ] | | |
| 2 Sequence preparation | [ ] | | |
| 3 Alignment | [ ] | | |
| 4 Conservation | [ ] | | |
| 5 Similarity | [ ] | | |
| 6 Phylogenetics | [ ] | | |
| 7 Interpretation | [ ] | | |
| 8 Reporting | [ ] | | |

## Deviations (anything that differs from the phase files, with justification)

## Open issues
```

`.gitignore`:

```
logs/
__pycache__/
*.pyc
.ipynb_checkpoints/
.DS_Store
*.log
```

`tests/test_phase0.py`:

```python
from pathlib import Path
from cytc.utils import REPO_ROOT, load_config


def test_config_loads():
    cfg = load_config()
    assert cfg["project"]["random_seed"] == 42


def test_dirs_exist():
    cfg = load_config()
    for key in ["raw", "interim", "processed", "alignments", "conservation",
                "similarity", "trees", "figures", "tables"]:
        assert (REPO_ROOT / cfg["paths"][key]).exists(), key
```

Run tests with `PYTHONPATH=src pytest -q` (or add a `conftest.py` that inserts `src` into `sys.path`, or `pip install -e .` with a minimal `pyproject.toml`; any of these is fine, pick one and document it).

## 4. Expected outputs

| File | Description |
|------|-------------|
| Full directory tree | as in `CLAUDE.md` section 3 |
| `environment.yml`, `config/config.yaml`, `Makefile`, `PROGRESS.md`, `.gitignore` | scaffolding |
| `src/cytc/utils.py`, `scripts/00_check_env.py` | helpers |
| `logs/00_check_env_*.log` | proof the toolchain works |
| `tests/test_phase0.py` | passing tests |

## 5. Validation checks

- [ ] `conda activate cytc && python scripts/00_check_env.py` exits 0.
- [ ] Toy-run of Clustal Omega, MUSCLE, IQ-TREE and FastTree all produced output files.
- [ ] `config.yaml` `tools:` section contains the *actually detected* binary names.
- [ ] `ncbi_email` in config has been set to a real address (ask the user if unknown; Entrez rejects the placeholder).
- [ ] `pytest -q` passes.
- [ ] `PROGRESS.md` Phase 0 row ticked with the tool versions recorded in Notes.

## 6. Common pitfalls

- **MUSCLE v3 vs v5 syntax mismatch.** v5: `muscle -align in.fa -output out.afa`. v3: `muscle -in in.fa -out out.afa`. This project standardizes on **v5**.
- **Python can't import `cytc`.** Use `sys.path.insert` in scripts as shown, or install the package in editable mode.
- **IQ-TREE writes files next to the prefix.** Always pass `--prefix` pointing into `results/trees/...`.
- **WSL path issues.** Keep the repo inside the Linux filesystem (`~/cytc-phylo`), not under `/mnt/c/...`. It is faster and avoids permission problems.

## 7. Done when

- [ ] Every item in section 5 is checked.
- [ ] Phase 0 committed (`git commit -m "Phase 0: environment and scaffolding"`).
- [ ] `PROGRESS.md` updated.
