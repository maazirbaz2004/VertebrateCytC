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

