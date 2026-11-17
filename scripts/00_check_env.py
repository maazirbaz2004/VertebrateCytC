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

