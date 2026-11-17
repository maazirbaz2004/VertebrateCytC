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

