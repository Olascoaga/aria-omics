"""
ARIA per-experiment file logging.

Every `logging.getLogger("aria.X")` is a child of the root "aria" logger.
attach_experiment_log() adds a FileHandler to that root for the duration
of one experiment, so anything an agent logs lands in
  ~/.aria/logs/exp_<experiment_id>.log
Detach when the experiment ends to prevent handler accumulation.
"""

from __future__ import annotations

import logging
from pathlib import Path


_LOG_DIR = Path.home() / ".aria" / "logs"


def attach_experiment_log(experiment_id: str,
                          level: int = logging.DEBUG,
                          log_dir: Path | None = None) -> logging.FileHandler:
    """
    Attach a FileHandler for this experiment_id to the "aria" root logger.
    Returns the handler so the caller can detach it via detach_experiment_log().
    """
    target_dir = Path(log_dir) if log_dir else _LOG_DIR
    target_dir.mkdir(parents=True, exist_ok=True)
    path = target_dir / f"exp_{experiment_id}.log"

    handler = logging.FileHandler(path, encoding="utf-8")
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s — %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    # Tag the handler so detach can find it again even if the caller drops
    # the reference (defensive).
    handler.set_name(f"aria_exp_{experiment_id}")

    root = logging.getLogger("aria")
    if root.level > level or root.level == 0:
        root.setLevel(level)
    root.addHandler(handler)

    root.info(f"=== Experiment {experiment_id} log started ({path}) ===")
    return handler


def detach_experiment_log(handler: logging.FileHandler | None) -> None:
    """Detach and close a previously-attached experiment log handler."""
    if handler is None:
        return
    root = logging.getLogger("aria")
    try:
        root.info(f"=== Experiment log closed ({handler.baseFilename}) ===")
    except Exception:
        pass
    root.removeHandler(handler)
    try:
        handler.close()
    except Exception:
        pass
