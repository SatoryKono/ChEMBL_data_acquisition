"""CLI proxy for validating post-processing pipeline configuration."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from library.config.postprocessing import main


if __name__ == "__main__":
    raise SystemExit(main())
