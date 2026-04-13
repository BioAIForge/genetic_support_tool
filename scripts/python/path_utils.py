from __future__ import annotations

import os
from pathlib import Path


def resolve_output_root(project_root: Path, configured_root: str) -> Path:
    override = os.environ.get("GENETIC_TOOL_OUTPUT_ROOT")
    if override:
        return Path(override)
    configured_path = Path(configured_root)
    if configured_path.is_absolute():
        return configured_path
    return project_root / configured_path
