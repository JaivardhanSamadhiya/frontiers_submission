#!/usr/bin/env python3
"""Stage the pinned public source repositories for StaphStudy validation."""
from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SOURCES = {
    ROOT / "external/upstream_vhip_tool": (
        "https://github.com/DuhaimeLab/VirusHostInteractionPredictor.git",
        "b62a7b6ba08a056b6e04b25ef4dd3a355ebd179d",
    ),
    ROOT / "external/upstream_vhip": (
        "https://github.com/DuhaimeLab/VHIP_analyses_Bastien_et_al_2023.git",
        "5c87aaf7164f78c11f598e0bb56681932ce8bbad",
    ),
}


def git(*args: str, cwd: Path | None = None) -> str:
    result = subprocess.run(
        ["git", *args], cwd=cwd, check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    return result.stdout.strip()


def main() -> None:
    for path, (url, commit) in SOURCES.items():
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            git("clone", url, str(path))
        observed = git(
            "-c", f"safe.directory={path.as_posix()}", "rev-parse", "HEAD",
            cwd=path,
        )
        if observed != commit:
            raise RuntimeError(
                f"{path} is at {observed}; expected pinned commit {commit}. "
                "Move or remove the directory before restaging."
            )
    print("PASS: pinned StaphStudy validation sources are staged")


if __name__ == "__main__":
    main()
