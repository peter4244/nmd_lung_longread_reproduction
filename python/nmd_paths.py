"""Path resolution for the NMD long-read analysis — Python mirror of R/load_config.R.

Same contract: environment variable NMD_<KEY> wins, then config/paths.yml, then
the value is interpreted relative to the repository root.

    from nmd_paths import P
    P["MODEL"] / "results_4ct" / "tx_summary.tsv"
"""
import os
from pathlib import Path


def repo_root(start: Path | None = None) -> Path:
    d = (start or Path.cwd()).resolve()
    for cand in (d, *d.parents):
        if (cand / "config" / "paths.yml").exists():
            return cand
    raise RuntimeError("repo root not found (no config/paths.yml above %s)" % d)


def nmd_paths(root: Path | None = None) -> dict[str, Path]:
    root = root or repo_root(Path(__file__).parent)
    out: dict[str, Path] = {}
    for line in (root / "config" / "paths.yml").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or ":" not in line:
            continue
        key, val = (s.strip() for s in line.split(":", 1))
        if not val:
            continue
        val = os.environ.get(f"NMD_{key}", val)
        p = Path(val).expanduser()
        out[key] = p if p.is_absolute() else (root / p)
    out["ROOT"] = root
    return {k: Path(os.path.normpath(v)) for k, v in out.items()}


P = nmd_paths()
