from __future__ import annotations
import os, sys, shutil
from pathlib import Path

def bundled_base_dir() -> Path:
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
    return Path(__file__).resolve().parent

def etminflux_home() -> Path:
    override = os.environ.get("ETMINFLUX_HOME")
    if override:
        return Path(override)
    return Path.home() / "Documents" / "etMINFLUX"

def ensure_user_tree() -> dict[str, Path]:
    home = etminflux_home()
    cfg_dir = home / "config_files"
    ana_dir = home / "analysis_pipelines"
    tr_dir  = home / "transform_pipelines"

    cfg_dir.mkdir(parents=True, exist_ok=True)
    ana_dir.mkdir(parents=True, exist_ok=True)
    tr_dir.mkdir(parents=True, exist_ok=True)

    defaults = bundled_base_dir() / "defaults"

    # config
    src_json = defaults / "etMINFLUX_setup.json"
    dst_json = cfg_dir / "etMINFLUX_setup.json"
    if src_json.exists() and not dst_json.exists():
        shutil.copy2(src_json, dst_json)

    # copy default guidetext if missing
    src_guidetext = defaults / "guidetext.md"
    dst_guidetext = cfg_dir / "guidetext.md"
    if src_guidetext.exists() and not dst_guidetext.exists():
        shutil.copy2(src_guidetext, dst_guidetext)

    # pipelines: populate if empty
    for src, dst in [(defaults / "analysis_pipelines", ana_dir),
                     (defaults / "transform_pipelines", tr_dir)]:
        if src.exists() and not any(dst.iterdir()):
            for item in src.rglob("*"):
                rel = item.relative_to(src)
                target = dst / rel
                if item.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    if not target.exists():
                        shutil.copy2(item, target)

    return {
        "home": home,
        "config_dir": cfg_dir,
        "analysis_dir": ana_dir,
        "transform_dir": tr_dir,
        "setup_json": dst_json,
        "guidetext": dst_guidetext,
    }