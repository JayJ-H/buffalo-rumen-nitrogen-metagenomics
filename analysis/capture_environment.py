"""Record software versions without local installation paths."""
from pathlib import Path
import importlib.metadata
import json
import platform

root = Path(__file__).resolve().parents[1]
environment = root / "environment"
environment.mkdir(exist_ok=True)
versions = {"Python": platform.python_version()}
for package in ["numpy", "pandas", "scipy", "pypdf", "Pillow"]:
    versions[package] = importlib.metadata.version(package)
(environment / "Python_versions.json").write_text(json.dumps(versions, indent=2))
