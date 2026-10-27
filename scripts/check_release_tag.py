#!/usr/bin/env python3
"""Check that a git ref is a version tag equal to the pyproject.toml version.

Usage: check_release_tag.py REF [PYPROJECT]   (REF like refs/tags/v0.1.5)
Exit status 0 if they match, 1 otherwise. Used by .github/workflows/publish.yml.
"""
import re
import sys


def pyproject_version(path="pyproject.toml"):
    with open(path, encoding="utf-8") as f:
        m = re.search(r'^version\s*=\s*"([^"]+)"', f.read(), re.M)
    if not m:
        raise ValueError(f"no version in {path}")
    return m.group(1)


def check(ref, path="pyproject.toml"):
    """Return (ok, message)."""
    m = re.fullmatch(r"refs/tags/v(\d+\.\d+\.\d+)", ref)
    if not m:
        return False, f"{ref!r} is not a refs/tags/vX.Y.Z tag"
    version = pyproject_version(path)
    if m.group(1) != version:
        return False, f"tag version {m.group(1)} != pyproject version {version}"
    return True, f"tag v{version} matches pyproject"


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    ok, msg = check(argv[1], argv[2] if len(argv) > 2 else "pyproject.toml")
    print(("OK: " if ok else "ERROR: ") + msg)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
