#!/usr/bin/env python3
"""Verify the release's recorded bytes. This is not a proof replay."""
from __future__ import annotations
import hashlib
import sys
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    manifest = root / "MANIFEST.sha256"
    checked = 0
    seen: set[str] = set()
    try:
        for line_number, line in enumerate(manifest.read_text().splitlines(), 1):
            if not line.strip():
                continue
            expected, separator, relative = line.partition("  ")
            if not separator or len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
                raise ValueError(f"Malformed manifest line {line_number}")
            if relative in seen:
                raise ValueError(f"Duplicate manifest entry: {relative}")
            seen.add(relative)
            candidate = root / relative
            candidate.resolve().relative_to(root.resolve())
            if Path(relative).is_absolute() or ".." in Path(relative).parts or not candidate.is_file():
                raise ValueError(f"Invalid or missing file: {relative}")
            digest = hashlib.sha256()
            with candidate.open("rb") as stream:
                for block in iter(lambda: stream.read(1 << 20), b""):
                    digest.update(block)
            if digest.hexdigest() != expected:
                raise ValueError(f"Hash mismatch: {relative}")
            checked += 1
    except (OSError, ValueError) as exc:
        print(f"REJECT_MANIFEST: {exc}", file=sys.stderr)
        return 1
    if checked == 0:
        print("REJECT_MANIFEST: empty manifest", file=sys.stderr)
        return 1
    print(f"PASS_MANIFEST: {checked} recorded files; integrity only, not a mathematical proof check.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
