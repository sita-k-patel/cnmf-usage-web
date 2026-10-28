#!/usr/bin/env python3
"""Compare frozen scientific files with repository pins, using only the standard library."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
data = json.loads((ROOT / "integrity" / "frozen.json").read_text())
errors = []
for rel, expected in data["files"].items():
    path = ROOT / rel
    if not path.is_file():
        errors.append(f"MISSING: {rel}")
        continue
    got = hashlib.sha256(path.read_bytes()).hexdigest()
    if got != expected:
        errors.append(f"SHA MISMATCH: {rel}")
    else:
        print(f"OK: {rel}")
if errors:
    raise SystemExit("\n".join(errors))
print("PASS: Frozen scientific file integrity")
