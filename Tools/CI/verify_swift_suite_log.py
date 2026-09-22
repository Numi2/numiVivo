#!/usr/bin/env python3
"""Reject filtered Swift Testing runs that execute no tests or omit their suite."""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys


def validate_log(text: str, suite: str) -> int:
    if not suite.strip():
        raise ValueError("A nonempty Swift Testing suite name is required")
    text = re.sub(r"\x1b\[[0-9;]*m", "", text)
    suite_pattern = r'\bSuite "?' + re.escape(suite) + r'"? passed(?: after|\s|$)'
    if not re.search(suite_pattern, text):
        raise ValueError("The requested Swift Testing suite did not report success")
    counts = re.findall(r"\bTest run with (\d+) tests?(?: in \d+ suites?)? passed\b", text)
    if len(counts) != 1 or int(counts[0]) <= 0:
        raise ValueError("Expected exactly one successful, nonempty Swift Testing run")
    # The dedicated filtered run must not pass by skipping its learning tests.
    # Match the status after the complete quoted or unquoted name. Failure
    # words inside a test name are harmless; quoted skip reasons are not.
    bad_status = (
        r'(?m)^[^\w\n]*(?:Test|Suite)\s+'
        r'(?:"(?:\\.|[^"\\\n])*"|[^"\n]*?)\s+(?:failed|skipped)\b'
    )
    if re.search(bad_status, text, re.IGNORECASE) or re.search(r"(?mi)^\s*error:", text):
        raise ValueError("The learning run contains skipped tests or failures")
    return int(counts[0])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--suite", required=True)
    args = parser.parse_args()
    try:
        count = validate_log(args.log.read_text(encoding="utf-8"), args.suite)
    except (OSError, UnicodeError, ValueError) as exc:
        print("Test execution check: " + str(exc), file=sys.stderr)
        return 1
    print(f"Verified {count} executed software tests in {args.suite}; not biological validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
