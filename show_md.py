#!/usr/bin/env python3
"""Render a Markdown file nicely in the terminal using rich."""

import sys
from pathlib import Path

from rich.console import Console
from rich.markdown import Markdown


def main() -> int:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <file.md>", file=sys.stderr)
        return 1

    path = Path(sys.argv[1])
    if not path.exists():
        print(f"Error: file not found: {path}", file=sys.stderr)
        return 1

    console = Console()
    console.print(Markdown(path.read_text()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
