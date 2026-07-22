#!/usr/bin/env python3
"""Extract Figure 2--3 source cells verbatim from the original R notebook.

This is a provenance utility, not a rerun notebook: cells retain their original
stateful assumptions and are therefore accompanied by a source-cell index.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
PROJECT_DIR = PACKAGE_DIR.parents[1]
NOTEBOOK = PROJECT_DIR / "batch6.ipynb"
OUT_DIR = PACKAGE_DIR / "provenance"


def expand(*parts: object) -> list[int]:
    result: list[int] = []
    for part in parts:
        if isinstance(part, range):
            result.extend(part)
        else:
            result.append(int(part))
    return result


# Cells were selected from the manuscript panel-to-code audit.  They are copied
# unchanged, so a reviewer can trace every retained statement back to batch6.
SECTIONS = {
    "00_data_loading_and_metadata.R": expand(
        range(93, 106), range(202, 226), range(428, 444), range(630, 644),
        range(925, 931), range(1003, 1010),
    ),
    "01_figure2_cells.R": expand(
        400, 480, 482, 484, range(787, 800), 933, 935, 936, 947,
        953, 963, 965, 966, 978, 980, 986, 987, 1014, 1017,
        1019, 1020, 1021, 1022,
    ),
    "02_figure3_cells.R": expand(
        180, 181, 182, 183, 184, 185, 186, 400, 403, 404,
        range(415, 424), range(453, 455), 587, 588, 589, 590, 591,
        592, 593, 594, 595, 664, 677, 679, 680, 693, 694, 699,
        700, 701, 702, 703, 704, 705, 706, 707, 708, 820, 822,
        823, 825, 826, 827, 839, 843, 844, 845, 848, 849, 850,
        851, 852, 889, 893, 894, 933, 935, 937, 938, 940, 941,
        944, 953, 954, 955, 956, 957, 961, 963, 965, 966, 976,
        977, 978, 980, 982, 983,
    ),
}


def cell_text(cell: dict) -> str:
    return "".join(cell.get("source", []))


def main() -> None:
    if not NOTEBOOK.exists():
        raise FileNotFoundError(f"Original notebook not found: {NOTEBOOK}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    cells = notebook["cells"]
    rows: list[dict[str, str]] = []

    for filename, requested_ids in SECTIONS.items():
        ids = sorted(set(requested_ids))
        bad = [idx for idx in ids if idx < 0 or idx >= len(cells)]
        if bad:
            raise IndexError(f"Invalid cell indexes for {filename}: {bad}")

        lines = [
            "# AUTO-GENERATED PROVENANCE FILE -- DO NOT EDIT\n",
            "# Source: ../../batch6.ipynb\n",
            "# Each block was copied verbatim; it may depend on interactive\n",
            "# objects created outside this extraction. Use the source index and\n",
            "# README before running individual blocks.\n\n",
        ]
        for idx in ids:
            cell = cells[idx]
            text = cell_text(cell)
            first = next((x.strip() for x in text.splitlines() if x.strip()), "<empty>")
            rows.append(
                {
                    "section": filename,
                    "cell_index": str(idx),
                    "cell_type": cell.get("cell_type", ""),
                    "first_line": first,
                }
            )
            lines.append(f"\n# ===== batch6.ipynb cell {idx} ({cell.get('cell_type', '')}) =====\n")
            if cell.get("cell_type") == "markdown":
                for line in text.splitlines():
                    lines.append(f"# {line}\n")
            else:
                lines.append(text if text.endswith("\n") else text + "\n")

        (OUT_DIR / filename).write_text("".join(lines), encoding="utf-8")

    with (OUT_DIR / "source_cell_index.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["section", "cell_index", "cell_type", "first_line"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Extracted {len(rows)} source-cell records to {OUT_DIR}")


if __name__ == "__main__":
    main()
