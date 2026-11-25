"""Reporting & Visualization Module.

Formats search results for terminal output and exports reports to JSON or CSV.
"""

import json
import csv
from pathlib import Path
from typing import List, Dict, Any


def format_terminal_report(
    sequence_name: str,
    metrics: Dict[str, Any],
    matches: List[Dict[str, Any]],
    pattern: str,
    max_mismatches: int,
) -> str:
    """Format a clean, professional terminal summary report.

    Args:
        sequence_name: Name or ID of the sequence record.
        metrics: Dictionary of sequence metrics.
        matches: List of match dictionaries.
        pattern: Target pattern searched.
        max_mismatches: Mismatch threshold used.

    Returns:
        Formatted multi-line report string.
    """
    lines = []
    lines.append("=" * 65)
    lines.append(f"🧬 DNA SEQUENCE MATCHER REPORT")
    lines.append("=" * 65)
    lines.append(f"Record ID      : {sequence_name}")
    lines.append(f"Sequence Length: {metrics['length']} bp")
    lines.append(f"GC Content     : {metrics['gc_content']}%")
    lines.append(f"Target Pattern : {pattern}")
    lines.append(f"Max Mismatches : {max_mismatches}")
    lines.append(f"Total Matches  : {len(matches)}")
    lines.append("-" * 65)

    if not matches:
        lines.append("No matches found matching the criteria.")
    else:
        lines.append(f"{'#':<4} {'Coord (1-based)':<16} {'Strand':<10} {'Mismatches':<12} {'Match / Visualization'}")
        lines.append("-" * 65)
        for idx, m in enumerate(matches, 1):
            coord_str = f"{m['start_coordinate_1based']}..{m['end_coordinate_1based']}"
            strand = m["strand"]
            mismatches = m["mismatches"]
            vis = m["visualization"]
            lines.append(f"{idx:<4} {coord_str:<16} {strand:<10} {mismatches:<12} {vis}")

    lines.append("=" * 65)
    return "\n".join(lines)


def export_report_json(report_data: Dict[str, Any], output_path: Path | str) -> None:
    """Export report dictionary to a JSON file.

    Args:
        report_data: Complete dictionary containing metrics and matches.
        output_path: Destination file path.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(report_data, f, indent=4)


def export_report_csv(matches: List[Dict[str, Any]], output_path: Path | str) -> None:
    """Export match records to a CSV file.

    Args:
        matches: List of match dictionaries.
        output_path: Destination file path.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    fieldnames = [
        "start_coordinate_1based",
        "end_coordinate_1based",
        "start_index_0based",
        "end_index_0based",
        "strand",
        "mismatches",
        "matched_sequence",
        "visualization",
    ]

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for match in matches:
            writer.writerow(match)
