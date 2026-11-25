#!/usr/bin/env python3
"""Command Line Interface for DNA Sequence Matcher."""

import argparse
import sys
from pathlib import Path

from dna_matcher.validator import validate_sequence, parse_fasta, InvalidSequenceError
from dna_matcher.processor import calculate_metrics
from dna_matcher.matcher import exact_search, fuzzy_search, double_stranded_search
from dna_matcher.reporter import format_terminal_report, export_report_json, export_report_csv


def main():
    parser = argparse.ArgumentParser(
        description="DNA Sequence Matcher: Slices, validates, and compares biological sequences."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sequence", type=str, help="Raw DNA sequence string.")
    group.add_argument("--fasta", type=str, help="Path to input FASTA (.fasta/.fa) file.")

    parser.add_argument("--pattern", type=str, required=True, help="Target nucleotide pattern/motif to search.")
    parser.add_argument("--mismatches", type=int, default=0, help="Maximum allowed mismatches for fuzzy search (default: 0).")
    parser.add_argument("--reverse-complement", action="store_true", help="Search both forward and reverse-complement strands.")
    parser.add_argument("--output-json", type=str, help="Path to export report as JSON.")
    parser.add_argument("--output-csv", type=str, help="Path to export matches as CSV.")

    args = parser.parse_args()

    try:
        if args.fasta:
            fasta_path = Path(args.fasta)
            records = parse_fasta(fasta_path)
        else:
            cleaned_seq = validate_sequence(args.sequence)
            records = {"CLI_Input_Sequence": cleaned_seq}

        for name, seq in records.items():
            metrics = calculate_metrics(seq)

            if args.reverse_complement:
                matches = double_stranded_search(seq, args.pattern, max_mismatches=args.mismatches)
            elif args.mismatches > 0:
                matches = fuzzy_search(seq, args.pattern, max_mismatches=args.mismatches)
            else:
                matches = exact_search(seq, args.pattern)

            report_data = {
                "record_name": name,
                "metrics": metrics,
                "pattern": args.pattern,
                "max_mismatches": args.mismatches,
                "reverse_complement": args.reverse_complement,
                "total_matches": len(matches),
                "matches": matches,
            }

            terminal_output = format_terminal_report(
                sequence_name=name,
                metrics=metrics,
                matches=matches,
                pattern=args.pattern,
                max_mismatches=args.mismatches,
            )
            print(terminal_output)

            if args.output_json:
                json_path = Path(args.output_json) if len(records) == 1 else Path(args.output_json).with_stem(f"{args.output_json}_{name}")
                export_report_json(report_data, json_path)
                print(f"Exported JSON report to: {json_path}")

            if args.output_csv:
                csv_path = Path(args.output_csv) if len(records) == 1 else Path(args.output_csv).with_stem(f"{args.output_csv}_{name}")
                export_report_csv(matches, csv_path)
                print(f"Exported CSV matches to: {csv_path}")

    except (InvalidSequenceError, FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
