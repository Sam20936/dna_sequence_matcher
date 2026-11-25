"""Sequence Validation & FASTA Parsing Module.

Handles validation of IUPAC nucleotide sequences and provides a pure Python
FASTA file parser.
"""

from pathlib import Path
from typing import Dict, Iterator, Tuple

# Valid IUPAC nucleotide characters (including 'N' for unknown/any base)
VALID_NUCLEOTIDES = set("ATCGN")


class InvalidSequenceError(ValueError):
    """Raised when a DNA sequence contains invalid characters."""
    pass


def validate_sequence(sequence: str) -> str:
    """Validate and normalize a raw DNA sequence string.

    Strips whitespace and newlines, converts to uppercase, and checks that all
    characters belong to the valid IUPAC alphabet (A, T, C, G, N).

    Args:
        sequence: Raw input DNA sequence string.

    Returns:
        Normalized uppercase DNA sequence string.

    Raises:
        InvalidSequenceError: If any invalid character is encountered.

    Time Complexity:
        O(N) where N is the length of the sequence.
    Space Complexity:
        O(N) for storing the normalized string.
    """
    if not isinstance(sequence, str):
        raise InvalidSequenceError("Sequence must be a string.")

    cleaned = sequence.strip().upper()
    cleaned = "".join(cleaned.split())

    if not cleaned:
        raise InvalidSequenceError("Sequence cannot be empty.")

    invalid_chars = set(cleaned) - VALID_NUCLEOTIDES
    if invalid_chars:
        raise InvalidSequenceError(
            f"Invalid nucleotide character(s) found: {sorted(invalid_chars)}. "
            f"Allowed characters: {sorted(VALID_NUCLEOTIDES)}"
        )

    return cleaned


def parse_fasta(file_path: Path | str) -> Dict[str, str]:
    """Parse a multi-record FASTA file into a dictionary of header -> sequence.

    A pure Python parser that handles multi-line sequences and validates each record.

    Args:
        file_path: Path to the .fasta or .fa file.

    Returns:
        A dictionary mapping sequence header IDs to their normalized sequences.

    Raises:
        FileNotFoundError: If the FASTA file does not exist.
        InvalidSequenceError: If any record contains invalid nucleotides.

    Time Complexity:
        O(N) where N is the total number of characters in the file.
    Space Complexity:
        O(N) to store sequences in memory.
    """
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    sequences: Dict[str, str] = {}
    current_header = None
    current_seq_lines = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    raw_seq = "".join(current_seq_lines)
                    sequences[current_header] = validate_sequence(raw_seq)
                current_header = line[1:].strip()
                current_seq_lines = []
            else:
                current_seq_lines.append(line)

        if current_header is not None:
            raw_seq = "".join(current_seq_lines)
            sequences[current_header] = validate_sequence(raw_seq)

    if not sequences:
        raise InvalidSequenceError(f"No valid FASTA records found in {path}")

    return sequences
