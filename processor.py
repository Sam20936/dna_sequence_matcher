"""String Slicing & Sequence Operations Module.

Provides sliding window generation, reverse complement computation, and core
sequence statistics (length, GC-content, base frequency).
"""

from typing import Dict, Iterator, Tuple

# Complement mapping for DNA bases (including IUPAC 'N')
COMPLEMENT_MAP = str.maketrans("ATCGN", "TAGCN")


def sliding_window(sequence: str, k: int, step: int = 1) -> Iterator[Tuple[int, str]]:
    """Generate overlapping or non-overlapping k-mers (windows) from a sequence.

    Args:
        sequence: Validated DNA sequence string.
        k: Window size (k-mer length).
        step: Step size (stride) between windows.

    Yields:
        Tuples of (0-based start index, window_substring).

    Time Complexity:
        O(N * k) where N is sequence length and k is window size.
    Space Complexity:
        O(k) per window generated.
    """
    if k <= 0:
        raise ValueError("Window size k must be greater than 0.")
    if step <= 0:
        raise ValueError("Step size must be greater than 0.")
    if len(sequence) < k:
        return

    for i in range(0, len(sequence) - k + 1, step):
        yield i, sequence[i : i + k]


def reverse_complement(sequence: str) -> str:
    """Compute the biological reverse complement of a DNA sequence.

    Replaces A<->T and C<->G, and reverses the string direction (5'->3' to 3'->5').

    Args:
        sequence: Validated DNA sequence string.

    Returns:
        Reverse complement DNA string.

    Time Complexity:
        O(N) where N is sequence length.
    Space Complexity:
        O(N) for the resulting string.
    """
    return sequence.translate(COMPLEMENT_MAP)[::-1]


def calculate_metrics(sequence: str) -> Dict[str, float | int | Dict[str, int]]:
    """Calculate key sequence metrics including length, GC-content, and base frequencies.

    Args:
        sequence: Validated DNA sequence string.

    Returns:
        Dictionary containing length, GC-content percentage, and counts per base.

    Time Complexity:
        O(N) where N is sequence length.
    Space Complexity:
        O(1) auxiliary space (fixed alphabet size).
    """
    length = len(sequence)
    if length == 0:
        return {"length": 0, "gc_content": 0.0, "base_counts": {}}

    base_counts = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "C": sequence.count("C"),
        "G": sequence.count("G"),
        "N": sequence.count("N"),
    }

    gc_count = base_counts["G"] + base_counts["C"]
    gc_content = (gc_count / length) * 100.0

    return {
        "length": length,
        "gc_content": round(gc_content, 2),
        "base_counts": base_counts,
    }
