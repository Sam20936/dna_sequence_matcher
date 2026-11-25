"""Pattern Matching Engine Module.

Supports exact sub-sequence searches, fuzzy/point-mutation matching (Hamming distance),
and simultaneous double-stranded searches.
"""

from typing import List, Dict, Any
from dna_matcher.validator import validate_sequence
from dna_matcher.processor import reverse_complement, sliding_window


def hamming_distance(seq1: str, seq2: str) -> int:
    """Calculate the Hamming distance between two strings of equal length.

    Args:
        seq1: First nucleotide string.
        seq2: Second nucleotide string.

    Returns:
        Number of mismatched positions.

    Raises:
        ValueError: If lengths are unequal.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length to compute Hamming distance.")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def format_mismatch_highlight(target: str, query: str) -> str:
    """Format mismatch comparison visually, e.g., AT[G->C]CA.

    Args:
        target: Sub-sequence window from the genome.
        query: Pattern query string.

    Returns:
        String with mismatched positions highlighted.
    """
    if len(target) != len(query):
        return f"{target} (length mismatch)"

    result = []
    for t_char, q_char in zip(target, query):
        if t_char == q_char:
            result.append(t_char)
        else:
            result.append(f"[{t_char}->{q_char}]")
    return "".join(result)


def exact_search(sequence: str, pattern: str) -> List[Dict[str, Any]]:
    """Find all exact occurrences of a pattern in a DNA sequence.

    Args:
        sequence: Main DNA sequence.
        pattern: Target motif/pattern.

    Returns:
        List of match dictionaries containing indices and coordinates.

    Time Complexity:
        O(N * M) where N is sequence length and M is pattern length.
    Space Complexity:
        O(K) where K is number of matches.
    """
    seq = validate_sequence(sequence)
    pat = validate_sequence(pattern)
    m = len(pat)
    matches = []

    if m == 0 or m > len(seq):
        return matches

    for start_idx, window in sliding_window(seq, k=m, step=1):
        if window == pat:
            matches.append({
                "start_index_0based": start_idx,
                "end_index_0based": start_idx + m,
                "start_coordinate_1based": start_idx + 1,
                "end_coordinate_1based": start_idx + m,
                "matched_sequence": window,
                "mismatches": 0,
                "visualization": window,
                "strand": "forward"
            })
    return matches


def fuzzy_search(sequence: str, pattern: str, max_mismatches: int = 1) -> List[Dict[str, Any]]:
    """Find fuzzy occurrences of a pattern within a maximum mismatch threshold (Hamming distance).

    Simulates Single Nucleotide Polymorphisms (SNPs) or point mutations.

    Args:
        sequence: Main DNA sequence.
        pattern: Target motif/pattern.
        max_mismatches: Maximum allowed mismatches (default 1).

    Returns:
        List of match dictionaries containing mismatch counts and visual highlights.

    Time Complexity:
        O(N * M) where N is sequence length and M is pattern length.
    Space Complexity:
        O(K) where K is number of matches.
    """
    seq = validate_sequence(sequence)
    pat = validate_sequence(pattern)
    m = len(pat)
    matches = []

    if m == 0 or m > len(seq):
        return matches

    for start_idx, window in sliding_window(seq, k=m, step=1):
        dist = hamming_distance(window, pat)
        if dist <= max_mismatches:
            vis = format_mismatch_highlight(window, pat) if dist > 0 else window
            matches.append({
                "start_index_0based": start_idx,
                "end_index_0based": start_idx + m,
                "start_coordinate_1based": start_idx + 1,
                "end_coordinate_1based": start_idx + m,
                "matched_sequence": window,
                "mismatches": dist,
                "visualization": vis,
                "strand": "forward"
            })
    return matches


def double_stranded_search(
    sequence: str, pattern: str, max_mismatches: int = 0
) -> List[Dict[str, Any]]:
    """Search for pattern on both forward and reverse-complement strands.

    Args:
        sequence: Main DNA sequence (forward strand).
        pattern: Target motif/pattern.
        max_mismatches: Maximum allowed mismatches.

    Returns:
        Combined list of match dictionaries with strand designation.
    """
    seq = validate_sequence(sequence)
    pat = validate_sequence(pattern)
    
    if max_mismatches == 0:
        forward_matches = exact_search(seq, pat)
    else:
        forward_matches = fuzzy_search(seq, pat, max_mismatches=max_mismatches)

    rev_seq = reverse_complement(seq)
    seq_len = len(seq)
    
    if max_mismatches == 0:
        rev_matches = exact_search(rev_seq, pat)
    else:
        rev_matches = fuzzy_search(rev_seq, pat, max_mismatches=max_mismatches)

    adjusted_rev_matches = []
    for match in rev_matches:
        rc_start_0 = match["start_index_0based"]
        rc_end_0 = match["end_index_0based"]
        orig_end_0 = seq_len - rc_start_0
        orig_start_0 = seq_len - rc_end_0

        adjusted_match = match.copy()
        adjusted_match["start_index_0based"] = orig_start_0
        adjusted_match["end_index_0based"] = orig_end_0
        adjusted_match["start_coordinate_1based"] = orig_start_0 + 1
        adjusted_match["end_coordinate_1based"] = orig_end_0
        adjusted_match["strand"] = "reverse"
        adjusted_rev_matches.append(adjusted_match)

    return forward_matches + adjusted_rev_matches
