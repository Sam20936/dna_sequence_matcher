"""Unit tests for DNA Sequence Matcher package."""

import unittest
from pathlib import Path
import tempfile

from dna_matcher.validator import validate_sequence, parse_fasta, InvalidSequenceError
from dna_matcher.processor import sliding_window, reverse_complement, calculate_metrics
from dna_matcher.matcher import exact_search, fuzzy_search, double_stranded_search


class TestValidator(unittest.TestCase):

    def test_validate_sequence_valid(self):
        seq = "atcgN ATCG"
        self.assertEqual(validate_sequence(seq), "ATCGNATCG")

    def test_validate_sequence_invalid(self):
        with self.assertRaises(InvalidSequenceError):
            validate_sequence("ATCGX")
        with self.assertRaises(InvalidSequenceError):
            validate_sequence("")

    def test_parse_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "test.fasta"
            fasta_path.write_text(">seq1\nATCG\n>seq2\nGGCC\n", encoding="utf-8")
            data = parse_fasta(fasta_path)
            self.assertIn("seq1", data)
            self.assertEqual(data["seq1"], "ATCG")
            self.assertEqual(data["seq2"], "GGCC")


class TestProcessor(unittest.TestCase):

    def test_sliding_window(self):
        windows = list(sliding_window("ATCG", k=2, step=1))
        self.assertEqual(windows, [(0, "AT"), (1, "TC"), (2, "CG")])

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATCG"), "CGAT")
        self.assertEqual(reverse_complement("AAATTT"), "AAATTT")

    def test_calculate_metrics(self):
        metrics = calculate_metrics("GCGCATAT")
        self.assertEqual(metrics["length"], 8)
        self.assertEqual(metrics["gc_content"], 50.0)
        self.assertEqual(metrics["base_counts"]["G"], 2)


class TestMatcher(unittest.TestCase):

    def test_exact_search(self):
        seq = "AAATCGATCGATCG"
        matches = exact_search(seq, "ATCG")
        self.assertEqual(len(matches), 3)
        self.assertEqual(matches[0]["start_coordinate_1based"], 4)

    def test_fuzzy_search(self):
        seq = "AAATCGATCG"
        matches = fuzzy_search(seq, "ATCC", max_mismatches=1)
        self.assertEqual(len(matches), 3)
        self.assertEqual(matches[0]["mismatches"], 1)

    def test_double_stranded_search(self):
        seq = "ATCG"
        matches = double_stranded_search(seq, "ATCG", max_mismatches=0)
        self.assertEqual(len(matches), 2)


if __name__ == "__main__":
    unittest.main()
