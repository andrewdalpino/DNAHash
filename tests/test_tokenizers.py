from unittest import TestCase

from dna_hash.tokenizers import Kmer, Canonical, Fragment


class TestKmer(TestCase):
    def test_tokenize(self):
        tokenizer = Kmer(k=6, skip_invalid=True)

        tokens = tokenizer.tokenize("CGGTTCAGCANG")

        expected = ["CGGTTC", "GGTTCA", "GTTCAG", "TTCAGC", "TCAGCA"]

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

    def test_tokenize_empty_string(self):
        tokenizer = Kmer(k=6, skip_invalid=True)

        tokens = list(tokenizer.tokenize(""))

        self.assertEqual(tokens, [])

    def test_tokenize_short_sequence(self):
        tokenizer = Kmer(k=6, skip_invalid=True)

        tokens = list(tokenizer.tokenize("ACGT"))

        self.assertEqual(tokens, [])

    def test_tokenize_with_skip_invalid_false(self):
        tokenizer = Kmer(k=3, skip_invalid=False)

        with self.assertRaises(ValueError):
            list(tokenizer.tokenize("ACNGT"))


class TestCanonical(TestCase):
    def test_tokenize(self):
        tokenizer = Canonical(Kmer(k=6, skip_invalid=True))

        tokens = tokenizer.tokenize("CGGTTCAGCANG")

        expected = ["CGGTTC", "GGTTCA", "CTGAAC", "GCTGAA", "TCAGCA"]

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

    def test_tokenize_all_complementary(self):
        tokenizer = Canonical(Kmer(k=3, skip_invalid=True))

        tokens = list(tokenizer.tokenize("GACTCAGT"))

        expected = ["GAC", "ACT", "CTC", "TCA", "CAG", "ACT"]

        self.assertEqual(tokens, expected)

    def test_nested_tokenizer_properties(self):
        nested_tokenizer = Kmer(k=4, skip_invalid=True)

        tokenizer = Canonical(nested_tokenizer)

        tokens = list(tokenizer.tokenize("ACGTNNNAAA"))

        self.assertEqual(tokens, ["ACGT"])

    def test_palindromic_sequence(self):
        tokenizer = Canonical(Kmer(k=6, skip_invalid=True))

        tokens = list(tokenizer.tokenize("TACGTA"))

        self.assertEqual(tokens, ["TACGTA"])


class TestFragment(TestCase):
    def test_tokenize(self):
        tokenizer = Fragment(n=4, skip_invalid=True)

        tokens = tokenizer.tokenize("CGGTTCAGCANGTAAT")

        expected = ["CGGT", "TCAG", "TAAT"]

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

    def test_tokenize_with_different_fragment_size(self):
        tokenizer = Fragment(n=3, skip_invalid=True)

        tokens = list(tokenizer.tokenize("ACGTACGT"))

        expected = ["ACG", "TAC", "GT"]

        self.assertEqual(tokens, expected)

    def test_tokenize_with_invalid_characters(self):
        tokenizer = Fragment(n=2, skip_invalid=True)

        tokens = list(tokenizer.tokenize("ACNNTGXX"))

        expected = ["AC", "TG"]

        self.assertEqual(tokens, expected)

    def test_tokenize_with_skip_invalid_false(self):
        tokenizer = Fragment(n=3, skip_invalid=False)

        with self.assertRaises(ValueError):
            list(tokenizer.tokenize("ACNGTX"))
