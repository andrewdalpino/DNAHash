import unittest

from dna_hash import tokenizers

class TestKmer(unittest.TestCase):
    def test_tokenize(self):
        tokenizer = tokenizers.Kmer(k=6)

        tokens = tokenizer.tokenize('CGGTTCAGCANG')

        expected = ['CGGTTC', 'GGTTCA', 'GTTCAG', 'TTCAGC', 'TCAGCA']

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

        self.assertEqual(tokenizer.dropped, 6)

class TestCanonical(unittest.TestCase):
    def test_tokenize(self):
        tokenizer = tokenizers.Canonical(tokenizers.Kmer(k=6))

        tokens = tokenizer.tokenize('CGGTTCAGCANG')

        expected = ['CGGTTC', 'GGTTCA', 'CTGAAC', 'GCTGAA', 'TCAGCA']

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

class TestFragment(unittest.TestCase):
    def test_tokenize(self):
        tokenizer = tokenizers.Fragment(n=4)

        tokens = tokenizer.tokenize('CGGTTCAGCANGTAAT')

        expected = ['CGGT', 'TCAG', 'TAAT']

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

        self.assertEqual(tokenizer.dropped, 1)
        