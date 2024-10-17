from unittest import TestCase

from dna_hash.tokenizers import Kmer, Canonical, Fragment

class TestKmer(TestCase):
    def test_tokenize(self):
        tokenizer = Kmer(k=6, skip_invalid=True)

        tokens = tokenizer.tokenize('CGGTTCAGCANG')

        expected = ['CGGTTC', 'GGTTCA', 'GTTCAG', 'TTCAGC', 'TCAGCA']

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

        self.assertEqual(tokenizer.dropped, 6)

class TestCanonical(TestCase):
    def test_tokenize(self):
        tokenizer = Canonical(Kmer(k=6, skip_invalid=True))

        tokens = tokenizer.tokenize('CGGTTCAGCANG')

        expected = ['CGGTTC', 'GGTTCA', 'CTGAAC', 'GCTGAA', 'TCAGCA']

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

class TestFragment(TestCase):
    def test_tokenize(self):
        tokenizer = Fragment(n=4, skip_invalid=True)

        tokens = tokenizer.tokenize('CGGTTCAGCANGTAAT')

        expected = ['CGGT', 'TCAG', 'TAAT']

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

        self.assertEqual(tokenizer.dropped, 1)
        