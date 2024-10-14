import unittest

from dna_hash import tokenizers

TEST_SEQUENCES = [
    'TGAGCATTCCCGATATGGTTCTCATTCATGCAGGCACTAATGATATGTTGCAAAACTATTCCGTCAATCATGCTGTAGATCATATTCGAGAAACGATTGATGTGCTACGCCAAAGAAATCCAAATACAGATTGGCAATTACATGAATTTCGCAGGTATCTCAAAAAAAATGGATACCCAAATCATGCAAATACGTTAATGGATTTATGGAAAAATCCATATCAATCAAACGCAAAAGAAA',
    'TATCCAAGTGATTTACATGCAGATCTCATATAGCAAAGATGTTAAATAATACAAAAAATAACCACTTTTCATCGTATCAATATAGGTCCGTCAAGACAATTTTGAAATATTGATGCTACGGATTTATATGTTTTAATATTCTTTTAAGGAAAGTCTGCGGGAGTCCTTGAGGTAGAAAATTGAAATTATTTTGAGGAGAAAAAGAATTTTGAACGAAGGAGTTTTGCCGTAGTAAAGGGA',
    'AAGAAAAGTTTTCAATTCTTAGTGCAATAAATGCTCATACCAAAATTCTTTTTTTTAATTAAGCTGATCAATTTTTACAAAAAACAACATCGAATGTCAACATCTTTTTCTCCAATTAGAAAAAACAATCATCACCAGCTAACGTCATAATCGGTTTCGACATTTTTGACATATGAACTTTAAAGGTGATACCCTTTAGATACTCAGTTTATACTAAAAATAACAACCTTATATCTATAA',
]

class TestKmer(unittest.TestCase):
    def test_basic(self):
        tokenizer = tokenizers.Kmer(k=10)

        tokens = []

        for sequence in TEST_SEQUENCES:
            for token in tokenizer.tokenize(sequence):
                tokens.append(token)

        self.assertEqual(690, len(tokens))
        self.assertEqual(0, tokenizer.dropped)

    def test_tokenize(self):
        tokenizer = tokenizers.Kmer(k=6)

        tokens = tokenizer.tokenize('CGGTTCAGCANG')

        expected = ['CGGTTC', 'GGTTCA', 'GTTCAG', 'TTCAGC', 'TCAGCA']

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

class TestCanonical(unittest.TestCase):
    def test_basic(self):
        tokenizer = tokenizers.Canonical(tokenizers.Kmer(k=10))

        tokens = []

        for sequence in TEST_SEQUENCES:
            for token in tokenizer.tokenize(sequence):
                tokens.append(token)

        self.assertEqual(690, len(tokens))

    def test_tokenize(self):
        tokenizer = tokenizers.Canonical(tokenizers.Kmer(k=6))

        tokens = tokenizer.tokenize('CGGTTCAGCANG')

        expected = ['CGGTTC', 'GGTTCA', 'CTGAAC', 'GCTGAA', 'TCAGCA']

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])

class TestFragment(unittest.TestCase):
    def test_basic(self):
        tokenizer = tokenizers.Fragment(n=10)

        tokens = []

        for sequence in TEST_SEQUENCES:
            for token in tokenizer.tokenize(sequence):
                tokens.append(token)

        self.assertEqual(69, len(tokens))

    def test_tokenize(self):
        tokenizer = tokenizers.Fragment(n=4)

        tokens = tokenizer.tokenize('CGGTTCAGCANGTAAT')

        expected = ['CGGT', 'TCAG', 'TAAT']

        for i, token in enumerate(tokens):
            self.assertEqual(token, expected[i])