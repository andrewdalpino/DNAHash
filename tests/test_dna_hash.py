import unittest

import dna_hash

class TestDNAHash(unittest.TestCase):
    def test_basic(self):
        hash_table = dna_hash.DNAHash(
            max_false_positive_rate=0.01,
            num_hashes=4,
            layer_size=32000000,
        )

        self.assertEqual(hash_table.num_singletons, 0)
        self.assertEqual(hash_table.num_sequences, 0)
        self.assertEqual(hash_table.num_unique_sequences, 0)

        hash_table.increment('ACTG')

        self.assertEqual(hash_table.num_singletons, 1)
        self.assertEqual(hash_table.num_sequences, 1)
        self.assertEqual(hash_table.num_unique_sequences, 1)
        self.assertEqual(hash_table['ACTG'], 1)

        hash_table.increment('ACTG')

        self.assertEqual(hash_table.num_singletons, 0)
        self.assertEqual(hash_table.num_sequences, 2)
        self.assertEqual(hash_table.num_unique_sequences, 1)
        self.assertEqual(hash_table['ACTG'], 2)

    def test_top_k(self):
        hash_table = dna_hash.DNAHash(
            max_false_positive_rate=0.01,
            num_hashes=4,
            layer_size=32000000,
        )

        hash_table['CTGA'] = 1
        hash_table['ACTG'] = 10
        hash_table['GCGC'] = 4
        hash_table['AAAA'] = 9

        top = list(hash_table.top(3))

        self.assertEqual(top[0], ('ACTG', 10))
        self.assertEqual(top[1], ('AAAA', 9))
        self.assertEqual(top[2], ('GCGC', 4))
