import unittest
import random

import dna_hash

class TestDNAHash(unittest.TestCase):
    BASES = ['A', 'C', 'T', 'G']

    @classmethod
    def random_read(cls, k: int) -> str:
        return ''.join(cls.BASES[random.randint(0, 3)] for i in range(0, k))

    def test_increment(self):
        hash_table = dna_hash.DNAHash()

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

        self.assertEqual(hash_table.max(), 2)
        self.assertEqual(hash_table.argmax(), 'ACTG')

    def test_top_k(self):
        hash_table = dna_hash.DNAHash()

        hash_table['CTGA'] = 1
        hash_table['ACTG'] = 10
        hash_table['GCGC'] = 4
        hash_table['AAAA'] = 9
        hash_table['AAAT'] = 2

        top = list(hash_table.top(3))

        self.assertEqual(len(top), 3)

        self.assertEqual(top[0], ('ACTG', 10))
        self.assertEqual(top[1], ('AAAA', 9))
        self.assertEqual(top[2], ('GCGC', 4))

    def test_advanced(self):
        random.seed(1)

        hash_table = dna_hash.DNAHash()

        for i in range(0, 100000):
            hash_table.increment(self.random_read(8))

        self.assertEqual(hash_table.num_sequences, 100000)
        self.assertEqual(hash_table.num_singletons, 21780)
        self.assertEqual(hash_table.num_unique_sequences, 51243)

        self.assertEqual(hash_table.max(), 9)
        self.assertEqual(hash_table.argmax(), 'AGACTAAA')
        