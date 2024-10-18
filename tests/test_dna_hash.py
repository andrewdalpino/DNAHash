import random

from unittest import TestCase

from dna_hash import DNAHash

class TestDNAHash(TestCase):
    BASES = ['A', 'C', 'T', 'G']

    @classmethod
    def random_read(cls, k: int) -> str:
        return ''.join(cls.BASES[random.randint(0, 3)] for i in range(0, k))

    def test_increment(self):
        hash_table = DNAHash()

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
        hash_table = DNAHash()

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

    def test_large_dataset(self):
        random.seed(1)

        hash_table = DNAHash()

        for i in range(0, 100000):
            hash_table.increment(self.random_read(8))

        self.assertEqual(hash_table.num_sequences, 100000)
        self.assertEqual(hash_table.num_singletons, 21780)
        self.assertEqual(hash_table.num_unique_sequences, 51243)

        self.assertEqual(hash_table.max(), 9)
        self.assertEqual(hash_table.argmax(), 'AGACTAAA')

    def test_long_sequences(self):
        random.seed(1)

        hash_table = DNAHash()

        sequence = self.random_read(500)

        hash_table.insert(sequence, 420)

        self.assertEqual(hash_table.num_sequences, 420)
        self.assertEqual(hash_table.num_singletons, 0)
        self.assertEqual(hash_table.num_unique_sequences, 1)

        argmax = hash_table.argmax()

        self.assertEqual(argmax, sequence)
        self.assertEqual(len(argmax), 500)
        