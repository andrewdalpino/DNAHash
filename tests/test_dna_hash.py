import random

from unittest import TestCase

from dna_hash import DNAHash


class TestDNAHash(TestCase):
    BASES = ["A", "C", "T", "G"]

    @classmethod
    def random_read(cls, k: int) -> str:
        return "".join(cls.BASES[random.randint(0, 3)] for i in range(0, k))

    def test_increment(self):
        hash_table = DNAHash()

        self.assertEqual(hash_table.num_singletons, 0)
        self.assertEqual(hash_table.num_sequences, 0)
        self.assertEqual(hash_table.num_unique_sequences, 0)

        hash_table.increment("ACTG")

        self.assertEqual(hash_table.num_singletons, 1)
        self.assertEqual(hash_table.num_sequences, 1)
        self.assertEqual(hash_table.num_unique_sequences, 1)
        self.assertEqual(hash_table["ACTG"], 1)

        hash_table.increment("ACTG")

        self.assertEqual(hash_table.num_singletons, 0)
        self.assertEqual(hash_table.num_sequences, 2)
        self.assertEqual(hash_table.num_unique_sequences, 1)
        self.assertEqual(hash_table["ACTG"], 2)

        self.assertEqual(hash_table.max(), 2)
        self.assertEqual(hash_table.argmax(), "ACTG")

    def test_top_k(self):
        hash_table = DNAHash()

        hash_table["CTGA"] = 1
        hash_table["ACTG"] = 10
        hash_table["GCGC"] = 4
        hash_table["AAAA"] = 9
        hash_table["AAAT"] = 2

        top = list(hash_table.top(3))

        self.assertEqual(len(top), 3)

        self.assertEqual(top[0], ("ACTG", 10))
        self.assertEqual(top[1], ("AAAA", 9))
        self.assertEqual(top[2], ("GCGC", 4))

    def test_large_dataset(self):
        random.seed(1)

        hash_table = DNAHash()

        for i in range(0, 100000):
            hash_table.increment(self.random_read(8))

        self.assertEqual(hash_table.num_sequences, 100000)
        self.assertEqual(hash_table.num_singletons, 21780)
        self.assertEqual(hash_table.num_unique_sequences, 51243)

        self.assertEqual(hash_table.max(), 9)
        self.assertEqual(hash_table.argmax(), "AGACTAAA")

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

        def test_encode_decode(self):
            """Test the encode and decode functionality."""
            sequences = ["A", "C", "G", "T", "ACTG", "GATTACA", "AAAAAA", "GGGCCC"]

            for seq in sequences:
                encoded = DNAHash.encode(seq)
                decoded = DNAHash.decode(encoded)
                self.assertEqual(seq, decoded)

    def test_get_method(self):
        """Test the get method for sequence retrieval."""
        hash_table = DNAHash()

        # Test get on empty hash
        self.assertEqual(hash_table.get("ACGT"), 0)

        # Test get on singleton
        hash_table.increment("ACGT")
        self.assertEqual(hash_table.get("ACGT"), 1)

        # Test get on non-singleton
        hash_table.increment("ACGT")
        self.assertEqual(hash_table.get("ACGT"), 2)

        # Test get on another sequence with direct insertion
        hash_table.insert("GGGG", 5)
        self.assertEqual(hash_table.get("GGGG"), 5)

        # Test get on non-existent sequence
        self.assertEqual(hash_table.get("TTTT"), 0)

    def test_insert_method(self):
        """Test the insert method with various counts."""
        hash_table = DNAHash()

        # Insert with count 1
        hash_table.insert("ACGT", 1)
        self.assertEqual(hash_table.num_singletons, 1)
        self.assertEqual(hash_table.num_sequences, 1)

        # Insert with count > 1
        hash_table.insert("GGGG", 3)
        self.assertEqual(hash_table.num_singletons, 1)
        self.assertEqual(hash_table.num_non_singletons, 3)
        self.assertEqual(hash_table.num_sequences, 4)

        # Update existing singleton to non-singleton
        hash_table.insert("ACGT", 5)
        self.assertEqual(hash_table.num_singletons, 0)
        self.assertEqual(hash_table.num_non_singletons, 8)
        self.assertEqual(hash_table.num_sequences, 8)

        # Update existing non-singleton
        hash_table.insert("GGGG", 10)
        self.assertEqual(hash_table.num_non_singletons, 15)
        self.assertEqual(hash_table.num_sequences, 15)

    def test_invalid_inputs(self):
        """Test error handling for invalid inputs."""
        hash_table = DNAHash()

        # Test insert with count < 1
        with self.assertRaises(ValueError):
            hash_table.insert("ACGT", 0)

        with self.assertRaises(ValueError):
            hash_table.insert("ACGT", -1)

        with self.assertRaises(ValueError):
            hash_table.encode("ACGTN")

    def test_dictionary_interface(self):
        """Test the dictionary-like behavior."""
        hash_table = DNAHash()

        # Test __setitem__
        hash_table["ACGT"] = 3
        self.assertEqual(hash_table.get("ACGT"), 3)

        # Test __getitem__
        self.assertEqual(hash_table["ACGT"], 3)
        self.assertEqual(hash_table["TTTT"], 0)  # Non-existent sequence

        # Test __len__
        self.assertEqual(len(hash_table), 1)

        hash_table["GGGG"] = 1
        self.assertEqual(len(hash_table), 2)

        hash_table.increment("AAAA")
        self.assertEqual(len(hash_table), 3)

    def test_empty_sequence_handling(self):
        """Test handling of empty sequences."""
        hash_table = DNAHash()

        # Test encode/decode with empty sequence
        encoded = DNAHash.encode("")
        self.assertEqual(DNAHash.decode(encoded), "")

        # Test insert and get with empty sequence
        hash_table.insert("", 1)
        self.assertEqual(hash_table.get(""), 1)

        hash_table.increment("")
        self.assertEqual(hash_table.get(""), 2)
