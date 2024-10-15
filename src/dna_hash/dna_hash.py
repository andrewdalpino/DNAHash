import math
from typing import Iterator, Tuple
import sys

import numpy as np
from nptyping import NDArray

import pybloomer

class DNAHash(object):
    """A specialized datastructure for counting short DNA sequences for use in Bioinformatics."""

    UP_BIT = 1

    BITS_PER_BASE = 2

    BASE_ENCODE_MAP = {
        'A': 0,
        'C': 1,
        'T': 2,
        'G': 3,
    }

    BASE_DECODE_MAP = {
        0: 'A',
        1: 'C',
        2: 'T',
        3: 'G',
    }

    MAX_SEQUENCE_LENGTH = math.ceil(math.log(sys.maxsize, 2) / BITS_PER_BASE) - 1

    @classmethod
    def _encode(cls, sequence: str) -> int:
        """Encode a variable-length sequence using the up2bit representation."""
        m = len(sequence)

        if m > cls.MAX_SEQUENCE_LENGTH:
            raise ValueError('Sequence length must be less than'
                + f' {cls.MAX_SEQUENCE_LENGTH}, {sequence} given.')

        hash = cls.UP_BIT

        for i in range(m - 1, -1, -1):
            base = sequence[i]

            if base not in cls.BASE_ENCODE_MAP:
                raise ValueError(f'Invalid base {base} given in sequence {sequence}.')

            hash <<= 2

            hash += cls.BASE_ENCODE_MAP[base]

        return hash

    @classmethod
    def _decode(cls, hash: int) -> str:
        """Decode an up2bit hash into a variable-length sequence."""
        sequence = ''

        for i in range(0, int(math.log(hash, 2)), 2):
            base = (hash >> i) & 3

            sequence += cls.BASE_DECODE_MAP[base]

        return sequence

    def __init__(self,
                max_false_positive_rate: float = 0.01,
                num_hashes: int = 4,
                layer_size: int = 32000000) -> None:

        self.filter = pybloomer.BloomFilter(
            max_false_positive_rate=max_false_positive_rate,
            num_hashes=num_hashes,
            layer_size=layer_size,
        )

        self.counts: dict[int, int] = {}
        self.num_singletons = 0

    @property
    def num_sequences(self) -> int:
        """Return the total number of sequences counted so far."""
        return self.num_non_singletons + self.num_singletons

    @property
    def num_unique_sequences(self) -> int:
        """Return the number of unique sequences stored in the hash table."""
        return len(self.counts) + self.num_singletons

    @property
    def num_non_singletons(self) -> int:
        """Return the total number of non-singleton sequences counted so far."""
        return sum(self.counts.values())

    def insert(self, sequence: str, count: int = 1) -> None:
        """Insert a sequence count into the hash table."""
        if count < 1:
            raise ValueError(f'Count cannot be less than 1, {count} given.')

        exists = self.filter.exists_or_insert(sequence)

        if count > 1:
            hash = self._encode(sequence)

            if exists and hash not in self.counts:
                self.num_singletons -= 1

            self.counts[hash] = count
        
        elif not exists:
            self.num_singletons += 1

    def increment(self, sequence: str) -> None:
        """Increment the count for a given sequence by 1."""
        exists = self.filter.exists_or_insert(sequence)

        if exists:
            hash = self._encode(sequence)

            if hash in self.counts:
                self.counts[hash] += 1
            
            else:
                self.num_singletons -= 1

                self.counts[hash] = 2

        else:
            self.num_singletons += 1

    def max(self) -> int:
        """Return the highest sequence count."""
        return max(self.counts.values())

    def argmax(self) -> str:
        """Return the sequence with the highest count."""
        hash, count = max(self.counts.items(), key=lambda item: item[1])

        return self._decode(hash)

    def get(self, sequence: str) ->int:
        """Return the count for a sequence."""
        exists = self.filter.exists(sequence)
        
        if not exists:
            raise ValueError('Sequence not found in hash table.')

        hash = self._encode(sequence)

        if hash in self.counts:
            return self.counts[hash]

        return 1

    def top(self, k: int = 10) -> Iterator:
        """ Return the k sequences with the highest counts."""
        counts = sorted(self.counts.items(), key=lambda item: item[1], reverse=True)

        for hash, count in counts[0:k]:
            sequence = self._decode(hash)

            yield (sequence, count)

    def histogram(self, bins: int = 10) -> Tuple[NDArray, NDArray]:
        """Return a histogram of sequences bucketed by their counts."""
        return np.histogram(list(self.counts.values()), bins=bins)

    def __setitem__(self, sequence: str, count: int) -> None:
        self.insert(sequence, count)

    def __getitem__(self, sequence: str) -> int:
        return self.get(sequence)
        
    def __len__(self) -> int:
        return self.num_sequences
        