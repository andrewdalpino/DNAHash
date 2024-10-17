import math
import sys

from typing import Iterator, Tuple

import okbloomer

from dna_hash.tokenizers import Fragment

class DNAHash(object):
    """A specialized datastructure for counting DNA sequences for use in Bioinformatics."""

    UP_BIT = 1

    BASE_ENCODE_MAP = {
        'A': 0,
        'C': 1,
        'T': 2,
        'G': 3,
    }

    BITS_PER_BASE = max(BASE_ENCODE_MAP.values()).bit_length()

    MAX_FRAGMENT_LENGTH = math.ceil(math.log(sys.maxsize, 2) / BITS_PER_BASE) - UP_BIT.bit_length()

    BASE_DECODE_MAP = {encoding: base for base, encoding in BASE_ENCODE_MAP.items()}

    def __init__(self,
                max_false_positive_rate: float = 0.01,
                num_hashes: int = 4,
                layer_size: int = 32000000) -> None:

        self.filter = okbloomer.BloomFilter(
            max_false_positive_rate=max_false_positive_rate,
            num_hashes=num_hashes,
            layer_size=layer_size,
        )

        self.counts: dict[Tuple[int, ...], int] = {}
        self.num_singletons = 0
        self.tokenizer = Fragment(self.MAX_FRAGMENT_LENGTH)

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

    def insert(self, sequence: str, count: int) -> None:
        """Insert a sequence count into the hash table."""
        if count < 1:
            raise ValueError(f'Count cannot be less than 1, {count} given.')

        exists = self.filter.exists_or_insert(sequence)

        if count > 1:
            hashes = self._encode(sequence)

            if exists and hashes not in self.counts:
                self.num_singletons -= 1

            self.counts[hashes] = count
        
        elif not exists:
            self.num_singletons += 1

    def increment(self, sequence: str) -> None:
        """Increment the count for a given sequence by 1."""
        exists = self.filter.exists_or_insert(sequence)

        if exists:
            hashes = self._encode(sequence)

            if hashes in self.counts:
                self.counts[hashes] += 1
            
            else:
                self.num_singletons -= 1

                self.counts[hashes] = 2

        else:
            self.num_singletons += 1

    def max(self) -> int:
        """Return the highest sequence count."""
        return max(self.counts.values())

    def argmax(self) -> str:
        """Return the sequence with the highest count."""
        hashes, count = max(self.counts.items(), key=lambda item: item[1])

        return self._decode(hashes)

    def get(self, sequence: str) ->int:
        """Return the count for a sequence."""
        exists = self.filter.exists(sequence)
        
        if not exists:
            raise ValueError('Sequence not found in hash table.')

        hashes = self._encode(sequence)

        if hashes in self.counts:
            return self.counts[hashes]

        return 1

    def top(self, k: int = 10) -> Iterator:
        """ Return the k sequences with the highest counts."""
        counts = sorted(self.counts.items(), key=lambda item: item[1], reverse=True)

        for hashes, count in counts[0:k]:
            sequence = self._decode(hashes)

            yield (sequence, count)

    def _encode(self, sequence: str) -> Tuple[int, ...]:
        """Encode a variable-length sequence using the up2bit representation."""
        hashes = []

        for fragment in self.tokenizer.tokenize(sequence):
            n = len(fragment)

            hash = self.UP_BIT

            for i in range(n - 1, -1, -1):
                base = fragment[i]

                hash <<= 2
                hash += self.BASE_ENCODE_MAP[base]

            hashes.append(hash)

        return tuple(hashes)

    def _decode(self, hashes: Tuple[int, ...]) -> str:
        """Decode an up2bit representation into a variable-length sequence."""
        sequence = ''
        
        for hash in hashes:
            if hash == self.UP_BIT:
                continue

            for i in range(0, int(math.log(hash, 2)), 2):
                encoding = (hash >> i) & 3

                sequence += self.BASE_DECODE_MAP[encoding]

        return sequence

    def __setitem__(self, sequence: str, count: int) -> None:
        self.insert(sequence, count)

    def __getitem__(self, sequence: str) -> int:
        return self.get(sequence)
        
    def __len__(self) -> int:
        return self.num_sequences
        