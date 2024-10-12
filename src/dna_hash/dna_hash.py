import math
from typing import Iterator

import pybloomer

UP_BIT = 1

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

class DNAHash(object):
    counts = {
        # 'ACTG': 0
    }

    num_singletons = 0

    @staticmethod
    def _encode(sequence: str) -> int:
        """Encode a sequence using the up2bit representation."""
        hash = UP_BIT

        for i in range(len(sequence) - 1, -1, -1):
            base = sequence[i]

            hash <<= 2

            hash += BASE_ENCODE_MAP[base]

        return hash

    @staticmethod
    def _decode(hash: int) -> str:
        """Decode an up2bit hash into a sequence."""
        sequence = ''

        for i in range(0, int(math.log(hash, 2)), 2):
            base = (hash >> i) & 3

            sequence += BASE_DECODE_MAP[base]

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

    @property
    def num_non_singletons(self) -> int:
        """Return the total number of non-singleton sequences counted so far."""
        return sum(self.counts.values())

    @property
    def num_sequences(self) -> int:
        """Return the total number of sequences counted so far."""
        return self.num_non_singletons + self.num_singletons

    @property
    def num_unique_sequences(self) -> int:
        """Return the number of unique sequences stored in the hash table."""
        return len(self.counts) + self.num_singletons

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
        return max(self.counts)

    def argmax(self) -> str:
        """Return the sequence with the highest count."""
        return max(self.counts, key=self.counts.get)

    def top(self, k: int = 10) -> Iterator:
        """ Return the k sequences with the highest counts."""
        counts = dict(sorted(self.counts.items(), key=lambda count: count[1], reverse=True))

        n = 0

        for sequence, count in counts.items():
            sequence = self._decode(sequence)

            yield (sequence, count)
            n += 1

            if n >= k:
                break


    def __getitem__(self, sequence: str) -> int:
        exists = self.filter.exists(sequence)
        
        if not exists:
            raise ValueError('Sequence not found in hash table.')

        hash = self._encode(sequence)

        if hash in self.counts:
            return self.counts[hash]

        return 1

    def __setitem__(self, sequence: str, count: int) -> None:
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
        
    def __len__(self) -> int:
        return self.num_sequences