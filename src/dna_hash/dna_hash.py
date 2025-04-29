import math

from typing import Iterator

from okbloomer import BloomFilter


class DNAHash(object):
    """
    A specialized datastructure for counting genetic sequences for use in Machine
    Learning and Bioinformatics.
    """

    UP_BIT = 1

    BASE_ENCODE_MAP = {
        "A": 0,
        "C": 1,
        "T": 2,
        "G": 3,
    }

    BASE_DECODE_MAP = {encoding: base for base, encoding in BASE_ENCODE_MAP.items()}

    @classmethod
    def encode(cls, sequence: str) -> int:
        """Encode a variable-length sequence using the up2bit representation."""

        hash = cls.UP_BIT

        for i in range(len(sequence) - 1, -1, -1):
            base = sequence[i]

            if base not in cls.BASE_ENCODE_MAP:
                raise ValueError(f"Invalid base '{base}' in sequence '{sequence}'.")

            hash <<= 2
            hash += cls.BASE_ENCODE_MAP[base]

        return hash

    @classmethod
    def decode(cls, hash: int) -> str:
        """Decode an up2bit representation into a variable-length sequence."""

        sequence = ""

        for i in range(0, int(math.log(hash, 2)), 2):
            encoding = (hash >> i) & 3

            sequence += cls.BASE_DECODE_MAP[encoding]

        return sequence

    def __init__(
        self,
        max_false_positive_rate: float = 0.01,
        num_hashes: int = 4,
        layer_size: int = 32000000,
    ) -> None:
        self._filter = BloomFilter(
            max_false_positive_rate=max_false_positive_rate,
            num_hashes=num_hashes,
            layer_size=layer_size,
        )

        self._counts: dict[int, int] = {}
        self._num_singletons = 0

    @property
    def counts(self) -> dict[int, int]:
        return self._counts

    @property
    def num_sequences(self) -> int:
        """Return the total number of sequences counted so far."""

        return self.num_non_singletons + self.num_singletons

    @property
    def num_unique_sequences(self) -> int:
        """Return the number of unique sequences stored in the hash table."""

        return len(self._counts) + self._num_singletons

    @property
    def num_singletons(self) -> int:
        return self._num_singletons

    @property
    def num_non_singletons(self) -> int:
        """Return the total number of non-singleton sequences counted so far."""

        return sum(self._counts.values())

    def insert(self, sequence: str, count: int) -> None:
        """Insert a sequence count into the hash table."""

        if count < 1:
            raise ValueError(f"Count cannot be less than 1, {count} given.")

        exists = self._filter.exists_or_insert(sequence)

        if count > 1:
            hash = self.encode(sequence)

            if exists and hash not in self.counts:
                self._num_singletons -= 1

            self._counts[hash] = count
        elif not exists:
            self._num_singletons += 1

    def increment(self, sequence: str) -> None:
        """Increment the count for a given sequence by 1."""

        exists = self._filter.exists_or_insert(sequence)

        if exists:
            hash = self.encode(sequence)

            if hash in self._counts:
                self._counts[hash] += 1
            else:
                self._num_singletons -= 1
                self._counts[hash] = 2
        else:
            self._num_singletons += 1

    def max(self) -> int:
        """Return the highest sequence count."""

        return max(self._counts.values())

    def argmax(self) -> str:
        """Return the sequence with the highest count."""

        hash, _ = max(self._counts.items(), key=lambda item: item[1])

        return self.decode(hash)

    def get(self, sequence: str) -> int:
        """Return the count for a sequence."""

        exists = self._filter.exists(sequence)

        if not exists:
            return 0

        hash = self.encode(sequence)

        if hash in self._counts:
            return self._counts[hash]

        return 1

    def top(self, k: int = 10) -> Iterator[tuple[str, int]]:
        """Return the k sequences with the highest counts."""

        counts = sorted(self._counts.items(), key=lambda item: item[1], reverse=True)

        for hash, count in counts[:k]:
            sequence = self.decode(hash)

            yield (sequence, count)

    def __setitem__(self, sequence: str, count: int) -> None:
        self.insert(sequence, count)

    def __getitem__(self, sequence: str) -> int:
        return self.get(sequence)

    def __len__(self) -> int:
        return self.num_unique_sequences
