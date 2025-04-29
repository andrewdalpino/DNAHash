import re

from abc import ABC, abstractmethod
from typing import Iterator

INVALID_BASE_REGEX = r"[^ACTG]"


class Tokenizer(ABC):
    """Base tokenizer class"""

    @abstractmethod
    def tokenize(self, sequence: str) -> Iterator[str]:
        pass


class Kmer(Tokenizer):
    """Generates tokens of length k from reads."""

    def __init__(self, k: int, skip_invalid: bool = False) -> None:
        if k < 1:
            raise ValueError(f"K cannot be less than 1, {k} given.")

        self._k = k
        self._skip_invalid = skip_invalid
        self._invalid_base = re.compile(INVALID_BASE_REGEX)

    def tokenize(self, sequence: str) -> Iterator[str]:
        """Tokenize the sequence."""

        i = 0

        while i <= len(sequence) - self._k:
            token = sequence[i : i + self._k]

            invalid_token = self._invalid_base.search(token)

            if invalid_token:
                if not self._skip_invalid:
                    offset = i + invalid_token.start()

                    raise ValueError(
                        f"Invalid base detected at offset {offset} in sequence."
                    )
                else:
                    i += 1 + invalid_token.start()

                    continue

            i += 1

            yield token


class Canonical(Tokenizer):
    """Tokenize sequences in their canonical form."""

    BASE_COMPLIMENT_MAP = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }

    @classmethod
    def reverse_complement(cls, sequence: str) -> str:
        """Return the reverse complement of a sequence."""

        complement = ""

        for base in reversed(sequence):
            if base not in cls.BASE_COMPLIMENT_MAP:
                raise ValueError("Invalid base {base} given.")

            complement += cls.BASE_COMPLIMENT_MAP[base]

        return complement

    def __init__(self, base: Tokenizer) -> None:
        self._base = base

    def tokenize(self, sequence: str) -> Iterator[str]:
        """Tokenize the sequence."""

        tokens = self._base.tokenize(sequence)

        for token in tokens:
            yield min(token, self.reverse_complement(token))


class Fragment(Tokenizer):
    """Generates a non-overlapping fragment of length n from a sequence."""

    def __init__(self, n: int, skip_invalid: bool = False) -> None:
        if n < 1:
            raise ValueError(f"N must be greater than 1, {n} given.")

        self._n = n
        self._skip_invalid = skip_invalid
        self._invalid_base = re.compile(INVALID_BASE_REGEX)

    def tokenize(self, sequence: str) -> Iterator[str]:
        """Tokenize the sequence."""

        m = len(sequence)

        if m < self._n:
            yield sequence

            return

        for i in range(0, m, self._n):
            token = sequence[i : i + self._n]

            invalid_token = self._invalid_base.search(token)

            if invalid_token:
                if not self._skip_invalid:
                    offset = i + invalid_token.start()

                    raise ValueError(
                        f"Invalid base detected at offset {offset} in sequence."
                    )
                else:
                    continue

            yield token
