import re

from abc import ABC, abstractmethod
from typing import Iterator

INVALID_BASE_REGEX = r'[^ACTG]'

class Tokenizer(ABC):
    """Base tokenizer class"""

    @abstractmethod
    def tokenize(self, sequence: str) -> Iterator[str]:
        pass

class Kmer(Tokenizer):
    """Generates tokens of length k from reads."""
    
    def __init__(self, k: int, skip_invalid: bool = False) -> None:
        if k < 1:
            raise ValueError(f'K cannot be less than 1, {k} given.')

        self.k = k
        self.skip_invalid = skip_invalid
        self.invalid_base = re.compile(INVALID_BASE_REGEX)
        self.dropped = 0

    def tokenize(self, sequence: str) -> Iterator[str]:
        """Tokenize the sequence."""
        i = 0

        while i < len(sequence) - self.k:
            token = sequence[i:i + self.k]

            invalid_token = self.invalid_base.search(token)

            if invalid_token:
                if not self.skip_invalid:
                    offset = i + invalid_token.start()
                
                    raise ValueError('Invalid base detected at'
                        + f' offset {offset} in sequence.')

                else:
                    skip = 1 + invalid_token.start()

                    i += skip

                    self.dropped += skip

                    continue

            i += 1
            
            yield token

class Canonical(Tokenizer):
    """Tokenize sequences in their canonical form."""

    BASE_COMPLIMENT_MAP = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }

    @classmethod
    def reverse_complement(cls, sequence: str) -> str:
        """Return the reverse complement of a sequence."""
        complement = ''

        for i in range(len(sequence) - 1, -1, -1):
            base = sequence[i]

            if base not in cls.BASE_COMPLIMENT_MAP:
                raise ValueError('Invalid base {base} given.')

            complement += cls.BASE_COMPLIMENT_MAP[base]

        return complement

    def __init__(self, base: Tokenizer) -> None:
        self.base = base

    def tokenize(self, sequence: str) -> Iterator[str]:
        """Tokenize the sequence."""
        tokens = self.base.tokenize(sequence)

        for token in tokens:
            yield min(token, self.reverse_complement(token))

class Fragment(Tokenizer):
    """Generates a non-overlapping fragment of length n from a sequence."""

    def __init__(self, n: int, skip_invalid: bool = False) -> None:
        if n < 1:
            raise ValueError(f'N must be greater than 1, {n} given.')

        self.n = n
        self.skip_invalid = skip_invalid
        self.invalid_base = re.compile(INVALID_BASE_REGEX)
        self.dropped = 0

    def tokenize(self, sequence: str) -> Iterator[str]:
        """Tokenize the sequence."""
        m = len(sequence)

        if m < self.n:
            yield sequence

            return

        for i in range(0, m, self.n):
            token = sequence[i:i + self.n]

            invalid_token = self.invalid_base.search(token)

            if invalid_token:
                if not self.skip_invalid:
                    offset = i + invalid_token.start()
                
                    raise ValueError('Invalid base detected at'
                        + f' offset {offset} in sequence.')

                else:
                    self.dropped += 1

                    continue

            yield token
