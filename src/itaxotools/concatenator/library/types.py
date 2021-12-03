
from typing import NamedTuple, Callable, List, Optional, Tuple
from collections import deque
from enum import Enum

import re


class Justification(Enum):
    Left = 'Left', str.ljust
    Right = 'Right', str.rjust
    Center = 'Center', str.center
    NoJust = 'None', None

    def __init__(self, description: str, method: Optional[Callable]):
        self.description = description
        self.method = method

    def apply(self, text: str, *args, **kwargs):
        if not self.method:
            return text
        return self.method(text, *args, **kwargs)


class TextCase(Enum):
    Unchanged = 'Unchanged', None
    Upper = 'Uppercase', str.upper
    Lower = 'Lowercase', str.lower

    def __init__(self, description: str, method: Optional[Callable]):
        self.description = description
        self.method = method

    def apply(self, text: str, *args, **kwargs):
        if not self.method:
            return text
        return self.method(text, *args, **kwargs)


class CodonSet(NamedTuple):
    name: str
    position: int
    position_end: int

    def __str__(self) -> str:
        return (
            f'{self.name} = {self.position}-{self.position_end}\\3;')


class Charset(NamedTuple):
    name: str
    position: int
    length: int
    frame: int
    codons: Tuple[str, str, str]

    def __str__(self) -> str:
        return (
            f'{self.name} = {self.position}-{self.position_end};')

    def codon_sets(self) -> List[CodonSet]:
        frame = self.frame
        codons = deque(
            re.sub(r'\*\*', self.name, codon)
            for codon in self.codons)
        if frame < 0:
            offset = self.length % 3
            frame = - frame
            codons.reverse()
            codons.rotate(offset)
        codons.rotate(frame - 1)
        result = list()
        for offset, codon in enumerate(codons):
            result.append(CodonSet(
                codon, self.position + offset, self.position_end))
        return result

    @property
    def position_end(self) -> int:
        return self.position + self.length - 1
