#!/usr/bin/env python3

from typing import Dict, Optional, Iterator, Tuple, Set, TypeVar, Iterable, List
import json
import sys
import os
import regex
import itertools
from collections import Counter
import pandas as pd

from .resources import get_resource

stop_codons_path = get_resource("stop_codons.json")
try:
    stop_codons: Optional[Dict[str, Set[int]]] = {
        codon: set(tables)
        for codon, tables in json.load(open(stop_codons_path)).items()
    }
except json.JSONDecodeError:
    stop_codons = None
else:
    assert stop_codons is not None
    stop_codons_set = set(stop_codons.keys()) if stop_codons else None

    table_set = set(table for tables in stop_codons.values() for table in tables)


def detect_stop_codons(
    sequence: str, stop_codons: Set[str]
) -> Iterator[Tuple[str, int]]:
    """
    Iterator over forward-sense stop codons in sequence.

    Yields codon with a reading frame (1, 2 or 3)

    Ignores stop codons in the end
    """
    stops_regex = regex.compile("|".join(stop_codons), regex.IGNORECASE)
    sequence_length = len(sequence)
    for stop_match in stops_regex.finditer(sequence, overlapped=True):
        if sequence_length - stop_match.start() < 6:
            continue
        frame = stop_match.start() % 3 + 1
        yield stop_match.group(), frame


def detect_reverse_stop_codons(
    sequence: str, stop_codons: Set[str]
) -> Iterator[Tuple[str, int]]:
    """
    Iterator over reverse-sense stop codons in sequence.

    Yields codon with a reading frame (-1, -2 or -3)

    Ignores stop codons in the end (from the reverse perspective, i.e. the beginning)
    """
    stops_regex = regex.compile(
        "|".join(codon[::-1] for codon in stop_codons), regex.IGNORECASE
    )
    end_frame = len(sequence)
    for stop_match in stops_regex.finditer(sequence, overlapped=True):
        if stop_match.start() < 3:
            continue
        frame = -((end_frame - stop_match.start()) % 3) - 1
        yield stop_match.group()[::-1], frame


T = TypeVar("T")


def collect_non_unique(iter: Iterable[T]) -> Set[T]:
    seen: Set[T] = set()
    result: Set[T] = set()
    for x in iter:
        if x in seen:
            result.add(x)
        seen.add(x)
    return result


def allowed_translations(
    stops: Set[Tuple[str, int]], stop_codons: Dict[str, Set[int]]
) -> Dict[int, Set[int]]:
    """
    Returns mapping from reading frame to possible translation tables
    """
    result: Dict[int, Set[int]] = {
        frame: table_set.copy() for frame in [-3, -2, -1, 1, 2, 3]
    }
    for stop_codon, frame in stops:
        result[frame] -= stop_codons[stop_codon]
    return result


def choose_frame(frame_to_table: Dict[int, Set[int]]) -> List[int]:
    return [frame for frame, tables in frame_to_table.items() if tables]


def detect_reading_frame(sequence: str) -> Optional[List[int]]:
    if not stop_codons_set:
        return None
    stops = set(
        itertools.chain(
            detect_stop_codons(sequence, stop_codons_set),
            detect_reverse_stop_codons(sequence, stop_codons_set),
        )
    )
    assert stop_codons is not None
    return choose_frame(allowed_translations(stops, stop_codons))


def column_reading_frame(column: pd.Series) -> List[int]:
    """
    Returns the list of reading frames detected in `column`.

    The result is ordered by decreasing occurences
    """
    frame_counter = Counter()
    for _, seq in column.items():
        frame_counter.update(detect_reading_frame(seq))
    return [frame for frame, _ in frame_counter.most_common()]


def split_codon_charsets(
    column: pd.Series, frame: int
) -> Tuple[pd.Series, pd.Series, pd.Series]:
    starts = {
        1: (0, 1, 2),
        2: (2, 0, 1),
        3: (1, 2, 0),
        -1: (-1, -2, -3),
        -2: (-3, -1, -2),
        -3: (-2, -3, -1),
    }[frame]
    if frame > 0:
        return tuple([column.str[start::3] for start in starts])
    else:
        return tuple([column.str[start::-3] for start in starts])