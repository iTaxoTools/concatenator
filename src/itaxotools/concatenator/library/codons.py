#!/usr/bin/env python3

from typing import Dict, Optional, Iterator, Tuple, Set, TypeVar, Iterable, List, Any
import json
import sys
import os
import regex
import itertools
from collections import Counter
from enum import IntEnum
from dataclasses import dataclass

import pandas as pd

from .resources import get_resource


class ReadingFrame(IntEnum):

    def __new__(cls, value: int, text: str):
        obj = int.__new__(cls, value)  # type: ignore
        obj._value_ = value
        obj.text = text
        return obj

    def __str__(self):
        return f'{self.label} ({self.text})'

    @property
    def label(self):
        return f'{self.value:+}'

    Unknown = 0, 'unknown'

    P1 = +1, 'starts with 1st codon position'
    P2 = +2, 'starts with 2nd codon position'
    P3 = +3, 'starts with 3rd codon position'
    N1 = -1, 'reverse complement of +1'
    N2 = -2, 'reverse complement of +2'
    N3 = -3, 'reverse complement of +3'


@dataclass(frozen=True)
class GeneticCodeDescription:
    name: str
    ncbieaa: str
    sncbieaa: str
    abbr_name: Optional[str] = None


def validate_genetic_code_table(table: Any) -> None:
    if not isinstance(table, dict):
        raise ValueError(
            "In genetic_codes.json, in 'Genetic_code_tables' "
            "one of the values is not a dictionary")
    try:
        for field, field_type in (("name", str),
                                  ("id", int), ("ncbieaa", str), ("sncbieaa", str)):
            if not isinstance(table[field], field_type):
                raise ValueError(
                    "In genetic_codes.json, in 'Genetic_code_tables' "
                    f"the field '{field}' of one of the values is not an {field_type}"
                )
    except KeyError as e:
        raise ValueError(
            "In genetic_codes.json, in 'Genetic_code_tables' "
            f"the field '{e.args}' of one of the values is missing"
        ) from e
    for data_field in ("ncbieaa", "sncbieaa"):
        if len(table[data_field]) != 64:
            raise ValueError(
                "In genetic_codes.json, in 'Genetic_code_tables' "
                f"the field '{data_field}' of all of the values should be 64 bytes"
            )
    if "abbr_name" in table:
        if not isinstance(table["abbr_name"], str):
            raise ValueError(
                "In genetic_codes.json, in 'Genetic_code_tables' "
                f"the optional field 'abbr_name' of one of the values is not an str"
            )


def load_genetic_codes() -> Dict[int, GeneticCodeDescription]:
    genetic_codes_path = get_resource("genetic_codes.json")
    with open(genetic_codes_path) as genetic_codes_file:
        try:
            json_content = json.load(genetic_codes_file)
        except json.JSONDecodeError as e:
            raise ValueError("genetic_codes.json is not a valid JSON") from e
    if not isinstance(json_content, dict):
        raise ValueError("genetic_codes.json is not a dictionary")
    try:
        tables = json_content["Genetic_code_tables"]
    except KeyError as e:
        raise ValueError(
            "genetic_codes.json is missing 'Genetic_code_tables' field") from e
    if not isinstance(tables, list):
        raise ValueError(
            "In genetic_codes.json 'Genetic_code_tables' is not an array") from e
    result: Dict[int, GeneticCodeDescription] = {}
    for table in tables:
        validate_genetic_code_table(table)
        assert isinstance(table, dict)
        gc_id: int = table.pop("id")
        result[gc_id] = GeneticCodeDescription(**table)
    return result


_GC_DESCRIPTIONS = load_genetic_codes()


class _GeneticCodePrototype(int):

    def __new__(cls, value: int):
        obj = int.__new__(cls, value)  # type: ignore
        obj._value_ = value
        return obj


class GeneticCode(IntEnum):

    def __new__(cls, text: str, stops: List[str]):
        value = len(cls.__members__)
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.text = text
        obj.stops = stops
        return obj

    Unknown = 'Unknown', []

    SGC0 = 'Standard', ['TAA', 'TAG', 'TGA']
    SGC1 = 'Vertebrate Mitochondrial', ['TAA', 'TAG', 'AGA', 'AGG']
    SGC2 = 'Yeast Mitochondrial', ['TAA', 'TAG', 'TGA']
    SGC3 = 'Mold/Protozoan/Coelenterate Mitochondrial; Mycoplasma; Spiroplasma', [
        'TAA', 'TAG']
    ...


class GeneData:
    def __init__(self, series: pd.Series):
        self.series: pd.Series = series
        self.genetic_code: GeneticCode = GeneticCode(0)
        self.reading_frame: ReadingFrame = ReadingFrame(0)
        self.codons: Tuple[str, str, str] = ('**_1st', '**_2nd', '**_3rd')
        self.missing: str = 'N?'
        self.gap: str = '-'


class BadReadingFrame(Exception):
    pass


# Placeholder, needs to be implemented, in library.codons?
def some_func(data: GeneData) -> GeneData:
    return GeneData(data.series)


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
