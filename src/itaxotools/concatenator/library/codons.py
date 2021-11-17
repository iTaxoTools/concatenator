#!/usr/bin/env python3

from __future__ import annotations
from typing import (Dict, Optional, Iterator, Tuple, Set,
                    TypeVar, Iterable, List, Any, Counter, DefaultDict)
import json
import sys
import os
import regex  # type: ignore
import itertools
from enum import IntEnum, Enum
from dataclasses import dataclass

import pandas as pd  # type: ignore

from .resources import get_resource


class ReadingFrame(IntEnum):

    def __new__(cls, value: int, text: str):
        obj = int.__new__(cls, value)  # type: ignore
        obj._value_ = value
        obj.text = text
        return obj

    def __str__(self):
        return f'{self.label} ({self.text})'

    @classmethod
    def from_int(cls, value: int) -> ReadingFrame:
        """
        Constructor for `ReadingFrame` that passes type-checking
        """
        return cls(value)  # type: ignore

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
class _GeneticCodeDescription:
    name: str
    ncbieaa: str
    sncbieaa: str
    abbr_name: Optional[str] = None

    # bases and codons in the same order as
    # in https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/gc.prt
    BASES = "TCAG"
    CODONS = tuple(c1 + c2 + c3
                   for c1, c2, c3 in itertools.product(BASES, repeat=3))

    def stops(self) -> List[str]:
        return [codon for seaa, codon in zip(self.sncbieaa, self.CODONS) if seaa == "*"]


def _validate_genetic_code_table(table: Any) -> None:
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
    if table["id"] == 0:
        raise ValueError("In genetic_codes.json, a table has an invalid 'id' 0")
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


def _load_genetic_codes() -> Dict[int, _GeneticCodeDescription]:
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
    result: Dict[int, _GeneticCodeDescription] = {}
    for table in tables:
        _validate_genetic_code_table(table)
        assert isinstance(table, dict)
        gc_id: int = table.pop("id")
        result[gc_id] = _GeneticCodeDescription(**table)
    return result


_GC_DESCRIPTIONS = _load_genetic_codes()


class _GeneticCodePrototype(int):

    def __new__(cls, value: int):
        obj = int.__new__(cls, value)  # type: ignore
        obj._value_ = value
        if value == 0:
            obj.text = 'Unknown'
            obj.stops = []
        else:
            obj.text = _GC_DESCRIPTIONS[value].name
            obj.stops = _GC_DESCRIPTIONS[value].stops()
        return obj


GeneticCode = Enum(  # type: ignore
    'GeneticCode',
    dict(
        **{"Unknown": 0},
        **{"SGC"+str(gc_id - 1): gc_id for gc_id in _GC_DESCRIPTIONS},
    ),
    type=_GeneticCodePrototype)


class NoReadingFrames(Exception):
    def __init__(self, gene_name: str):
        self.gene_name = gene_name
        super().__init__(
            f'No possible reading frames exist for gene {repr(gene_name)}')


class BadReadingFrame(Exception):
    def __init__(self, gene_name: str, reading_frame: ReadingFrame):
        self.gene_name = gene_name
        self.reading_frame = reading_frame
        super().__init__((
            f'Bad reading frame for gene {repr(gene_name)}: '
            f'multiple stop codons detected for {repr(reading_frame)}'))


class AmbiguousReadingFrame(Exception):
    def __init__(self, gene_name: str, reading_frame_set: Set[ReadingFrame]):
        self.gene_name = gene_name
        self.reading_frame_set = reading_frame_set
        super().__init__((
            f'Ambiguous reading frame for gene {repr(gene_name)}: '
            f'possible values: {repr(reading_frame_set)}'))


def collect_stop_codons() -> Dict[str, Set[int]]:
    stop_codons: Dict[str, Set[int]] = DefaultDict(set)
    for gc_table in GeneticCode:
        if gc_table == GeneticCode(0):
            continue
        for stop in gc_table.stops:  # type: ignore
            stop_codons[stop].add(gc_table._value_)
    return stop_codons


STOP_CODONS: Dict[str, Set[int]] = collect_stop_codons()
STOP_CODONS_SET = set(STOP_CODONS.keys())

TABLE_SET = set(
    gc_table._value_ for gc_table in GeneticCode if gc_table != GeneticCode.Unknown)


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
    seq_len = len(sequence)
    for stop_match in stops_regex.finditer(sequence, overlapped=True):
        if stop_match.start() < 3:
            continue
        frame = - ((seq_len - stop_match.end() - 1) % 3 + 1)
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
        frame: TABLE_SET.copy() for frame in [-3, -2, -1, 1, 2, 3]
    }
    for stop_codon, frame in stops:
        result[frame] -= stop_codons[stop_codon]
    return result


def choose_frame(frame_to_table: Dict[int, Set[int]]) -> List[int]:
    return [frame for frame, tables in frame_to_table.items() if tables]


def detect_reading_frame(sequence: str,
                         gc_table: GeneticCode = GeneticCode(0)) -> List[int]:
    if gc_table:
        stop_codons_set = set(gc_table.stops)  # type: ignore
    else:
        stop_codons_set = STOP_CODONS_SET
    stops = set(
        itertools.chain(
            detect_stop_codons(sequence, stop_codons_set),
            detect_reverse_stop_codons(sequence, stop_codons_set),
        )
    )
    if gc_table:
        disallowed_frames = {frame for frame, count in
                             Counter((frame for codon, frame in stops)).items()
                             if count == 1}
        return [frame for frame in [1, 2, 3, -1, -2, -3]
                if frame not in disallowed_frames]
    else:
        return choose_frame(allowed_translations(stops, STOP_CODONS))


def detect_reading_combinations(sequence: str,
                                gc_table: GeneticCode = GeneticCode(0)
                                ) -> Set[Tuple[GeneticCode, int]]:
    if gc_table:
        stop_codons_set = set(gc_table.stops)  # type: ignore
    else:
        stop_codons_set = STOP_CODONS_SET
    stops = set(
        itertools.chain(
            detect_stop_codons(sequence, stop_codons_set),
            detect_reverse_stop_codons(sequence, stop_codons_set),
        )
    )
    if gc_table:
        disallowed_frames = {frame for frame, count in
                             Counter((frame for codon, frame in stops)).items()
                             if count == 1}
        return {(gc_table, frame) for frame in [1, 2, 3, -1, -2, -3]
                if frame not in disallowed_frames}
    else:
        return {(GeneticCode(gc_id), frame) for frame, gcs in allowed_translations(
            stops, STOP_CODONS).items() for gc_id in gcs}


def last_codon_slice(seq: str, reading_frame: int) -> Optional[slice]:
    read_offset = abs(reading_frame) - 1
    if len(seq) < read_offset + 3:
        return None
    leftover = (len(seq) - read_offset) % 3
    stop = len(seq) - leftover
    start = stop - 3
    if reading_frame < 0:
        stop = - stop - 1
        start = - start - 1
    return slice(start, stop)


def last_codon(seq: str, reading_frame: int) -> Optional[str]:
    """
    Return last codon in `seq` according to `reading_frame`
    """
    the_slice = last_codon_slice(seq, reading_frame)
    if the_slice is None:
        return None
    return seq[the_slice]


ReadingCombination = Tuple[GeneticCode, int]


def filter_by_last_codons(column: pd.Series,
                          reading_combinations: Set[ReadingCombination]) -> None:
    """
    Removes elements of `reading_combinations`,
    which don't correspond to the last codon
    of at least one sequence from `column`
    """
    for _, seq in column.items():
        for gc, frame in reading_combinations.copy():
            if last_codon(seq, frame) not in gc.stops:  # type: ignore
                reading_combinations.remove((gc, frame))


def column_reading_combinations(
    column: pd.Series,
    gc_table: GeneticCode = GeneticCode(0)
) -> Set[ReadingCombination]:
    """
    Returns the set of reading frames that are valid for all sequences in `column`.

    """
    reading_combinations: Optional[Set[Tuple[GeneticCode, int]]] = None
    for _, seq in column.items():
        seq_reading_combinations = detect_reading_combinations(seq, gc_table)
        if reading_combinations is None:
            reading_combinations = seq_reading_combinations
        else:
            reading_combinations = reading_combinations.intersection(
                seq_reading_combinations)
    assert reading_combinations is not None
    return reading_combinations


def extract_frames(reading_combinations: Set[ReadingCombination]) -> Set[ReadingFrame]:
    return {ReadingFrame.from_int(frame) for _, frame in reading_combinations}


def final_column_reading_frame(
    column: pd.Series,
    genetic_code: GeneticCode = GeneticCode(0),
    reading_frame: ReadingFrame = ReadingFrame.from_int(0),
) -> ReadingFrame:
    """
    Determine the correct reading frame for a given Series,
    under optional assumptions for the genetic code and reading frame.
    """
    if column.empty:
        raise AmbiguousReadingFrame(column.name, set())
    reading_combinations = column_reading_combinations(column, genetic_code)
    possible_frames = extract_frames(reading_combinations)
    if not possible_frames:
        raise NoReadingFrames(column.name)
    elif not reading_frame:
        if len(possible_frames) > 1:
            # Could still determine a reading frame here
            # if a singular possible frame ends with a stop codon
            filter_by_last_codons(column, reading_combinations)
            possible_frames = extract_frames(reading_combinations)
            if not possible_frames:
                raise NoReadingFrames(column.name)
            elif len(possible_frames) > 1:
                raise AmbiguousReadingFrame(column.name, possible_frames)
        return ReadingFrame.from_int(possible_frames.pop())
    elif reading_frame:
        if reading_frame not in possible_frames:
            raise BadReadingFrame(column.name, reading_frame)
        return reading_frame
    else:
        assert False, "Unreachable 'else' branch in 'final_column_reading_frame'"


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
        return tuple([column.str[start::3] for start in starts])  # type: ignore
    else:
        return tuple([column.str[start::-3] for start in starts])  # type: ignore
