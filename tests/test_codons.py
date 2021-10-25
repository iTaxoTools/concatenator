from typing import List, Tuple
from enum import IntEnum

import pandas as pd
import pytest


# Classes defined here for speculation, should be moved to library

class ReadingFrame(IntEnum):

    def __new__(cls, value: int, text: str):
        obj = int.__new__(cls, value)
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
    SGC3 = 'Mold/Protozoan/Coelenterate Mitochondrial; Mycoplasma; Spiroplasma', ['TAA', 'TAG']
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


@pytest.fixture
def series_good() -> pd.Series:
    series = pd.Series({
        'seqid1': 'GCAGTACCCTAA',
        'seqid2': 'GCAGTATAA',
        })
    series.name = 'gene1'
    return series


@pytest.fixture
def series_bad() -> pd.Series:
    series = pd.Series({
        'seqid1': 'GCAGTATAATAA',  # has two stop codons
        'seqid2': 'GCAGTATAA',
        })
    series.name = 'gene1'
    return series


@pytest.fixture
def gene_good_standard_1st(series_good) -> GeneData:
    data = GeneData(series_good)
    data.genetic_code = GeneticCode(1)
    data.reading_frame = ReadingFrame(1)
    return data


@pytest.fixture
def gene_good_standard_unknown(series_good) -> GeneData:
    data = GeneData(series_good)
    data.genetic_code = GeneticCode(1)
    data.reading_frame = ReadingFrame(0)
    return data


@pytest.fixture
def gene_good_unknown_unknown(series_good) -> GeneData:
    data = GeneData(series_good)
    data.genetic_code = GeneticCode(0)
    data.reading_frame = ReadingFrame(0)
    return data


@pytest.fixture
def gene_bad_standard_1st(series_bad) -> GeneData:
    data = GeneData(series_bad)
    data.genetic_code = GeneticCode(1)
    data.reading_frame = ReadingFrame(1)
    return data


def test_good_standard_1st(gene_good_standard_1st):
    result = some_func(gene_good_standard_1st)
    assert result.reading_frame == 1


def test_good_standard_unknown(gene_good_standard_unknown):
    result = some_func(gene_good_standard_unknown)
    assert result.reading_frame == 1


def test_good_unknown_unknown(gene_good_unknown_unknown):
    result = some_func(gene_good_unknown_unknown)
    assert result.reading_frame == 1


def test_bad_standard_1st(gene_bad_standard_1st):
    with pytest.raises(BadReadingFrame):
        some_func(gene_bad_standard_1st)


def test_model_reading_frame():
    assert not bool(ReadingFrame(0))
    assert ReadingFrame(1) == 1
    assert ReadingFrame(1).label == '+1'
    assert ReadingFrame(-2).label == '-2'


def test_model_genetic_code():
    assert not bool(GeneticCode(0))
    assert not GeneticCode(0).stops
    assert GeneticCode(1) == 1
    assert GeneticCode(1).stops
    assert GeneticCode(1).text == 'Standard'
    assert GeneticCode(0).text == 'Unknown'
