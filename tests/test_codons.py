from typing import List, Tuple
from enum import IntEnum

import pandas as pd
import pytest

from itaxotools.concatenator.library.codons import GeneData, GeneticCode, ReadingFrame, BadReadingFrame


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
    result = gene_good_standard_1st.detect_reading_frame()
    assert result.reading_frame == 1


def test_good_standard_unknown(gene_good_standard_unknown):
    result = gene_good_standard_unknown.detect_reading_frame()
    assert result.reading_frame == 1


def test_good_unknown_unknown(gene_good_unknown_unknown):
    result = gene_good_unknown_unknown.detect_reading_frame()
    assert result.reading_frame == 1


def test_bad_standard_1st(gene_bad_standard_1st):
    with pytest.raises(BadReadingFrame):
        gene_bad_standard_1st.detect_reading_frame()


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
