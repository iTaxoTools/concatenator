from typing import List, Tuple
from enum import IntEnum

import pandas as pd
import pytest

from itaxotools.concatenator.library.model import GeneSeries
from itaxotools.concatenator.library.codons import (
    GeneticCode, ReadingFrame, NoReadingFrames,
    BadReadingFrame, AmbiguousReadingFrame)
from itaxotools.concatenator.library.operators import OpDetectReadingFrame

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
def gene_good_standard_1st(series_good) -> GeneSeries:
    data = GeneSeries(series_good)
    data.genetic_code = GeneticCode(1)
    data.reading_frame = ReadingFrame(1)
    return data


@pytest.fixture
def gene_good_standard_unknown(series_good) -> GeneSeries:
    data = GeneSeries(series_good)
    data.genetic_code = GeneticCode(1)
    data.reading_frame = ReadingFrame(0)
    return data


@pytest.fixture
def gene_good_unknown_unknown(series_good) -> GeneSeries:
    data = GeneSeries(series_good)
    data.genetic_code = GeneticCode(0)
    data.reading_frame = ReadingFrame(0)
    return data


@pytest.fixture
def gene_bad_standard_1st(series_bad) -> GeneSeries:
    data = GeneSeries(series_bad)
    data.genetic_code = GeneticCode(1)
    data.reading_frame = ReadingFrame(1)
    return data


def test_good_standard_1st(gene_good_standard_1st):
    result = OpDetectReadingFrame()(gene_good_standard_1st)
    assert result.reading_frame == 1


def test_good_standard_unknown(gene_good_standard_unknown):
    result = OpDetectReadingFrame()(gene_good_standard_unknown)
    assert result.reading_frame == 1


def test_good_unknown_unknown(gene_good_unknown_unknown):
    result = OpDetectReadingFrame()(gene_good_unknown_unknown)
    assert result.reading_frame == 1


def test_bad_standard_1st(gene_bad_standard_1st):
    with pytest.raises(BadReadingFrame):
        OpDetectReadingFrame()(gene_bad_standard_1st)
