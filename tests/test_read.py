#!/usr/bin/env python3

from pathlib import Path
from typing import Dict
from dataclasses import dataclass

import pytest
import pandas as pd

from itaxotools.concatenator import (
    FileType, FileFormat, GeneSeries, GeneStream,
    autodetect, get_reader, get_loader, load_from_path)
from itaxotools.concatenator.library.operators import OpCheckValid

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


@dataclass
class File:
    name: str
    type: FileType
    format: FileFormat


@dataclass
class ReadTest:
    input: File
    kwargs: Dict
    missing: str
    gap: str


@pytest.fixture
def series_simple() -> pd.Series:
    series = pd.Series({
        'seqid1': 'GCAGTACCCTAA',
        'seqid2': 'GCAGTATAA',
        }, name='simple')
    return series


@pytest.fixture
def dataframe_multi() -> pd.DataFrame:
    df = pd.DataFrame({
        'simple': ['GCAGTACCCTAA', 'GCAGTATAA'],
        'sample': ['GTACCCTAA', 'GTATAA'],
        }, index = ['seqid1', 'seqid2'])
    return df


simple_tests = [
    ReadTest(File('simple.tsv', FileType.File, FileFormat.Tab), {}, 'Nn?', '-'),
    ReadTest(File('simple.nex', FileType.File, FileFormat.Nexus), {}, 'Nn?', '-'),
    ReadTest(File('simple.fas', FileType.File, FileFormat.Fasta), {}, 'Nn?', '-'),
    ReadTest(File('simple.phy', FileType.File, FileFormat.Phylip), {}, 'Nn?', '-'),
    ReadTest(File('simple.ali', FileType.File, FileFormat.Ali), {}, '?', '*'),
    ReadTest(File('with_spaces/simple.tsv', FileType.File, FileFormat.Tab), {}, 'Nn?', '-'),
    ReadTest(File('with_spaces/simple.nex', FileType.File, FileFormat.Nexus), {}, 'Nn?', '-'),
    ReadTest(File('with_spaces/simple.fas', FileType.File, FileFormat.Fasta), {}, 'Nn?', '-'),
    ReadTest(File('with_spaces/simple.phy', FileType.File, FileFormat.Phylip), {}, 'Nn?', '-'),
    ReadTest(File('with_spaces/simple.ali', FileType.File, FileFormat.Ali), {}, '?', '*'),
    ReadTest(File('altered.nex', FileType.File, FileFormat.Nexus), {}, '?', '*'),
]

multi_tests = [
    ReadTest(File('multi.tsv', FileType.File, FileFormat.Tab), {}, 'Nn?', '-'),
    ReadTest(File('multi.nex', FileType.File, FileFormat.Nexus), {}, 'Nn?', '-'),
    ReadTest(File('with_spaces/multi.tsv', FileType.File, FileFormat.Tab), {}, 'Nn?', '-'),
    ReadTest(File('with_spaces/multi.nex', FileType.File, FileFormat.Nexus), {}, 'Nn?', '-'),
]


def get_input_path(test: ReadTest) -> Path:
    return TEST_DATA_DIR / test.input.name


def get_stream(test: ReadTest) -> GeneStream:
    input_path = get_input_path(test)
    type, format = autodetect(input_path)
    assert type == test.input.type
    assert format == test.input.format
    reader = get_reader(test.input.type, test.input.format).update(**test.kwargs)
    return reader(input_path).pipe(OpCheckValid())


def assert_series_equal(s1: pd.Series, s2: pd.Series):
    for (k1, v1), (k2, v2) in zip(s1.items(), s2.items()):
        assert k1 == k2
        assert v1 == v2

def assert_matches(test: ReadTest, gene: GeneSeries, series: pd.Series) -> None:
    assert isinstance(gene, GeneSeries)
    assert_series_equal(gene.series, series)
    assert gene.name == series.name
    assert set(gene.missing) == set(test.missing)
    assert set(gene.gap) == set(test.gap)


@pytest.mark.parametrize("test", simple_tests)
def test_read_simple(test: ReadTest, series_simple: pd.Series) -> None:
    stream = get_stream(test)
    assert stream.source.type == test.input.type
    assert stream.source.format == test.input.format
    assert stream.source.path == get_input_path(test)
    genes = list(gene for gene in stream)
    assert len(genes) == 1
    for gene in genes:
        assert_matches(test, gene, series_simple)


@pytest.mark.parametrize("test", multi_tests)
def test_read_multi(test: ReadTest, dataframe_multi: pd.DataFrame) -> None:
    stream = get_stream(test)
    assert stream.source.type == test.input.type
    assert stream.source.format == test.input.format
    assert stream.source.path == get_input_path(test)
    genes = list(gene for gene in stream)
    assert len(genes) == 2
    for gene in genes:
        assert gene.name in dataframe_multi.columns
        assert_matches(test, gene, dataframe_multi[gene.name])


@pytest.mark.parametrize("test", simple_tests)
def test_load_simple(test: ReadTest, series_simple: pd.Series) -> None:
    input_path = TEST_DATA_DIR / test.input.name
    gdf = load_from_path(input_path)
    for col in gdf.dataframe:
        assert_series_equal(gdf.dataframe[col], series_simple)


@pytest.mark.parametrize("test", multi_tests)
def test_load_multi(test: ReadTest, dataframe_multi: pd.DataFrame) -> None:
    input_path = TEST_DATA_DIR / test.input.name
    gdf = load_from_path(input_path)
    df = gdf.dataframe
    assert len(df.columns) == 2
    for col in df:
        assert col in dataframe_multi.columns
        assert_series_equal(df[col], dataframe_multi[col])
