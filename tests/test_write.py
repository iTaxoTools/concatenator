#!/usr/bin/env python3

from pathlib import Path
from typing import Dict
from dataclasses import dataclass
from itertools import product

import pytest
import pandas as pd

from itaxotools.concatenator import (
    FileType, FileFormat, GeneSeries, GeneStream, get_writer)
from itaxotools.concatenator.library.file_writers import file_writers
from itaxotools.concatenator.library.codons import ReadingFrame
from itaxotools.concatenator.library.file_utils import ZipPath


TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


@dataclass
class WriteTest:
    type: FileType
    format: FileFormat
    kwargs: Dict
    stream: GeneStream
    file: str


@pytest.fixture
def stream_simple() -> GeneStream:
    series = pd.Series({
            'seq1': 'gc-Nn?TAA',
        }, name='gene')
    series.index.name = 'seqid'
    gene = GeneSeries(series)
    return GeneStream(iter([gene]))


@pytest.fixture
def stream_altered() -> GeneStream:
    series = pd.Series({
            'seq1': 'gc*???TAA',
        }, name='gene')
    series.index.name = 'seqid'
    gene = GeneSeries(series, missing='?', gap='*')
    return GeneStream(iter([gene]))


@pytest.fixture
def stream_reading_frames() -> GeneStream:
    series0 = pd.Series({'seq1': 'GCTA'}, name='gene0')
    series0.index.name = 'seqid'
    gene0 = GeneSeries(series0, reading_frame=ReadingFrame(1))

    series1 = pd.Series({'seq1': 'GCTA'}, name='gene1')
    series1.index.name = 'seqid'
    gene1 = GeneSeries(
        series1, reading_frame=ReadingFrame(1),
        codon_names=('**_A', '**_B', 'genus1_C'))

    series2 = pd.Series({'seq1': 'GCTA'}, name='gene2')
    series2.index.name = 'seqid'
    gene2 = GeneSeries(series2, reading_frame=ReadingFrame(2))

    series3 = pd.Series({'seq1': 'GCTA'}, name='gene3')
    series3.index.name = 'seqid'
    gene3 = GeneSeries(series3, reading_frame=ReadingFrame(3))

    series4 = pd.Series({'seq1': 'GCTA'}, name='gene4')
    series4.index.name = 'seqid'
    gene4 = GeneSeries(series4, reading_frame=ReadingFrame(-1))

    series5 = pd.Series({'seq1': 'GCTA'}, name='gene5')
    series5.index.name = 'seqid'
    gene5 = GeneSeries(series5, reading_frame=ReadingFrame(-2))

    series6 = pd.Series({'seq1': 'GCTA'}, name='gene6')
    series6.index.name = 'seqid'
    gene6 = GeneSeries(series6, reading_frame=ReadingFrame(-3))

    series7 = pd.Series({'seq1': 'GCTA'}, name='gene7')
    series7.index.name = 'seqid'
    gene7 = GeneSeries(series7, reading_frame=ReadingFrame(0))

    return GeneStream(iter([gene0, gene1, gene2, gene3, gene4, gene5, gene6, gene7]))


write_tests = [
    WriteTest(FileType.File, FileFormat.Nexus, {}, 'stream_simple', 'simple.nex'),
    WriteTest(FileType.File, FileFormat.Nexus, {}, 'stream_altered', 'altered.nex'),
    WriteTest(FileType.File, FileFormat.Nexus, dict(adjust_frames=False, include_full_markers_with_codons=False),
        'stream_reading_frames', 'reading_frames.nex'),
    WriteTest(FileType.File, FileFormat.Nexus, dict(adjust_frames=False, include_full_markers_with_codons=True),
        'stream_reading_frames', 'reading_frames_full.nex'),
    WriteTest(FileType.Directory, FileFormat.PartitionFinder, {}, 'stream_simple', 'partition_finder_simple'),
    WriteTest(FileType.Directory, FileFormat.IQTree, {}, 'stream_simple', 'iqtree_simple'),
    WriteTest(FileType.Directory, FileFormat.IQTree, dict(adjust_frames=False, include_full_markers_with_codons=False),
        'stream_reading_frames', 'iqtree_reading_frames'),
    WriteTest(FileType.Directory, FileFormat.IQTree, dict(adjust_frames=False, include_full_markers_with_codons=True),
        'stream_reading_frames', 'iqtree_reading_frames_full'),
]


def assert_eq_files(type: FileType, file1: Path, file2: Path) -> None:
    if type == FileType.Directory:
        files1 = sorted(file1.iterdir())
        files2 = sorted(file2.iterdir())
        assert len(files1) == len(files2)
        for (part1, part2) in zip(files1, files2):
            assert part1.name == part2.name
            assert part1.read_text() == part2.read_text()
    elif type == FileType.ZipArchive:
        files1 = sorted(ZipPath(file1).iterdir())
        files2 = sorted(ZipPath(file2).iterdir())
        assert len(files1) == len(files2)
        for (part1, part2) in zip(files1, files2):
            assert part1.name == part2.name
            assert part1.read_text() == part2.read_text()
    else:
        assert file1.read_text() == file2.read_text()


@pytest.mark.parametrize("test", write_tests)
def test_write(test: WriteTest, tmp_path: Path, request) -> None:
    stream = request.getfixturevalue(test.stream)
    output_path = TEST_DATA_DIR / test.file
    test_path = tmp_path / test.file
    writer = get_writer(test.type, test.format).update(**test.kwargs)
    writer(stream, test_path)
    assert_eq_files(test.type, test_path, output_path)


valid_type_formats = [
    (type, format) for type, format in product(list(FileType), list(FileFormat))
    if format in file_writers[type]
]


class SimpleFilterTest:
    def __init__(self):
        self.filtered = False

    def filter(self, stream: GeneStream) -> GeneStream:
        self.filtered = True
        return stream


@pytest.mark.parametrize("type, format", valid_type_formats)
def test_write_filter(
    type: FileType,
    format: FileFormat,
    stream_simple: GeneStream,
    tmp_path: Path
) -> None:

    simple_filter_test = SimpleFilterTest()
    out_path = tmp_path / 'out'
    writer = get_writer(type, format)
    writer.filters.append(simple_filter_test.filter)
    writer(stream_simple, out_path)
    assert simple_filter_test.filtered
    assert out_path.exists()
