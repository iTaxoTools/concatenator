#!/usr/bin/env python3

from pathlib import Path
from typing import Dict
from dataclasses import dataclass

import pytest
import pandas as pd

from itaxotools.concatenator import (
    FileType, FileFormat, GeneSeries, GeneStream,
    autodetect, get_writer)
from itaxotools.concatenator.library.file_writers import (
    FileWriter)

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


write_tests = [
    WriteTest(FileType.File, FileFormat.Nexus, {}, 'stream_simple', 'simple.nex'),
    WriteTest(FileType.File, FileFormat.Nexus, {}, 'stream_altered', 'altered.nex'),
]


def assert_eq_files(type: FileType, file1: Path, file2: Path) -> None:
    if type == FileType.Directory:
        for (part1, part2) in zip(sorted(file1.iterdir()), sorted(file2.iterdir())):
            assert part1.name == part2.name
            assert part1.read_text() == part2.read_text()
    elif type == FileType.ZipArchive:
        for (part1, part2) in zip(sorted(ZipPath(file1).iterdir()), sorted(ZipPath(file2).iterdir())):
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
