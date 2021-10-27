#!/usr/bin/env python3

from pathlib import Path
from typing import Dict
from dataclasses import dataclass

import pytest

from itaxotools.concatenator import (
    FileType, FileFormat, get_reader, get_writer)
from itaxotools.concatenator.library.file_utils import ZipPath
from itaxotools.concatenator.library.utils import Justification

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


@dataclass
class File:
    name: str
    type: FileType
    format: FileFormat


@dataclass
class SelfTest:
    file: File
    reader_kwds: Dict
    writer_kwds: Dict


@dataclass
class CrossTest:
    input: File
    output: File
    reader_kwds: Dict
    writer_kwds: Dict


self_test_data = [
    SelfTest(File('sequences_ali', FileType.Directory, FileFormat.Ali), {}, {}),
    SelfTest(File('sequences_ali.zip', FileType.ZipArchive, FileFormat.Ali), {}, {}),
    SelfTest(File('sequences.phy', FileType.File, FileFormat.Phylip), {}, {}),
    SelfTest(File('sequences_pad.fas', FileType.File, FileFormat.Fasta), {}, {}),
    SelfTest(File('sequences_no_pad.fas', FileType.File, FileFormat.Fasta), {}, dict(padding='')),
    SelfTest(File('sequences.tab', FileType.File, FileFormat.Tab), {}, {}),
    SelfTest(File('sequences_left_space.nex', FileType.File, FileFormat.Nexus), {},
        dict(justification=Justification.Left, separator=' ')),
    SelfTest(File('sequences_right_tab.nex', FileType.File, FileFormat.Nexus), {},
        dict(justification=Justification.Right, separator='\t')),
]

test_data = [
    CrossTest(test.file, test.file, test.reader_kwds, test.writer_kwds)
    for test in self_test_data
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


@pytest.mark.parametrize("test", test_data)
def test_read_write(test: CrossTest, tmp_path: Path) -> None:
    input_path = TEST_DATA_DIR / test.input.name
    output_path = TEST_DATA_DIR / test.output.name
    test_path = tmp_path / test.output.name
    reader = get_reader(test.input.type, test.input.format)(**test.reader_kwds)
    writer = get_writer(test.output.type, test.output.format)(**test.writer_kwds)
    writer(reader(input_path), test_path)
    assert_eq_files(test.output.type, test_path, output_path)
