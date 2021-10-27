#!/usr/bin/env python3

from pathlib import Path
from typing import Dict

import pytest

from itaxotools.concatenator import (
    autodetect, read_from_path, write_to_path, FileType, FileFormat,
    get_reader, get_writer)
from itaxotools.concatenator.library.file_utils import ZipPath

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


self_test_data = [
    ('sequences_ali', FileType.Directory, FileFormat.Ali, {}, {}),
    ('sequences_ali.zip', FileType.ZipArchive, FileFormat.Ali, {}, {}),
]

test_data = self_test_data


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


@pytest.mark.parametrize("file,type,format,reader_kwds,writer_kwds", test_data)
def test_read_write(
    file: Path,
    type: FileType,
    format: FileFormat,
    reader_kwds: Dict,
    writer_kwds: Dict,
    tmp_path: Path
) -> None:
    file_path = TEST_DATA_DIR / file
    test_path = tmp_path / file
    reader = get_reader(type, format)(*reader_kwds)
    writer = get_writer(type, format)(*writer_kwds)
    writer(reader(file_path), test_path)
    assert_eq_files(type, test_path, file_path)
