#!/usr/bin/env python3

from pathlib import Path

import pytest

from itaxotools.concatenator import autodetect, read_from_path, write_to_path, FileType
from itaxotools.concatenator.library.file_utils import ZipPath

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


def assert_eq_files(type: FileType, file1: Path, file2: Path) -> None:
    if type == FileType.Directory:
        for (part1, part2) in zip(sorted(file1.iterdir()), sorted(file2.iterdir())):
            assert part1.read_text() == part2.read_text()
    elif type == FileType.ZipArchive:
        for (part1, part2) in zip(sorted(ZipPath(file1).iterdir()), sorted(ZipPath(file2).iterdir())):
            assert part1.read_text() == part2.read_text()
    else:
        assert file1.read_text() == file2.read_text()


@pytest.mark.parametrize("expected_output", TEST_DATA_DIR.iterdir())
def test_read_write(expected_output: Path, tmp_path: Path) -> None:
    expected_type, expected_format = autodetect(expected_output)
    test_path = tmp_path / f"test_{expected_type}_{expected_format}"
    write_to_path(read_from_path(expected_output), test_path,
                  expected_type, expected_format)
    assert_eq_files(expected_type, test_path, expected_output)
