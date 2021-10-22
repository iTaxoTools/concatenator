#!/usr/bin/env python3

from pathlib import Path

import pytest

from itaxotools.concatenator import autodetect, read_from_path, write_to_path, FileType

TEST_DATA_DIR = Path("tests") / "data"  # I'm not sure if this is a good idea


def assert_eq_files(type: FileType, file1: Path, file2: Path) -> None:
    if type == FileType.Directory:
        for (part1, part2) in zip(file1.iterdir(), file2.iterdir()):
            assert part1.read_text() == part2.read_text()
    elif type == FileType.ZipArchive:
        for (part1, part2) in zip(file1.iterdir(), file2.iterdir()):
            assert part1.read_text() == part2.read_text()
    else:
        assert file1.read_text() == file2.read_text()


@pytest.mark.parametrize("data_dir", TEST_DATA_DIR.iterdir())
def test_read_write(data_dir: Path, tmp_path: Path) -> None:
    for input_file in data_dir.iterdir():
        for expected_output in data_dir.iterdir():
            expected_type, expected_format = autodetect(expected_output)
            test_path = tmp_path / f"test_{expected_type}_{expected_format}"
            write_to_path(read_from_path(input_file), test_path,
                          expected_type, expected_format)
            assert_eq_files(expected_type, test_path, expected_output)
