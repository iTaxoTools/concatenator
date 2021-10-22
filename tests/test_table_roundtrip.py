#!/usr/bin/env python3

from pathlib import Path

import pytest

from itaxotools.concatenator import (
    FileType, FileFormat, write_to_path, read_from_path, get_writer)
from itaxotools.concatenator.library.file_writers import WriterNotFound
from .gen_stream import gen_table, stream_table


@pytest.mark.parametrize("filetype", list(FileType))
@pytest.mark.parametrize("format", list(FileFormat))
def test_table_roundtrip(
        filetype: FileType, format: FileFormat, tmp_path: Path) -> None:
    try:
        get_writer(filetype, format)
    except WriterNotFound:
        return
    table = gen_table(filetype, format)
    tmp_file = tmp_path / (str(filetype) + str(format))
    write_to_path(stream_table(table), tmp_file, filetype, format)
    for col_original, col_roundtrip in zip(
            stream_table(table), read_from_path(tmp_file)):
        assert (col_original == col_roundtrip).all()
