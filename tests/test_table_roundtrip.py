#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
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
    columns_original = sorted(stream_table(table), key=lambda s: s.name)
    columns_roundtrip = sorted(read_from_path(tmp_file), key=lambda s: s.name)
    assert len(columns_original) == len(columns_roundtrip)
    for column_original, column_roundtrip in zip(columns_original, columns_roundtrip):
        assert column_original.name == column_roundtrip.name
        assert column_original.equals(column_roundtrip)
