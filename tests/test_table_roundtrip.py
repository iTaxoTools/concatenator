#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import pytest

from itaxotools.concatenator import (
    FileType, FileFormat, get_extension,
    write_to_path, read_from_path, get_writer)
from itaxotools.concatenator.library.file_writers import WriterNotFound
from gen_stream import gen_table, stream_table


@pytest.mark.parametrize("filetype", list(FileType))
@pytest.mark.parametrize("format", list(FileFormat))
def test_table_roundtrip(
        filetype: FileType, format: FileFormat, tmp_path: Path) -> None:
    try:
        get_writer(filetype, format)
    except WriterNotFound:
        return
    table = gen_table(filetype, format)
    table.to_csv(tmp_path / 'generated.tsv', sep="\t", line_terminator="\n")
    tmp_file = tmp_path / f'output{get_extension(filetype, format)}'
    write_to_path(stream_table(table), tmp_file, filetype, format)
    columns_original = sorted(stream_table(table), key=lambda s: s.name)
    columns_roundtrip = sorted(read_from_path(tmp_file), key=lambda s: s.name)
    assert len(columns_original) == len(columns_roundtrip)
    for column_original, column_roundtrip in zip(columns_original, columns_roundtrip):
        assert column_original.name == column_roundtrip.name
        assert type(column_original.index) == type(column_roundtrip.index)
        # for (k1, v1), (k2, v2) in zip(column_original.iteritems(), column_roundtrip.iteritems()):
        #     assert k1 == k2
        #     assert v1 == v2
        assert column_original.equals(column_roundtrip)
