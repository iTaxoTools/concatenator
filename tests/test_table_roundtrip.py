#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import pytest

from itaxotools.concatenator import (
    FileType, FileFormat, get_extension,
    write_to_path, read_from_path, get_writer)
from itaxotools.concatenator.library.file_writers import WriterNotFound
from gen_stream import gen_table, stream_table


def assert_series_equal(s1: pd.Series, s2: pd.Series):
    for (k1, v1), (k2, v2) in zip(s1.iteritems(), s2.iteritems()):
        assert k1 == k2
        assert v1 == v2


FORMATS = list(FileFormat)
FORMATS.remove(FileFormat.PartitionFinder)
FORMATS.remove(FileFormat.IQTree)

@pytest.mark.parametrize("filetype", list(FileType))
@pytest.mark.parametrize("format", FORMATS)
def test_table_roundtrip(
        filetype: FileType, format: FileFormat, tmp_path: Path) -> None:
    try:
        get_writer(filetype, format)
    except WriterNotFound:
        return
    table, name = gen_table(filetype, format)
    table.to_csv(tmp_path / 'generated.tsv', sep="\t", line_terminator="\n")
    tmp_file = tmp_path / f'{name}{get_extension(filetype, format)}'
    write_to_path(stream_table(table), tmp_file, filetype, format)
    genes_original = sorted(stream_table(table), key=lambda s: s.name)
    genes_roundtrip = sorted(read_from_path(tmp_file), key=lambda s: s.name)
    assert len(genes_original) == len(genes_roundtrip)
    for gene_original, gene_roundtrip in zip(genes_original, genes_roundtrip):
        assert gene_original.name == gene_roundtrip.name
        assert type(gene_original.series.index) == type(gene_roundtrip.series.index)
        assert_series_equal(gene_original.series, gene_roundtrip.series)
