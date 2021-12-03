#!/usr/bin/env python3

from typing import Iterator, BinaryIO, TextIO, Callable
import os

import pandas as pd

from .file_utils import ZipOutput, ZipInput


class ColumnWriter:
    """Type of the `column_writer` argument of `write_zip`"""

    def __init__(
        self, extension: str, write_column: Callable[[pd.DataFrame, str, TextIO], None]
    ) -> None:
        self.extension = extension
        self.write_column = write_column


def write_zip(
    column_writer: ColumnWriter,
    columns: Iterator[pd.DataFrame],
    out_archive: BinaryIO,
) -> None:
    zip_out = ZipOutput(out_archive)
    for column in columns:
        gene_name = column.columns[-1]
        with zip_out.open(gene_name + column_writer.extension) as outfile:
            column_writer.write_column(column, gene_name, outfile)


ColumnReader = Callable[[TextIO], pd.Series]


def read_zip(column_reader: ColumnReader, in_archive: BinaryIO) -> Iterator[pd.Series]:
    archive = ZipInput(in_archive)
    for filename, file in archive.files():
        column = column_reader(file)
        column.name, _ = os.path.splitext(filename)
        yield column
