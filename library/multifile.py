#!/usr/bin/env python3

from typing import Iterator, BinaryIO, TextIO, Callable

import pandas as pd

from library.file_utils import ZipOutput


def write_zip(
    write_column: Callable[[pd.DataFrame, str, TextIO], None],
    extension: str,
    columns: Iterator[pd.DataFrame],
    out_archive: BinaryIO,
) -> None:
    zip_out = ZipOutput(out_archive)
    for column in columns:
        gene_name = column.columns[-1]
        with zip_out.open(gene_name + extension) as outfile:
            write_column(column, gene_name, outfile)
