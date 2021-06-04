#!/usr/bin/env python3

from typing import Iterator, BinaryIO, TextIO

import pandas as pd

from library.utils import *
from library.file_utils import ZipOutput


def write_column(column: pd.DataFrame, gene_name: str, outfile: TextIO) -> None:
    column["seqid"] = into_seqids(column)
    for seqid, sequence in column[["seqid", gene_name]].itertuples(
        index=False, name=None
    ):
        print(">", seqid, sep="", file=outfile)
        print(sequence, file=outfile)


def write_fasta_zip(columns: Iterator[pd.DataFrame], out_archive: BinaryIO) -> None:
    zip_out = ZipOutput(out_archive)
    for column in columns:
        gene_name = column.columns[-1]
        with zip_out.open(gene_name + ".fas") as outfile:
            write_column(column, gene_name, outfile)
