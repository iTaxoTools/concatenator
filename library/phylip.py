#!/usr/bin/env python3

import logging
from typing import Iterator, BinaryIO, TextIO

import pandas as pd

from library.utils import *
from library.file_utils import ZipOutput


def write_column(column: pd.DataFrame, gene_name: str, outfile: TextIO) -> None:
    column["seqid"] = into_seqids(column.iloc[:, :-1].copy())
    max_length = max_length_if_not_same(column[gene_name])
    if max_length:
        logging.warning(
            f"Column 'sequence_{gene_name}' has sequences of unequal length.\n They will be padded with to the same length with 'N'"
        )
        column[gene_name] = make_equal_length(column[gene_name], max_length)

    # write Phylip heading
    seq_length = len(column[gene_name].iat[0])
    print(len(column), seq_length, file=outfile)
    for seqid, sequence in column[["seqid", gene_name]].itertuples(
        index=False, name=None
    ):
        print(seqid, sequence, file=outfile)


def write_phylip_zip(columns: Iterator[pd.DataFrame], out_archive: BinaryIO) -> None:
    zip_out = ZipOutput(out_archive)
    for column in columns:
        gene_name = column.columns[-1]
        with zip_out.open(gene_name + ".phy") as outfile:
            write_column(column, gene_name, outfile)
