#!/usr/bin/env python3

import logging
from typing import Iterator, BinaryIO, TextIO

import pandas as pd

from library.utils import *
from library.file_utils import ZipOutput


def write_column(column: pd.DataFrame, gene_name: str, outfile: TextIO) -> None:
    column["seqid"] = into_seqids(column.iloc[:, :-1].copy())

    # pad sequences if needed
    max_length = max_length_if_not_same(column[gene_name])
    if max_length:
        logging.warning(
            f"Column '{gene_name}' has sequences of unequal length.\n They will be padded with to the same length with '?'"
        )
        column[gene_name] = make_equal_length(
            column[gene_name], max_length, fillchar="?"
        )

    # do ali character translation
    trans_dict = str.maketrans("Nn-", "??*")
    column[gene_name] = column[gene_name].str.translate(trans_dict)

    # remove missing data
    column = column.loc[column[gene_name].str.contains(r"[^?*]", regex=True)]

    # write ali heading
    pos_num = len(column[gene_name].iat[0])
    otu_num = len(column)
    missing_count = column[gene_name].str.count(r"\?").sum()
    missing_percent = missing_count / (pos_num * otu_num) * 100

    print("#Number of positions:", pos_num, file=outfile)
    print("#Number of OTUs:", otu_num, file=outfile)
    print("#Percent of ?:", missing_percent, file=outfile)
    print("#", file=outfile)

    for seqid, sequence in column[["seqid", gene_name]].itertuples(
        index=False, name=None
    ):
        print(">", seqid, sep="", file=outfile)
        print(sequence, file=outfile)


def write_ali_zip(columns: Iterator[pd.DataFrame], out_archive: BinaryIO) -> None:
    zip_out = ZipOutput(out_archive)
    for column in columns:
        gene_name = column.columns[-1]
        with zip_out.open(gene_name + ".ali") as outfile:
            write_column(column, gene_name, outfile)
