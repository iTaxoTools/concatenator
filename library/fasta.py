#!/usr/bin/env python3

from library.multifile import ColumnWriter
from typing import TextIO

import pandas as pd

from library.utils import *


def write_column(column: pd.DataFrame, gene_name: str, outfile: TextIO) -> None:
    column["seqid"] = into_seqids(column.iloc[:, :-1].copy())
    column = column.loc[column[gene_name].str.contains(r"[^-Nn?]", regex=True)]
    for seqid, sequence in column[["seqid", gene_name]].itertuples(
        index=False, name=None
    ):
        print(">", seqid, sep="", file=outfile)
        print(sequence, file=outfile)


column_writer = ColumnWriter(".fas", write_column)
