#!/usr/bin/env python3

from library.multifile import ColumnWriter
import logging
from typing import TextIO

import pandas as pd

from library.utils import *


def write_column(column: pd.DataFrame, gene_name: str, outfile: TextIO) -> None:
    column["seqid"] = into_seqids(column.iloc[:, :-1].copy())
    max_length = max_length_if_not_same(column[gene_name])
    if max_length:
        logging.warning(
            f"Column '{gene_name}' has sequences of unequal length.\n They will be padded with to the same length with 'N'"
        )
        column[gene_name] = make_equal_length(column[gene_name], max_length)

    # write Phylip heading
    seq_length = len(column[gene_name].iat[0])
    print(len(column), seq_length, file=outfile)
    for seqid, sequence in column[["seqid", gene_name]].itertuples(
        index=False, name=None
    ):
        print(seqid, sequence, file=outfile)


column_writer = ColumnWriter(".phy", write_column)
