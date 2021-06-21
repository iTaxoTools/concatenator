#!/usr/bin/env python3

from typing import BinaryIO, Iterator, List
import logging

import pandas as pd

import library.phylip as phylip
from library.utils import *
from library.file_utils import ZipOutput

CHARSET_START = """
## ALIGNMENT FILE ##
alignment = Cophyla21Jan2019noninterleaved.phy;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = mrbayes;

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = aicc;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]
"""

CHARSET_END = """
## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = greedy;
"""


def write(rows: Iterator[List[str]], output: BinaryIO) -> None:
    header = next(rows)
    row_table = pd.DataFrame(rows, columns=header)
    if len(row_table) < 1:
        raise ValueError("No data in input file")
    description_columns = [
        column for column in row_table.columns if not column.startswith("sequence")
    ]
    sequence_columns = [
        column for column in row_table.columns if column.startswith("sequence")
    ]

    charset = {}
    for column in sequence_columns:
        max_length = max_length_if_not_same(row_table[column])
        if max_length:
            logging.warning(
                f"Column '{column}' has sequences of unequal length.\n They will be padded with to the same length with 'N'"
            )
            row_table[column] = make_equal_length(row_table[column], max_length)
        charset[column] = len(row_table[column].iat[0])
    content = row_table[description_columns].copy()
    content["sequence"] = row_table[sequence_columns[0]].str.cat(
        row_table[sequence_columns[1:]]
    )

    archive = ZipOutput(output)
    with archive.open("partition_finder.phy") as outfile:
        phylip.write_column(content, "sequence", outfile)

    with archive.open("partition_finder.cfg") as charset_file:
        print(CHARSET_START, file=charset_file)
        position = 1
        for name, length in charset.items():
            print(
                name,
                " = ",
                position,
                "-",
                position + length - 1,
                ";",
                sep="",
                file=charset_file,
            )
            position += length
        print(CHARSET_END, file=charset_file)
