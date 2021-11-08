#!/usr/bin/env python3

from typing import BinaryIO, Iterator, List
import logging

import pandas as pd

from . import phylip
from .utils import *
from .model import PathLike
from .types import Charset
from .file_utils import ZipOutput

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

CHARSET_FORMAT = """
## ALIGNMENT FILE ##
alignment = {alignment};

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = {branchlengths};

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = {models};

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = {model_selection};

## DATA BLOCKS: see manual for how to define ##
[data_blocks]

{data_blocks}

## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = {search};
"""


def write_rows(rows: Iterator[List[str]], output: BinaryIO) -> None:
    header = next(rows)
    row_table = pd.DataFrame(rows, columns=header)
    write_table(row_table, output)


def write_table(row_table: pd.DataFrame, output: BinaryIO) -> None:
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


def write_cfg(
    path: PathLike,
    charsets: List[Charset],
    alignment: str,
    branchlengths = 'linked',
    models = 'mrbayes',
    model_selection = 'aicc',
    search = 'greedy',
) -> None:

    data_blocks = list()
    for charset in charsets:
        if charset.frame:
            for codon in charset.codon_sets():
                data_blocks.append(f'{str(codon)}\n')
        else:
            data_blocks.append(f'{str(charset)}\n')

    with path.open('w', encoding='utf-8') as file:
        file.write(CHARSET_FORMAT.format(
            alignment = alignment,
            branchlengths = branchlengths,
            models = models,
            model_selection = model_selection,
            data_blocks = ''.join(data_blocks),
            search = search,
        ))
