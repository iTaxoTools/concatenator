#!/usr/bin/env python

from typing import TextIO, Iterator

import pandas as pd


def read(file: TextIO) -> Iterator[pd.DataFrame]:
    columns = file.readline().rstrip().split("\t")
    description_columns = [
        column for column in columns if not column.startswith("sequence_")
    ]
    sequence_columns = [column for column in columns if column.startswith("sequence_")]
    for sequence_column in sequence_columns:
        file.seek(0)
        columns = description_columns + [sequence_column]
        table = pd.read_table(file, usecols=columns, na_filter=False)[columns]
        yield table
