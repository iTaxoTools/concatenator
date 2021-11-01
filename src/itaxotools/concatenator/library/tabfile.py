#!/usr/bin/env python

from typing import TextIO, Iterator, TypeVar, List

import pandas as pd

from .model import GeneStream, GeneDataFrame, PathLike
from .utils import removeprefix


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


def write_from_columns(columns: Iterator[pd.Series], outfile: TextIO) -> None:
    table = pd.DataFrame()
    for column in columns:
        table = table.join(column, how="outer")
    table.index.name = "seqid"
    table.to_csv(outfile, sep="\t", line_terminator="\n")


T = TypeVar("T")


def pop_many(l: List[T], indices: List[int]) -> List[T]:
    """
    Pops the elements of `l` at `indices` and returns the list of them.

    Sorts `indices` in reverse order
    """
    indices.sort()
    indices.reverse()
    result = [l.pop(i) for i in indices]
    result.reverse()
    return result


def read_rows(input: TextIO) -> Iterator[List[str]]:
    for line in input:
        yield line.rstrip().split("\t")


def concat_sequences(rows: Iterator[List[str]]) -> Iterator[List[str]]:
    # read the header
    try:
        header = next(rows)
    except StopIteration:
        return

    # find the columns with description (i.e. not sequences) values
    description_indices = [
        i for i, column in enumerate(header) if "sequence" not in column.lower()
    ]
    description_columns = pop_many(header, description_indices)

    yield description_columns + ["sequence"]

    for row in rows:
        # separate description and sequences
        description_values = pop_many(row, description_indices)
        # yield description and the sequence
        yield description_values + [("".join(row))]


def write_rows(rows: Iterator[List[str]], output: TextIO) -> None:
    for row in rows:
        print(*row, sep="\t", file=output)


def dataframe_from_path(
    path: PathLike,
    sequence_prefix: str = 'sequence_'
) -> GeneDataFrame:
    data = pd.read_csv(path, sep='\t', dtype=str)
    indices = [x for x in data.columns if not x.startswith(sequence_prefix)]
    data.set_index(indices, inplace=True)
    data.columns = [
        removeprefix(col, sequence_prefix) for col in data.columns]
    gdf = GeneDataFrame(data, missing='?N', gap='-')
    return gdf

def stream_from_path(
    path: PathLike,
    sequence_prefix: str = 'sequence_'
) -> GeneStream:
    gdf = dataframe_from_path(path, sequence_prefix)
    return gdf.stream()
