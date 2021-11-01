#!/usr/bin/env python3

from typing import TextIO, Iterator, List

import pandas as pd

from .model import GeneSeries, PathLike
from .file_types import FileType
from .utils import *
from .multifile import ColumnWriter


def write_column(column: pd.DataFrame, gene_name: str, outfile: TextIO) -> None:
    column["seqid"] = into_seqids(column.iloc[:, :-1].copy())
    column = column.loc[column[gene_name].str.contains(r"[^-Nn?]", regex=True)]
    for seqid, sequence in column[["seqid", gene_name]].itertuples(
        index=False, name=None
    ):
        print(">", seqid, sep="", file=outfile)
        print(sequence, file=outfile)


column_writer = ColumnWriter(".fas", write_column)


def split_file(file: TextIO) -> Iterator[List[str]]:
    """
    Returns iterator that yield records as lists of lines
    """
    # find the beginning of the first record
    line = " "
    while line[0] != ">":
        line = file.readline()

    # chunk contains the already read lines of the current record
    chunk = []
    # put the first line of the first record into chunk
    chunk.append(line.rstrip())

    for line in file:
        # skip the blank lines
        if line == "" or line.isspace():
            continue

        # yield the chunk if the new record has begun
        if line[0] == ">":
            yield chunk
            chunk = []

        # put the first line of the new record into chunk
        chunk.append(line.rstrip())

    # yield the last record
    yield chunk


def column_reader(infile: TextIO) -> pd.Series:
    series = pd.Series({
        chunk[0][1:]: "".join(chunk[1:])
        for chunk in split_file(infile)})
    series.index.name = 'seqid'
    return series


def gene_from_path(path: PathLike) -> GeneSeries:
    with path.open() as file:
        series = pd.Series({
            chunk[0][1:]: ''.join(chunk[1:])
            for chunk in split_file(file)})
    series.index.name = 'seqid'
    series.name = path.stem
    return GeneSeries(series, missing='?N', gap='-')


def fasta_writer(gene: GeneSeries, outfile: TextIO) -> None:
    for index, sequence in gene.series.iteritems():
        if isinstance(index, tuple):
            # this should never happen
            index = '_'.join([str(x) for x in index if x is not None])
        outfile.write('>' + index + '\n')
        outfile.write(sequence + '\n')
