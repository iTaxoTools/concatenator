#!/usr/bin/env python3

import logging
from typing import TextIO, Iterator, List, NamedTuple

import pandas as pd

from .model import GeneSeries, PathLike
from .utils import *
from .multifile import ColumnWriter


class AliTuple(NamedTuple):
    species: str
    ali_tag: str
    sequence: str


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


column_writer = ColumnWriter(".ali", write_column)


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
        # skip the blank lines or comments
        if line == "" or line[0] == "#" or line.isspace():
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
    return pd.Series(
        {
            chunk[0][1:].replace("@", "_"): "".join(
                line.replace("*", "_") for line in chunk[1:]
            )
            for chunk in split_file(infile)
        }
    )


def chunk_reader(infile: TextIO) -> Iterator[AliTuple]:
    for chunk in split_file(infile):
        parts = chunk[0].split('@')
        yield AliTuple(
            species=parts[0][1:],
            ali_tag=''.join(parts[1:]),
            sequence=chunk[1])


def ali_reader(infile: TextIO) -> pd.Series:
    series = pd.Series({
        (x.species, x.ali_tag): x.sequence for x in chunk_reader(infile)})
    series.index.names = ['seqid', 'ali_tag']
    if series.index.to_frame()['ali_tag'].eq('').all():
        series.reset_index(level='ali_tag', drop=True, inplace=True)
    return series


def gene_from_path(path: PathLike) -> GeneSeries:
    with path.open(encoding='utf-8', errors='surrogateescape') as file:
        series = pd.Series({
            (x.species, x.ali_tag): x.sequence for x in chunk_reader(file)})
    series.index.names = ['seqid', 'ali_tag']
    if series.index.to_frame()['ali_tag'].eq('').all():
        series.reset_index(level='ali_tag', drop=True, inplace=True)
    series.name = path.stem
    return GeneSeries(series, missing='?', gap='*')


def ali_writer(gene: GeneSeries, outfile: TextIO) -> None:
    series = gene.series
    assert has_uniform_length(series)

    trans_dict = str.maketrans("Nn-", "??*")
    series = series.str.translate(trans_dict)

    pos_num = len(series.iat[0])
    otu_num = len(series)
    missing_count = series.str.count(r"\?").sum()
    missing_percent = missing_count / (pos_num * otu_num) * 100

    outfile.write(f'#Number of positions: {pos_num}\n')
    outfile.write(f'#Number of OTUs: {otu_num}\n')
    outfile.write(f'#Percent of ?: {missing_percent}\n')
    outfile.write('#\n')

    for index, sequence in series.items():
        if isinstance(index, tuple):
            index = '_'.join([str(x) for x in index if x is not None])
        outfile.write('>' + index + '\n')
        outfile.write(sequence + '\n')


def gene_to_path(gene: GeneSeries, path: PathLike) -> None:
    with path.open('w', encoding='utf-8', errors='surrogateescape') as file:
        ali_writer(gene, file)
