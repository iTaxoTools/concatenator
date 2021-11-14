#!/usr/bin/env python3

import logging
from typing import TextIO, Dict

import pandas as pd

from .model import GeneSeries, PathLike
from .multifile import ColumnWriter
from .utils import *


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


def column_reader(infile: TextIO) -> pd.Series:
    # skip the first line
    infile.readline()

    sequences: Dict[str, str] = {}

    for line in infile:
        # skip blank lines
        if line == "" or line.isspace():
            continue
        # separate name and sequence
        name, sequence = line.rstrip().rsplit(' ', 1)
        name = name.strip()
        # return the record
        sequences[name] = sequence
    series = pd.Series(sequences)
    series.index.name = 'seqid'
    return series


def gene_from_path(path: PathLike) -> GeneSeries:
    with path.open(encoding='utf-8') as file:
        series = column_reader(file)
    series.index.name = 'seqid'
    series.name = path.stem
    return GeneSeries(series, missing='Nn?', gap='-')


def phylip_writer(gene: GeneSeries, outfile: TextIO, relaxed: bool = True) -> None:
    series = gene.series
    assert has_uniform_length(series)
    seq_length = len(series.iat[0])
    outfile.write(f'{len(series)} {seq_length}\n')
    for index, sequence in series.iteritems():
        if isinstance(index, tuple):
            # this should never happen
            index = '_'.join([str(x) for x in index if x is not None])
        if not relaxed:
            assert len(index) > 10
        index = index.ljust(10)
        outfile.write(f'{index} {sequence}\n')


def gene_to_path(gene: GeneSeries, path: PathLike, relaxed: bool = True) -> None:
    with path.open('w', encoding='utf-8') as file:
        phylip_writer(gene, file, relaxed)
