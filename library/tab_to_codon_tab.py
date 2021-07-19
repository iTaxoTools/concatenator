#!/usr/bin/env python3

from typing import TextIO, Iterator
import logging

import pandas as pd

import library.tabfile as tabfile
import library.codons as codons


def split_charsets(columns: Iterator[pd.Series]) -> Iterator[pd.Series]:
    for column in columns:
        reading_frames = codons.column_reading_frame(column)
        if not reading_frames:
            raise ValueError(f"{column.name} has no valid reading frame")
        if len(reading_frames) > 1:
            logging.warning(
                f"{column.name} has ambiguous reading frames: {reading_frames}"
            )
        frame = reading_frames[0]
        for i, new_column in enumerate(codons.split_codon_charsets(column, frame)):
            new_column.name = str(column.name) + "_" + str(i)
            yield new_column


def process(infile: TextIO, outfile: TextIO) -> None:
    columns = (gene_table.iat[:, -1] for gene_table in tabfile.read(infile))
    tabfile.write_from_columns(split_charsets(columns), outfile)
