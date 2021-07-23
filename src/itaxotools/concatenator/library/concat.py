#!/usr/bin/env python

from typing import TextIO

from . import tabfile


def process(input: TextIO, output: TextIO) -> None:
    rows = tabfile.read_rows(input)
    concatenated_rows = tabfile.concat_sequences(rows)
    tabfile.write_rows(concatenated_rows, output)
