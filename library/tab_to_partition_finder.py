#!/usr/bin/env python3

from typing import BinaryIO, TextIO
import library.tabfile as tabfile
import library.partition_finder as partition_finder


def process(input: TextIO, output: BinaryIO) -> None:
    rows = tabfile.read_rows(input)
    partition_finder.write(rows, output)
