#!/usr/bin/env python3

from typing import BinaryIO, TextIO
from . import tabfile
from . import partition_finder


def process(input: TextIO, output: BinaryIO) -> None:
    rows = tabfile.read_rows(input)
    partition_finder.write_rows(rows, output)
