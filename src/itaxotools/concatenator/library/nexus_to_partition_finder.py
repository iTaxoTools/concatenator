#!/usr/bin/env python3

from typing import TextIO, BinaryIO

from . import nexus
from . import partition_finder


def process(input: TextIO, output: BinaryIO) -> None:
    table = nexus.read(input)
    partition_finder.write_table(table, output)
