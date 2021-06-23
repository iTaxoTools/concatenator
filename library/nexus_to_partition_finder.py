#!/usr/bin/env python3

from typing import TextIO, BinaryIO

import library.nexus as nexus
import library.partition_finder as partition_finder


def process(input: TextIO, output: BinaryIO) -> None:
    table = nexus.read(input)
    partition_finder.write_table(table, output)
