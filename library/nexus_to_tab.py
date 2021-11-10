#!/usr/bin/env python

from typing import TextIO

import library.nexus as nexus


def process(input: TextIO, output: TextIO) -> None:
    table = nexus.read(input)
    table.to_csv(output, sep="\t", index=False, line_terminator='\n')
