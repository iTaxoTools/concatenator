#!/usr/bin/env python

from typing import TextIO

from . import tabfile
from . import nexus


def process(input: TextIO, output: TextIO) -> None:
    nexus.write(tabfile.read(input), output)
