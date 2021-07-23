#!/usr/bin/env python

from typing import TextIO

import library.tabfile as tabfile
import library.nexus as nexus

def process(input: TextIO, output: TextIO) -> None:
    nexus.write(tabfile.read(input), output)
