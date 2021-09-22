#!/usr/bin/env python3

from typing import TextIO

import sys
import os
import tempfile
import pathlib

from itaxotools.DNAconvert import convertDNA, parse_format

def _dnaconvert(dna_input: str, dna_output: str, format: str) -> None:
    tabformat = parse_format("tab", ("", ""))
    assert tabformat is not None
    outformat = parse_format(format, ("", ""))
    assert outformat is not None
    with open(dna_input) as dna_infile, open(dna_output, mode="w") as dna_outfile:
        convertDNA(dna_infile, dna_outfile, tabformat, outformat,
                   allow_empty_sequences=False, disable_automatic_renaming=False)

def process(infile: TextIO, outfile: TextIO, format: str) -> None:
    with tempfile.TemporaryDirectory(prefix='dnaconvert_') as dir:
        path = pathlib.Path(dir)
        dna_input = str(path / 'input')
        dna_output = str(path / 'output')

        with open(dna_input, mode="w") as dna_infile:
            for line in infile:
                dna_infile.write(line)

        _dnaconvert(dna_input, dna_output, format)

        with open(dna_output) as dna_outfile:
            for line in dna_outfile:
                outfile.write(line)
