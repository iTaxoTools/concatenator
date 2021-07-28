#!/usr/bin/env python3

from typing import TextIO

import sys
import os
import tempfile

from itaxotools.DNAconvert import convertDNA, parse_format


def process(infile: TextIO, outfile: TextIO, format: str) -> None:
    _, dna_input = tempfile.mkstemp()
    _, dna_output = tempfile.mkstemp()
    with open(dna_input, mode="w") as dna_infile:
        for line in infile:
            dna_infile.write(line)
    tabformat = parse_format("tab", ("", ""))
    assert tabformat is not None
    outformat = parse_format(format, ("", ""))
    assert outformat is not None
    with open(dna_input) as dna_infile, open(dna_output, mode="w") as dna_outfile:
        convertDNA(dna_infile, dna_outfile, tabformat, outformat,
                   allow_empty_sequences=False, disable_automatic_renaming=False)

    with open(dna_output) as dna_outfile:
        for line in dna_outfile:
            outfile.write(line)
    os.unlink(dna_input)
    os.unlink(dna_output)
