#!/usr/bin/env python3

from typing import TextIO, BinaryIO

import library.tabfile as tabfile
import library.fasta as fasta
import library.phylip as phylip


def process_fasta(infile: TextIO, outfile: BinaryIO) -> None:
    columns = tabfile.read(infile)
    fasta.write_fasta_zip(columns, outfile)


def process_phylip(infile: TextIO, outfile: BinaryIO) -> None:
    columns = tabfile.read(infile)
    phylip.write_phylip_zip(columns, outfile)
