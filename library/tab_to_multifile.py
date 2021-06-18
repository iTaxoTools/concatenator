#!/usr/bin/env python3

from typing import TextIO, BinaryIO

import library.tabfile as tabfile
import library.multifile as multifile
import library.fasta as fasta
import library.phylip as phylip
import library.ali as ali


def process_fasta(infile: TextIO, outfile: BinaryIO) -> None:
    columns = tabfile.read(infile)
    multifile.write_zip(fasta.write_column, ".fas", columns, outfile)


def process_phylip(infile: TextIO, outfile: BinaryIO) -> None:
    columns = tabfile.read(infile)
    multifile.write_zip(phylip.write_column, ".phy", columns, outfile)


def process_ali(infile: TextIO, outfile: BinaryIO) -> None:
    columns = tabfile.read(infile)
    multifile.write_zip(ali.write_column, ".ali", columns, outfile)
