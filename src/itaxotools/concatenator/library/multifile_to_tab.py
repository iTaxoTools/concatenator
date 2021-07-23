#!/usr/bin/env python3

from typing import TextIO, BinaryIO

from . import tabfile
from . import multifile
from . import fasta
from . import phylip
from . import ali


def process_fasta(infile: BinaryIO, outfile: TextIO) -> None:
    tabfile.write_from_columns(multifile.read_zip(fasta.column_reader, infile), outfile)


def process_phylip(infile: BinaryIO, outfile: TextIO) -> None:
    tabfile.write_from_columns(
        multifile.read_zip(phylip.column_reader, infile), outfile
    )


def process_ali(infile: BinaryIO, outfile: TextIO) -> None:
    tabfile.write_from_columns(multifile.read_zip(ali.column_reader, infile), outfile)
