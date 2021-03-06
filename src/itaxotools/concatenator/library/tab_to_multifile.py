#!/usr/bin/env python3

from typing import TextIO, BinaryIO

from . import tabfile
from . import multifile
from . import fasta
from . import phylip
from . import ali


def process_fasta(infile: TextIO, outfile: BinaryIO) -> None:
    columns = tabfile.read(infile)
    multifile.write_zip(fasta.column_writer, columns, outfile)


def process_phylip(infile: TextIO, outfile: BinaryIO) -> None:
    columns = tabfile.read(infile)
    multifile.write_zip(phylip.column_writer, columns, outfile)


def process_ali(infile: TextIO, outfile: BinaryIO) -> None:
    columns = tabfile.read(infile)
    multifile.write_zip(ali.column_writer, columns, outfile)
