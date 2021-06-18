#!/usr/bin/env python3

from typing import TextIO, BinaryIO

import library.tabfile as tabfile
import library.multifile as multifile
import library.fasta as fasta
import library.phylip as phylip
import library.ali as ali


def process_fasta(infile: BinaryIO, outfile: TextIO) -> None:
    tabfile.write_from_columns(multifile.read_zip(fasta.column_reader, infile), outfile)


# def process_phylip(infile: BinaryIO, outfile: TextIO) -> None:
#     tabfile.write_from_columns(
#         multifile.read_zip(phylip.column_reader, infile), outfile
#     )


def process_ali(infile: BinaryIO, outfile: TextIO) -> None:
    tabfile.write_from_columns(multifile.read_zip(ali.column_reader, infile), outfile)
