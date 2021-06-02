#!/usr/bin/env python3

import io
import zipfile
from typing import TextIO, Callable, BinaryIO


def make_binary(
    process: Callable[[TextIO, TextIO], None]
) -> Callable[[BinaryIO, BinaryIO], None]:
    def process_bin(infile: BinaryIO, outfile: BinaryIO) -> None:
        infile_t = io.TextIOWrapper(infile, errors="replace")
        outfile_t = io.TextIOWrapper(outfile)
        process(infile_t, outfile_t)

    return process_bin


def make_binary_in(
    process: Callable[[TextIO, BinaryIO], None]
) -> Callable[[BinaryIO, BinaryIO], None]:
    def process_bin(infile: BinaryIO, outfile: BinaryIO) -> None:
        infile_t = io.TextIOWrapper(infile, errors="replace")
        process(infile_t, outfile)

    return process_bin


def make_binary_out(
    process: Callable[[BinaryIO, TextIO], None]
) -> Callable[[BinaryIO, BinaryIO], None]:
    def process_bin(infile: BinaryIO, outfile: BinaryIO) -> None:
        outfile_t = io.TextIOWrapper(outfile)
        process(infile, outfile_t)

    return process_bin


class ZipOutput:
    """Class for writing text output files to a zip archive"""

    def __init__(self, file: BinaryIO):
        self.archive = zipfile.ZipFile(file, mode="w")

    def open(self, name: str) -> TextIO:
        """Opens a file in the archive for writing"""
        file = self.archive.open(name, mode="w")
        return io.TextIOWrapper(file)
