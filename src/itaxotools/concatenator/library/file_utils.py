#!/usr/bin/env python3

import io
import zipfile
from pathlib import Path
from zipfile import ZipFile
from zipp import Path as ZipPath_
from typing import TextIO, Callable, BinaryIO, Iterator, Tuple, Union

__all__ = ['ZipFile', 'ZipPath', ...]


class ZipPath(ZipPath_):
    # zipp.Path does not subclass pathlib.Path, so we must implement these.
    # When we stop supporting Python3.8 for Win7, these will have to go.
    def __eq__(self, other):
        return (self.root == other.root) and (self.at == other.at)

    def __lt__(self, other):
        return self.at < other.at


PathLike = Union[Path, ZipPath]


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


class ZipInput:
    """Class for reading input files from a zip archive"""

    def __init__(self, file: BinaryIO):
        self.archive = zipfile.ZipFile(file, mode="r")

    def files(self) -> Iterator[Tuple[str, TextIO]]:
        """Yields names and the files in the archive"""
        for filename in self.archive.namelist():
            with self.archive.open(filename, mode="r") as file:
                yield (filename, io.TextIOWrapper(file, errors="replace"))
