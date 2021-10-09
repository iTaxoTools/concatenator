
from typing import Callable, Dict, Iterator
from pathlib import Path

import pandas as pd

from .utils import removeprefix
from .file_types import FileType, FileFormat
from .detect_file_type import autodetect

from .fasta import fasta_writer
from .phylip import phylip_writer
from .ali import ali_writer


CallableWriter = Callable[[Iterator[pd.Series], Path], None]
CallableWriterDecorator = Callable[[CallableWriter], CallableWriter]

file_writers: Dict[FileType, Dict[FileFormat, CallableWriter]] = {
    type: dict() for type in FileType}


def file_writer(
    type: FileType, format: FileFormat
) -> CallableWriterDecorator:
    def decorator(func: CallableWriter) -> CallableWriter:
        file_writers[type][format] = func
        return func
    return decorator


@file_writer(FileType.File, FileFormat.Fasta)
def writeFastaFile(iterator: Iterator[pd.Series], path: Path) -> None:
    # Pending filters: remove empty
    joined = pd.concat(iterator, axis=1)
    data = joined.apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    with path.open('w') as file:
        fasta_writer(data, file)


@file_writer(FileType.File, FileFormat.Phylip)
def writePhylipFile(iterator: Iterator[pd.Series], path: Path) -> None:
    # Pending filters: remove empty, replace Nn- with ?, pad with ?
    joined = pd.concat(iterator, axis=1)
    data = joined.apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    with path.open('w') as file:
        phylip_writer(data, file)


@file_writer(FileType.File, FileFormat.Ali)
def writeAliFile(iterator: Iterator[pd.Series], path: Path) -> None:
    # Pending filters: remove empty, pad with N
    joined = pd.concat(iterator, axis=1)
    data = joined.apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    with path.open('w') as file:
        ali_writer(data, file)


class WriterNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__(f'No writer for {str(type)} and {str(format)}')


def format_file_name(
    base: str,
    type: FileType,
    format: FileFormat,
) -> str:
    if type.extension is None:
        return base + format.extension
    return base + type.extension


def write_from_iterator(
    series: Iterator[pd.Series],
    path: Path,
    type: FileType,
    format: FileFormat,
) -> None:
    """Species as index, sequences as columns"""
    if format not in file_writers[type]:
        raise WriterNotFound(type)
    return file_writers[type][format](series, path)
