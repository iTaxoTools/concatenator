
from typing import Callable, Dict, Iterator, TextIO
from pathlib import Path

import pandas as pd

from .utils import removeprefix
from .file_types import FileFormat, FileType
from .file_utils import iterateZipArchive
from .detect_file_type import autodetect

from .nexus import read as nexus_read
from .ali import column_reader as ali_reader
from .fasta import column_reader as fasta_reader
from .phylip import column_reader as phylip_reader


CallableIterator = Callable[[Path], Iterator[pd.Series]]
CallableIteratorDecorator = Callable[[CallableIterator], CallableIterator]

file_iterators: Dict[FileType, Dict[FileFormat, CallableIterator]] = {
    type: dict() for type in FileType}


def file_iterator(
    type: FileType, format: FileFormat
) -> CallableIteratorDecorator:
    def decorator(func: CallableIterator) -> CallableIterator:
        file_iterators[type][format] = func
        return func
    return decorator


def readNexusFile(path: Path) -> pd.DataFrame:
    with path.open() as file:
        data = nexus_read(file)
    data.set_index('seqid', inplace=True)
    data.index.name = None
    return data


@file_iterator(FileType.File, FileFormat.Nexus)
def iterateNexus(path: Path) -> Iterator[pd.Series]:
    data = readNexusFile(path)
    for col in data:
        yield data[col]


def _readSeries(
    path: Path,
    func: Callable[[TextIO], pd.Series]
) -> pd.Series:
    with path.open() as file:
        series = func(file)
    series.name = path.stem
    return series


def readAliSeries(path: Path) -> pd.Series:
    series = _readSeries(path, ali_reader)
    return series.str.replace('_', '-', regex=False)  # to be filtered


def readFastaSeries(path: Path) -> pd.Series:
    return _readSeries(path, fasta_reader)


def readPhylipSeries(path: Path) -> pd.Series:
    return _readSeries(path, phylip_reader)


def _register_multifile_iterator(
    format: FileFormat,
    reader: Callable[[Path], pd.Series]
) -> None:

    @file_iterator(FileType.File, format)
    def _iterateSingleFile(path: Path) -> Iterator[pd.Series]:
        yield reader(path)

    @file_iterator(FileType.ZipArchive, format)
    def _iterateMultifileZip(path: Path) -> Iterator[pd.Series]:
        for part in iterateZipArchive(path):
            yield reader(part)


for format, reader in {
    FileFormat.Fasta: readFastaSeries,
    FileFormat.Phylip: readPhylipSeries,
    FileFormat.Ali: readAliSeries,
}.items():
    _register_multifile_iterator(format, reader)


@file_iterator(FileType.File, FileFormat.Tab)
def iterateTabFile(path: Path) -> Iterator[pd.Series]:
    with path.open() as file:
        columns = file.readline().rstrip().split("\t")
        sequences = [col for col in columns if col.startswith("sequence_")]
        file.seek(0)
        index = pd.read_table(
            file, usecols=['species'], dtype=str, na_filter=False).iloc[:, 0]
        for sequence in sequences:
            file.seek(0)
            table = pd.read_table(
                file, usecols=[sequence], dtype=str, na_filter=False)
            data = table.join(index)
            data.set_index('species', inplace=True)
            series = pd.Series(data.iloc[:, 0])
            series.name = removeprefix(sequence, 'sequence_')
            series.index.name = None
            yield series


class IteratorNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__((f'No iterator for {str(type)} and {str(format)}'))


def iterator_from_path(path: Path) -> Iterator[pd.Series]:
    """Species as index, sequence as name"""
    type, format = autodetect(path)
    if format in file_iterators[type]:
        for series in file_iterators[type][format](path):
            yield series
        return
    raise IteratorNotFound(type, format)
