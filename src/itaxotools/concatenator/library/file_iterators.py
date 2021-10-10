
from typing import Callable, Dict, Iterator, TextIO
from pathlib import Path

import pandas as pd

from .utils import removeprefix
from .file_types import FileFormat, FileType
from .file_utils import iterateZipArchive, iterateDirectory
from .detect_file_type import autodetect
from .operators import check_valid, index_to_multi

from .nexus import read as nexus_read
from .ali import ali_reader
from .fasta import fasta_reader
from .phylip import phylip_reader


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
    data.index.name = 'species'
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
    series.index.name = 'species'
    return index_to_multi(series)


def readAliSeries(path: Path) -> pd.Series:
    return _readSeries(path, ali_reader)


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

    @file_iterator(FileType.Directory, format)
    def _iterateMultifileDir(path: Path) -> Iterator[pd.Series]:
        for part in iterateDirectory(path):
            yield reader(part)

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
        # index_columns = ['species', 'specimen-voucher', 'locality']
        index_columns = ['species']
        index = pd.read_table(
            file, usecols=index_columns, dtype=str, na_filter=False)
        for sequence in sequences:
            file.seek(0)
            table = pd.read_table(
                file, usecols=[sequence], dtype=str, na_filter=False)
            data = table.join(index)
            data.set_index(index_columns, inplace=True)
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
            yield check_valid(series)
        return
    raise IteratorNotFound(type, format)
