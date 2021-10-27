
from typing import Callable, Dict, Iterator, TextIO
from pathlib import Path

import pandas as pd

from .utils import ConfigurableCallable, Param, removeprefix
from .file_types import FileFormat, FileType
from .file_utils import ZipPath
from .file_identify import autodetect
from .operators import OpCheckValid, OpIndexToMulti

from .nexus import read as nexus_read
from .ali import ali_reader
from .fasta import fasta_reader
from .phylip import phylip_reader

from . import SPECIES, SEQUENCE_PREFIX


class FileReader(ConfigurableCallable):
    type: FileType = None
    format: FileFormat = None

    def call(self, path: Path) -> Iterator[pd.Series]:
        raise NotImplementedError


file_readers: Dict[FileType, Dict[FileFormat, FileReader]] = {
    type: dict() for type in FileType}


def file_reader(
    type: FileType, format: FileFormat
) -> Callable[[FileReader], FileReader]:
    def decorator(reader: FileReader) -> FileReader:
        file_readers[type][format] = reader
        reader.type = type
        reader.format = format
        return callable
    return decorator


def readNexusFile(path: Path) -> pd.DataFrame:
    with path.open() as file:
        data = nexus_read(file, sequence_prefix='')
    data.set_index('seqid', inplace=True)
    data.index.name = SPECIES
    return OpIndexToMulti()(data)


@file_reader(FileType.File, FileFormat.Nexus)
class NexusReader(FileReader):
    def call(self, path: Path) -> Iterator[pd.Series]:
        data = readNexusFile(path)
        for col in data:
            yield data[col]


def _readSeries(
    path: Path,
    part_reader: Callable[[TextIO], pd.Series]
) -> pd.Series:
    with path.open() as file:
        series = part_reader(file)
    series.name = path.stem
    series.index.name = SPECIES
    return OpIndexToMulti()(series)


def readAliSeries(path: Path) -> pd.Series:
    return _readSeries(path, ali_reader)


def readFastaSeries(path: Path) -> pd.Series:
    return _readSeries(path, fasta_reader)


def readPhylipSeries(path: Path) -> pd.Series:
    return _readSeries(path, phylip_reader)


def _register_multifile_reader(
    format: FileFormat,
    reader: Callable[[Path], pd.Series]
) -> None:

    @file_reader(FileType.File, format)
    class _SingleFileReader(FileReader):
        def call(self, path: Path) -> Iterator[pd.Series]:
            yield reader(path)

    @file_reader(FileType.Directory, format)
    class _MultiDirReader(FileReader):
        def call(self, path: Path) -> Iterator[pd.Series]:
            for part in path.iterdir():
                yield reader(part)

    @file_reader(FileType.ZipArchive, format)
    class _MultiZipReader(FileReader):
        def call(self, path: Path) -> Iterator[pd.Series]:
            for part in ZipPath(path).iterdir():
                yield reader(part)


for format, reader in {
    FileFormat.Fasta: readFastaSeries,
    FileFormat.Phylip: readPhylipSeries,
    FileFormat.Ali: readAliSeries,
}.items():
    _register_multifile_reader(format, reader)


@file_reader(FileType.File, FileFormat.Tab)
class TabFileReaderSlow(FileReader):
    def call(self, path: Path) -> Iterator[pd.Series]:
        with path.open() as file:
            columns = file.readline().rstrip().split("\t")
            sequences = [x for x in columns if x.startswith(SEQUENCE_PREFIX)]
            indices = [x for x in columns if not x.startswith(SEQUENCE_PREFIX)]
            file.seek(0)
            index = pd.read_table(
                file, usecols=indices, dtype=str)
            for sequence in sequences:
                file.seek(0)
                table = pd.read_table(
                    file, usecols=[sequence], dtype=str)
                data = table.join(index)
                data.set_index(indices, inplace=True)
                series = pd.Series(data.iloc[:, 0])
                series.name = removeprefix(sequence, SEQUENCE_PREFIX)
                yield OpIndexToMulti()(series)


def readTab(path: Path, sequence_prefix: str) -> pd.DataFrame:
    data = pd.read_csv(path, sep='\t', dtype=str)
    indices = [x for x in data.columns if not x.startswith(sequence_prefix)]
    data.set_index(indices, inplace=True)
    data.columns = [
        removeprefix(col, sequence_prefix) for col in data.columns]
    return OpIndexToMulti()(data)


# Defined last takes precedence
@file_reader(FileType.File, FileFormat.Tab)
class TabFileReader(FileReader):
    sequence_prefix = Param(SEQUENCE_PREFIX)

    def call(self, path: Path) -> Iterator[pd.Series]:
        data = readTab(path, sequence_prefix=self.sequence_prefix)
        for col in data:
            yield data[col]


class ReaderNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__((f'No iterator for {str(type)} and {str(format)}'))


def get_reader(type: FileType, format: FileFormat):
    if format not in file_readers[type]:
        raise ReaderNotFound(type, format)
    return file_readers[type][format]


def read_from_path(path: Path) -> Iterator[pd.Series]:
    type, format = autodetect(path)
    reader = get_reader(type, format)
    for series in reader()(path):
        yield OpCheckValid()(series)
