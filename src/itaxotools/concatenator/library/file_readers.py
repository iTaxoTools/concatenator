
from typing import Callable, Dict, Iterator, TextIO
from pathlib import Path

import pandas as pd

from .model import GeneSeries, GeneStream
from .utils import ConfigurableCallable, Param, removeprefix
from .file_types import FileFormat, FileType
from .file_utils import ZipPath
from .file_identify import autodetect
from .operators import OpCheckValid

from .nexus import read as nexus_read
from .ali import ali_reader
from .fasta import fasta_reader
from .phylip import phylip_reader


class FileReader(ConfigurableCallable):
    type: FileType = None
    format: FileFormat = None

    def call(self, path: Path) -> GeneStream:
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


def readNexus(path: Path) -> pd.DataFrame:
    with path.open() as file:
        data = nexus_read(file, sequence_prefix='')
    data.set_index('seqid', inplace=True)
    return data


@file_reader(FileType.File, FileFormat.Nexus)
class NexusReader(FileReader):
    def call(self, path: Path) -> GeneStream:
        data = readNexus(path)
        return GeneStream.from_dataframe(data)


def _readSeries(
    path: Path,
    part_reader: Callable[[TextIO], pd.Series]
) -> pd.Series:
    with path.open() as file:
        series = part_reader(file)
    series.name = path.stem
    return series


def readAliGene(path: Path) -> GeneSeries:
    series =  _readSeries(path, ali_reader)
    return GeneSeries(series, missing='?', gap='*')


def readFastaGene(path: Path) -> GeneSeries:
    series = _readSeries(path, fasta_reader)
    return GeneSeries(series, missing='?N', gap='-')


def readPhylipGene(path: Path) -> GeneSeries:
    series = _readSeries(path, phylip_reader)
    return GeneSeries(series, missing='?N', gap='-')


def _register_multifile_reader(
    format: FileFormat,
    reader: Callable[[Path], GeneSeries]
) -> None:

    @file_reader(FileType.File, format)
    class _SingleFileReader(FileReader):
        def call(self, path: Path) -> GeneStream:
            return GeneStream(gene for gene in [reader(path)])

    @file_reader(FileType.Directory, format)
    class _MultiDirReader(FileReader):
        def call(self, path: Path) -> GeneStream:
            return GeneStream(
                (reader(part) for part in path.iterdir()))

    @file_reader(FileType.ZipArchive, format)
    class _MultiZipReader(FileReader):
        def call(self, path: Path) -> GeneStream:
            return GeneStream(
                (reader(part) for part in ZipPath(path).iterdir()))


for format, reader in {
    FileFormat.Fasta: readFastaGene,
    FileFormat.Phylip: readPhylipGene,
    FileFormat.Ali: readAliGene,
}.items():
    _register_multifile_reader(format, reader)


@file_reader(FileType.File, FileFormat.Tab)
class TabFileReaderSlow(FileReader):
    sequence_prefix = Param('sequence_')

    def iter(self, path: Path) -> Iterator[GeneSeries]:
        with path.open() as file:
            columns = file.readline().rstrip().split("\t")
            sequences = [x for x in columns if x.startswith(self.sequence_prefix)]
            indices = [x for x in columns if not x.startswith(self.sequence_prefix)]
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
                series.name = removeprefix(sequence, self.sequence_prefix)
                yield GeneSeries(series)

    def call(self, path: Path) -> GeneStream:
        return GeneStream(self.iter(path))


def readTab(path: Path, sequence_prefix: str = 'sequence_') -> pd.DataFrame:
    data = pd.read_csv(path, sep='\t', dtype=str)
    indices = [x for x in data.columns if not x.startswith(sequence_prefix)]
    data.set_index(indices, inplace=True)
    data.columns = [
        removeprefix(col, sequence_prefix) for col in data.columns]
    return data


# Defined last takes precedence
@file_reader(FileType.File, FileFormat.Tab)
class TabFileReader(FileReader):
    sequence_prefix = Param('sequence_')

    def call(self, path: Path) -> GeneStream:
        data = readTab(path, sequence_prefix=self.sequence_prefix)
        return GeneStream.from_dataframe(data)


class ReaderNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__((f'No iterator for {str(type)} and {str(format)}'))


def get_reader(type: FileType, format: FileFormat):
    if format not in file_readers[type]:
        raise ReaderNotFound(type, format)
    return file_readers[type][format]()


def read_from_path(path: Path) -> GeneStream:
    type, format = autodetect(path)
    reader = get_reader(type, format)
    return reader(path).pipe(OpCheckValid())
