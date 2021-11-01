
from typing import Callable, Dict, Iterator, TextIO
from pathlib import Path

import pandas as pd

from .model import GeneSeries, GeneStream, GeneIO
from .utils import ConfigurableCallable, Param, removeprefix
from .file_types import FileFormat, FileType
from .file_utils import ZipPath
from .file_identify import autodetect
from .operators import OpCheckValid

from .nexus import read as nexus_read
from . import ali, fasta, phylip


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


class _GeneReader(FileReader):
    geneIO: GeneIO = None

    def read(self, path: Path) -> GeneStream:
        raise NotImplementedError

    def call(self, path: Path) -> GeneStream:
        return self.read(path).pipe(OpCheckValid())


class _SingleFileReader(_GeneReader):
    def read(self, path: Path) -> GeneStream:
        return GeneStream(gene for gene in [self.geneIO.gene_from_path(path)])


class _MultiDirReader(_GeneReader):
    def read(self, path: Path) -> GeneStream:
        return GeneStream(
            self.geneIO.gene_from_path(part)
            for part in path.iterdir())


class _MultiZipReader(_GeneReader):
    def read(self, path: Path) -> GeneStream:
        return GeneStream(
            self.geneIO.gene_from_path(part)
            for part in ZipPath(path).iterdir())


def _register_type_reader(
    ftype: FileType,
    reader: _GeneReader,
) -> None:

    @file_reader(ftype, FileFormat.Fasta)
    class FastaReader(reader):
        geneIO = fasta

    @file_reader(ftype, FileFormat.Phylip)
    class PhylipReader(reader):
        geneIO = phylip

    @file_reader(ftype, FileFormat.Ali)
    class AliReader(reader):
        geneIO = ali


for ftype, reader in {
    FileType.File: _SingleFileReader,
    FileType.Directory: _MultiDirReader,
    FileType.ZipArchive: _MultiZipReader,
}.items():
    _register_type_reader(ftype, reader)


def readNexus(path: Path) -> pd.DataFrame:
    with path.open() as file:
        data = nexus_read(file, sequence_prefix='')
    data.set_index('seqid', inplace=True)
    return data


@file_reader(FileType.File, FileFormat.Nexus)
class NexusReader(FileReader):
    def call(self, path: Path) -> GeneStream:
        data = readNexus(path)
        return GeneStream.from_dataframe(data).pipe(OpCheckValid())


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
        return GeneStream(self.iter(path)).pipe(OpCheckValid())


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
        return GeneStream.from_dataframe(data).pipe(OpCheckValid())


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
    return reader(path)
