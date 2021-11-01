
from typing import Callable, Dict, Iterator, TextIO
from pathlib import Path

import pandas as pd

from .model import GeneSeries, GeneStream, GeneIO
from .utils import ConfigurableCallable, Param
from .file_types import FileFormat, FileType
from .file_utils import ZipPath
from .file_identify import autodetect
from .operators import OpCheckValid

from . import ali, fasta, phylip
from . import nexus, tabfile


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
        return GeneStream(
            iter([self.geneIO.gene_from_path(path)]))


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


@file_reader(FileType.File, FileFormat.Nexus)
class NexusReader(FileReader):
    def call(self, path: Path) -> GeneStream:
        stream = nexus.stream_from_path(path)
        return stream.pipe(OpCheckValid())


@file_reader(FileType.File, FileFormat.Tab)
class TabFileReader(FileReader):
    sequence_prefix = Param('sequence_')

    def call(self, path: Path) -> GeneStream:
        stream = tabfile.stream_from_path(path)
        return stream.pipe(OpCheckValid())


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
