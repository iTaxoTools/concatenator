
from typing import Callable, Dict, TextIO
from pathlib import Path

import pandas as pd

from .model import GeneSeries, GeneDataFrame, GeneStream, GeneIO
from .types import Justification
from .utils import ConfigurableCallable, Field
from .file_utils import ZipFile, ZipPath
from .file_types import FileType, FileFormat, get_extension
from .operators import (
    OpCheckValid, OpIndexMerge, OpPadRight, OpDropEmpty, OpApplyToSeries)

from . import ali, fasta, phylip
from . import nexus, tabfile


class FileWriter(ConfigurableCallable):
    type: FileType = None
    format: FileFormat = None

    def filter(self, stream: GeneStream) -> GeneStream:
        """Stream operations, obeys inheritance"""
        return stream.pipe(OpCheckValid())

    def call(self, stream: GeneStream, path: Path) -> None:
        """Write to path after applying filter operations"""
        raise NotImplementedError


file_writers: Dict[FileType, Dict[FileFormat, FileWriter]] = {
    type: dict() for type in FileType}


def file_writer(
    type: FileType, format: FileFormat
) -> Callable[[FileWriter], FileWriter]:
    def decorator(writer: FileWriter) -> FileWriter:
        file_writers[type][format] = writer
        writer.type = type
        writer.format = format
        return writer
    return decorator


class _GeneWriter(FileWriter):
    geneIO: GeneIO = None

    def write(self, stream: GeneStream, path: Path) -> None:
        raise NotImplementedError

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        self.write(stream, path)


class _ConcatenatedWriter(_GeneWriter):
    padding = Field('padding', value='')

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = (
            super().filter(stream)
            .pipe(OpIndexMerge())
            .pipe(OpPadRight(self.padding)))
        return stream

    def write(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        joined = GeneDataFrame.from_stream(stream)
        data = joined.dataframe.apply(
            lambda row: ''.join(row.values.astype(str)), axis=1)
        gene = GeneSeries(data, missing='', gap='')
        self.geneIO.gene_to_path(gene, path)


class _MultiFileWriter(_GeneWriter):
    padding = Field('padding', value='')

    @staticmethod
    def create(path: Path) -> Path:
        raise NotImplementedError

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = (
            super().filter(stream)
            .pipe(OpDropEmpty())
            .pipe(OpIndexMerge())
            .pipe(OpPadRight(self.padding)))
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        container = self.create(path)
        for gene in stream:
            name = gene.series.name + get_extension(FileType.File, self.format)
            part = container / name
            self.geneIO.gene_to_path(gene, part)


class _MultiDirWriter(_MultiFileWriter):
    @staticmethod
    def create(path: Path) -> Path:
        path.mkdir(exist_ok=True)
        return path


class _MultiZipWriter(_MultiFileWriter):
    @staticmethod
    def create(path: Path) -> Path:
        archive = ZipFile(path, 'w')
        return ZipPath(archive)


def _register_type_writer(
    ftype: FileType,
    writer: _GeneWriter,
) -> None:

    @file_writer(ftype, FileFormat.Fasta)
    class FastaWriter(writer):
        geneIO = fasta

    @file_writer(ftype, FileFormat.Phylip)
    class PhylipWriter(writer):
        geneIO = phylip

    @file_writer(ftype, FileFormat.Ali)
    class AliWriter(writer):
        geneIO = ali


for ftype, writer in {
    FileType.File: _ConcatenatedWriter,
    FileType.Directory: _MultiDirWriter,
    FileType.ZipArchive: _MultiZipWriter,
}.items():
    _register_type_writer(ftype, writer)


@file_writer(FileType.File, FileFormat.Nexus)
class NexusWriter(FileWriter):
    padding = Field('padding', value='-')
    justification = Field('justification', value=Justification.Left)
    separator = Field('separator', value=' ')

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = super().filter(stream)
        stream = (
            GeneDataFrame.from_stream(stream).to_stream()
            .pipe(OpIndexMerge())
            .pipe(OpApplyToSeries(lambda x: x.fillna('')))
            .pipe(OpPadRight(self.padding)))
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        nexus.stream_to_path(stream, path, self.justification, self.separator)


@file_writer(FileType.File, FileFormat.Tab)
class TabWriter(FileWriter):
    sequence_prefix = Field('sequence_prefix', value='sequence_')

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        tabfile.stream_to_path(stream, path, self.sequence_prefix)


class WriterNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__(f'No writer for {str(type)} and {str(format)}')


def get_writer(type: FileType, format: FileFormat):
    if format not in file_writers[type]:
        raise WriterNotFound(type, format)
    return file_writers[type][format]()


def write_to_path(
    series: GeneStream,
    path: Path,
    type: FileType,
    format: FileFormat,
) -> None:
    writer = get_writer(type, format)
    return writer(series, path)
