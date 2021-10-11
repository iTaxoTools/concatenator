
from typing import Callable, Dict, TextIO
from pathlib import Path

import pandas as pd

from .file_utils import createDirectory, createZipArchive, PathLike
from .file_types import FileType, FileFormat
from .utils import Stream, ConfigurableCallable, Param
from .operators import (
    OpIndexMerge, OpPadRight, OpIndexToMulti, OpDropEmpty, join_any, chain)

from .fasta import fasta_writer
from .phylip import phylip_writer
from .ali import ali_writer
from .nexus import write_from_series

from . import SEQUENCE_PREFIX


PartWriter = Callable[[pd.Series, TextIO], None]


class FileWriter(ConfigurableCallable):
    def call(self, stream: Stream, path: Path) -> None:
        raise NotImplementedError


file_writers: Dict[FileType, Dict[FileFormat, FileWriter]] = {
    type: dict() for type in FileType}


# Pending filters
filtered_fasta_writer = fasta_writer
filtered_phylip_writer = phylip_writer
filtered_ali_writer = ali_writer


def file_writer(
    type: FileType, format: FileFormat
) -> Callable[[FileWriter], FileWriter]:
    def decorator(writer: FileWriter) -> FileWriter:
        file_writers[type][format] = writer
        return writer
    return decorator


def _writeConcatenatedFormat(
    stream: Stream,
    path: Path,
    writer: PartWriter,
) -> None:
    joined = join_any(stream)
    data = joined.apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    with path.open('w') as file:
        writer(data, file)


def _register_concatenated_writer(
    format: FileFormat,
    writer: PartWriter,
) -> None:

    @file_writer(FileType.File, format)
    class _ConcatenatedFileWriter(FileWriter):
        padding = Param('-')

        def call(self, stream: Stream, path: Path) -> None:
            filters = chain([
                OpIndexMerge().to_filter,
                OpIndexToMulti().to_filter,
                OpPadRight(self.padding).to_filter,
                ])
            _writeConcatenatedFormat(filters(stream), path, writer)


for format, (writer, filters) in {
    FileFormat.Fasta: (filtered_fasta_writer, []),
    FileFormat.Phylip: (filtered_phylip_writer, []),
    FileFormat.Ali: (filtered_ali_writer, []),
}.items():
    _register_concatenated_writer(format, writer)


def _register_multifile_writer(
    type: FileType,
    format: FileFormat,
    creator: Callable[[Path], PathLike],
    writer: PartWriter,
) -> None:

    @file_writer(type, format)
    class _MultiFileWriter(FileWriter):
        padding = Param('-')

        def call(self, stream: Stream, path: Path) -> None:
            container = creator(path)
            for series in stream:
                name = format_file_name(series.name, FileType.File, format)
                part = container / name
                operator = chain([
                    OpDropEmpty(),
                    OpIndexMerge(),
                    OpIndexToMulti(),
                    OpPadRight(self.padding),
                    ])
                with part.open('w') as file:
                    writer(operator(series), file)


for type, creator in {
    FileType.Directory: createDirectory,
    FileType.ZipArchive: createZipArchive,
}.items():
    for format, writer in {
        FileFormat.Fasta: filtered_fasta_writer,
        FileFormat.Phylip: filtered_phylip_writer,
        FileFormat.Ali: filtered_ali_writer,
    }.items():
        _register_multifile_writer(type, format, creator, writer)


@file_writer(FileType.File, FileFormat.Nexus)
class NexusFileWriter(FileWriter):
    padding = Param('*')

    def call(self, stream: Stream, path: Path) -> None:
        data = OpIndexMerge()(join_any(stream))
        generator = (data[col].fillna('') for col in data)
        with path.open('w') as file:
            write_from_series(
                OpPadRight(self.padding).to_filter(generator), file)


@file_writer(FileType.File, FileFormat.Tab)
class TabFileWriter(FileWriter):
    def call(self, stream: Stream, path: Path) -> None:
        data = join_any(stream)
        data.columns = [SEQUENCE_PREFIX + col for col in data.columns]
        with path.open('w') as file:
            data.to_csv(file, sep="\t", line_terminator="\n")


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


def write_from_stream(
    series: Stream,
    path: Path,
    type: FileType,
    format: FileFormat,
) -> None:
    """Species as index, sequences as columns"""
    if format not in file_writers[type]:
        raise WriterNotFound(type)
    return file_writers[type][format]()(series, path)
