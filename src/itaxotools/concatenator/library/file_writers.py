
from typing import Callable, Dict, TextIO
from pathlib import Path

import pandas as pd

from .file_utils import createDirectory, createZipArchive, PathLike
from .file_types import FileType, FileFormat, get_extension
from .utils import Stream, ConfigurableCallable, Param, Justification
from .operators import (
    OpIndexMerge, OpPadRight, OpDropEmpty, join_any, chain)

from .fasta import fasta_writer
from .phylip import phylip_writer
from .ali import ali_writer
from .nexus import write_from_series as nexus_write


PartWriter = Callable[[pd.Series, TextIO], None]


class FileWriter(ConfigurableCallable):
    type: FileType = None
    format: FileFormat = None

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
        writer.type = type
        writer.format = format
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
    class _SingleFileWriter(FileWriter):
        padding = Param('')

        def call(self, stream: Stream, path: Path) -> None:
            filters = chain([
                OpIndexMerge().to_filter,
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
        padding = Param('')

        def call(self, stream: Stream, path: Path) -> None:
            container = creator(path)
            for series in stream:
                name = series.name + get_extension(FileType.File, format)
                part = container / name
                operator = chain([
                    OpDropEmpty(),
                    OpIndexMerge(),
                    OpPadRight(self.padding),
                    ])
                with part.open('w') as file:
                    writer(operator(series), file)


for file_type, creator in {
    FileType.Directory: createDirectory,
    FileType.ZipArchive: createZipArchive,
}.items():
    for file_format, writer in {
        FileFormat.Fasta: filtered_fasta_writer,
        FileFormat.Phylip: filtered_phylip_writer,
        FileFormat.Ali: filtered_ali_writer,
    }.items():
        _register_multifile_writer(file_type, file_format, creator, writer)


@file_writer(FileType.File, FileFormat.Nexus)
class NexusFileWriter(FileWriter):
    padding = Param('-')
    justification = Param(Justification.Left)
    separator = Param(' ')

    def call(self, stream: Stream, path: Path) -> None:
        data = OpIndexMerge()(join_any(stream))
        generator = (data[col].fillna('') for col in data)
        with path.open('w') as file:
            nexus_write(
                OpPadRight(self.padding).to_filter(generator),
                file,
                self.justification,
                self.separator)


@file_writer(FileType.File, FileFormat.Tab)
class TabFileWriter(FileWriter):
    sequence_prefix = Param('sequence_')

    def call(self, stream: Stream, path: Path) -> None:
        data = join_any(stream)
        data.columns = [self.sequence_prefix + col for col in data.columns]
        with path.open('w') as file:
            data.to_csv(file, sep="\t", line_terminator="\n")


class WriterNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__(f'No writer for {str(type)} and {str(format)}')


def get_writer(type: FileType, format: FileFormat):
    if format not in file_writers[type]:
        raise WriterNotFound(type, format)
    return file_writers[type][format]


def write_to_path(
    series: Stream,
    path: Path,
    type: FileType,
    format: FileFormat,
) -> None:
    writer = get_writer(type, format)
    return writer()(series, path)
