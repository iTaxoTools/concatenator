
from typing import Callable, Dict, Iterator, TextIO
from pathlib import Path

import pandas as pd

from .file_utils import createDirectory, createZipArchive, PathLike
from .file_types import FileType, FileFormat

from .fasta import fasta_writer
from .phylip import phylip_writer
from .ali import ali_writer
from .nexus import write_from_series

from . import SPECIES, SEQUENCE_PREFIX


CallableWriter = Callable[[Iterator[pd.Series], Path], None]
CallableWriterDecorator = Callable[[CallableWriter], CallableWriter]

file_writers: Dict[FileType, Dict[FileFormat, CallableWriter]] = {
    type: dict() for type in FileType}


# Pending filters
# merge multi-index
# Fasta: remove empty
# Phylip: remove empty, replace Nn- with ?, pad with ?
# Ali: remove empty, pad with N
filtered_fasta_writer = fasta_writer
filtered_phylip_writer = phylip_writer
filtered_ali_writer = ali_writer


def file_writer(
    type: FileType, format: FileFormat
) -> CallableWriterDecorator:
    def decorator(func: CallableWriter) -> CallableWriter:
        file_writers[type][format] = func
        return func
    return decorator


def _writeConcatenatedFormat(
    iterator: Iterator[pd.Series],
    path: Path,
    writer: Callable[[pd.Series, TextIO], None]
) -> None:
    joined = pd.concat(iterator, axis=1)
    data = joined.apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    with path.open('w') as file:
        writer(data, file)


def _register_concatenated_writer(
    format: FileFormat,
    writer: Callable[[pd.Series, TextIO], None],
    filters: Callable[[Iterator[pd.Series]], Iterator[pd.Series]] = list(),
) -> None:

    @file_writer(FileType.File, format)
    def _writeConcatenatedFile(it: Iterator[pd.Series], path: Path) -> None:
        _writeConcatenatedFormat(it, path, writer)


for format, (writer, filters) in {
    FileFormat.Fasta: (filtered_fasta_writer, []),
    FileFormat.Phylip: (filtered_phylip_writer, []),
    FileFormat.Ali: (filtered_ali_writer, []),
}.items():
    _register_concatenated_writer(format, writer, filters)


def _register_multifile_writer(
    type: FileType,
    format: FileFormat,
    creator: Callable[[Path], PathLike],
    writer: Callable[[pd.Series, TextIO], None],
) -> None:

    @file_writer(type, format)
    def _writeMultifile(iterator: Iterator[pd.Series], path: Path) -> None:
        container = creator(path)
        for series in iterator:
            name = format_file_name(series.name, FileType.File, format)
            part = container / name
            with part.open('w') as file:
                writer(series, file)


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
def writeNexusFile(iterator: Iterator[pd.Series], path: Path) -> None:
    data = pd.concat(iterator, axis=1)
    generator = (data[col] for col in data)
    with path.open('w') as file:
        write_from_series(generator, file)


@file_writer(FileType.File, FileFormat.Tab)
def writeTabFile(iterator: Iterator[pd.Series], path: Path) -> None:
    data = pd.concat(iterator, axis=1)
    data.columns = [SEQUENCE_PREFIX + col for col in data.columns]
    data.index.name = SPECIES
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
