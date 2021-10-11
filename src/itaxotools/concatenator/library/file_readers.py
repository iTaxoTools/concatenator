
from typing import Callable, Dict
from pathlib import Path

import pandas as pd

from .utils import removeprefix
from .file_types import FileType, FileFormat
from .detect_file_type import autodetect
from .operators import OpIndexToMulti, join_any
from .file_iterators import (
    file_iterators, iterate_path,
    readAliSeries, readFastaSeries, readPhylipSeries,
    readNexusFile as _readNexusFile,
    )

from . import SEQUENCE_PREFIX


CallableReader = Callable[[Path], pd.DataFrame]
CallableReaderDecorator = Callable[[CallableReader], CallableReader]

file_readers: Dict[FileType, Dict[FileFormat, CallableReader]] = {
    type: dict() for type in FileType}


def file_reader(
    type: FileType, format: FileFormat
) -> CallableReaderDecorator:
    def decorator(func: CallableReader) -> CallableReader:
        file_readers[type][format] = func
        return func
    return decorator


@file_reader(FileType.File, FileFormat.Tab)
def readTabFile(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep='\t', dtype=str)
    indices = [x for x in data.columns if not x.startswith(SEQUENCE_PREFIX)]
    data.set_index(indices, inplace=True)
    data.columns = [removeprefix(col, SEQUENCE_PREFIX) for col in data.columns]
    return OpIndexToMulti()(data)


@file_reader(FileType.File, FileFormat.Nexus)
def readNexusFile(path: Path) -> pd.DataFrame:
    return _readNexusFile(path)


@file_reader(FileType.File, FileFormat.Ali)
def readAliFile(path: Path) -> pd.DataFrame:
    return pd.DataFrame(readAliSeries(path))


@file_reader(FileType.File, FileFormat.Fasta)
def readFastaFile(path: Path) -> pd.DataFrame:
    return pd.DataFrame(readFastaSeries(path))


@file_reader(FileType.File, FileFormat.Phylip)
def readPhylipFile(path: Path) -> pd.DataFrame:
    return pd.DataFrame(readPhylipSeries(path))


class ReaderNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__(f'No reader for {str(type)} and {str(format)}')


def dataframe_from_path(path: Path) -> pd.DataFrame:
    """Species as index, sequences as columns"""
    type, format = autodetect(path)
    if format in file_readers[type]:
        return file_readers[type][format](path)
    elif format in file_iterators[type]:
        return join_any(iterate_path(path))
    raise ReaderNotFound(type)
