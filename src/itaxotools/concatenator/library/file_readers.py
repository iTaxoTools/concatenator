
from typing import Callable, Dict, TextIO
from pathlib import Path

import pandas as pd

from .file_types import FileType
from .detect_file_type import autodetect
from .nexus import read as nexus_read
from .ali import column_reader as ali_reader
from .fasta import column_reader as fasta_reader
from .phylip import column_reader as phylip_reader


CallableReader = Callable[[Path], pd.DataFrame]

type_readers: Dict[FileType, CallableReader] = dict()


CallableReaderDecorator = Callable[[CallableReader], CallableReader]
type_reader(type: FileType) -> CallableReaderDecorator:
    def decorator(func: CallableReader) -> CallableReader:
        type_readers[type] = func
        return func
    return decorator


@type_reader(FileType.TabFile)
def readTabFile(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    data.drop(columns=['specimen-voucher', 'locality'], inplace=True)
    data.set_index(data.loc[:, 'species'])
    data.drop(columns=['species'], inplace=True)
    data.columns = [c.removeprefix('sequence_') for c in data.columns]
    return data


@type_reader(FileType.NexusFile)
def readNexusFile(path: Path) -> pd.DataFrame:
    with path.open() as file:
        data = nexus_read(file)
    data.set_index('seqid', inplace=True)
    data.index.name = None
    return data


def _readSeries(
    path: Path,
    func: Callable[[TextIO], pd.Series]
) -> pd.Series:
    with path.open() as file:
        series = func(file)
    series.name = path.stem
    return series


def readAliSeries(path: Path) -> pd.Series:
    series = _readSeries(path, ali_reader)
    return series.str.replace('_', '-', regex=False)


def readFastaSeries(path: Path) -> pd.Series:
    return _readSeries(path, fasta_reader)


def readPhylipSeries(path: Path) -> pd.Series:
    return _readSeries(path, phylip_reader)


@type_reader(FileType.AliFile)
def readAliFile(path: Path) -> pd.DataFrame:
    return pd.DataFrame(readAliSeries(path))


@type_reader(FileType.FastaFile)
def readFastaFile(path: Path) -> pd.DataFrame:
    return pd.DataFrame(readFastaSeries(path))


@type_reader(FileType.PhylipFile)
def readPhylipFile(path: Path) -> pd.DataFrame:
    return pd.DataFrame(readPhylipSeries(path))


class ReaderNotFound(Exception):
    def __init__(self, type: FileType):
        self.type = type
        super().__init__(f'No reader for FileType: {str(type)}')


def dataframe_from_path(path: Path) -> pd.DataFrame:
    """Species as index, sequences as columns"""
    type = autodetect(path)
    if type not in type_readers:
        raise ReaderNotFound(type)
    data = type_readers[type](path)
    return data
