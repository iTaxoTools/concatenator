
from typing import Callable, Dict
from pathlib import Path

import pandas as pd

from .utils import ConfigurableCallable, removeprefix
from .file_types import FileType, FileFormat
from .detect_file_type import autodetect
from .operators import OpIndexToMulti, join_any
from .file_iterators import (
    file_iterators, iterate_path,
    readAliSeries, readFastaSeries, readPhylipSeries,
    readNexusFile as _readNexusFile,
    )

from . import SEQUENCE_PREFIX


class FileLoader(ConfigurableCallable):
    def call(self, path: Path) -> pd.DataFrame:
        raise NotImplementedError


file_loaders: Dict[FileType, Dict[FileFormat, FileLoader]] = {
    type: dict() for type in FileType}


def file_loader(
    type: FileType, format: FileFormat
) -> Callable[[FileLoader], FileLoader]:
    def decorator(loader: FileLoader) -> FileLoader:
        file_loaders[type][format] = loader
        return loader
    return decorator


@file_loader(FileType.File, FileFormat.Tab)
class TabFileLoader(FileLoader):
    def call(self, path: Path) -> pd.DataFrame:
        data = pd.read_csv(path, sep='\t', dtype=str)
        indices = [x for x in data.columns if not x.startswith(SEQUENCE_PREFIX)]
        data.set_index(indices, inplace=True)
        data.columns = [
            removeprefix(col, SEQUENCE_PREFIX) for col in data.columns]
        return OpIndexToMulti()(data)


@file_loader(FileType.File, FileFormat.Nexus)
class NexusFileLoader(FileLoader):
    def call(self, path: Path) -> pd.DataFrame:
        return _readNexusFile(path)


@file_loader(FileType.File, FileFormat.Ali)
class AliFileLoader(FileLoader):
    def call(self, path: Path) -> pd.DataFrame:
        return pd.DataFrame(readAliSeries(path))


@file_loader(FileType.File, FileFormat.Fasta)
class FastaFileLoader(FileLoader):
    def call(self, path: Path) -> pd.DataFrame:
        return pd.DataFrame(readFastaSeries(path))


@file_loader(FileType.File, FileFormat.Phylip)
class PhylipFileLoader(FileLoader):
    def call(self, path: Path) -> pd.DataFrame:
        return pd.DataFrame(readPhylipSeries(path))


class LoaderNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__(f'No loader for {str(type)} and {str(format)}')


def load_from_path(path: Path) -> pd.DataFrame:
    type, format = autodetect(path)
    if format in file_loaders[type]:
        return file_loaders[type][format]()(path)
    elif format in file_iterators[type]:
        return join_any(iterate_path(path))
    raise LoaderNotFound(type)
