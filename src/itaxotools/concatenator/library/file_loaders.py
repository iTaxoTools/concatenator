
from typing import Callable, Dict
from pathlib import Path

import pandas as pd

from .utils import ConfigurableCallable
from .file_types import FileType, FileFormat
from .file_identify import autodetect
from .operators import join_any
from .file_readers import (
    file_readers, read_from_path, readNexusFile, readTab,
    readAliSeries, readFastaSeries, readPhylipSeries,
    )


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
        return readTab(path)


@file_loader(FileType.File, FileFormat.Nexus)
class NexusFileLoader(FileLoader):
    def call(self, path: Path) -> pd.DataFrame:
        return readNexusFile(path)


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
    elif format in file_readers[type]:
        return join_any(read_from_path(path))
    raise LoaderNotFound(type)
