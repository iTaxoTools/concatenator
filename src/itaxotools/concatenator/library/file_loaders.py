
from typing import Callable, Dict
from pathlib import Path

import pandas as pd

from .model import GeneDataFrame
from .utils import ConfigurableCallable
from .file_types import FileType, FileFormat
from .file_identify import autodetect
from .operators import join_any
from .file_readers import file_readers, read_from_path, readNexus, readTab


class FileLoader(ConfigurableCallable):
    def call(self, path: Path) -> GeneDataFrame:
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
    def call(self, path: Path) -> GeneDataFrame:
        df = readTab(path)
        return GeneDataFrame(df, missing='?N', gap='-')


@file_loader(FileType.File, FileFormat.Nexus)
class NexusFileLoader(FileLoader):
    def call(self, path: Path) -> GeneDataFrame:
        df = readNexus(path)
        return GeneDataFrame(df, missing='?N', gap='-')


class LoaderNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__(f'No loader for {str(type)} and {str(format)}')


def get_loader(type: FileType, format: FileFormat):
    if format not in file_loaders[type]:
        raise LoaderNotFound(type, format)
    return file_loaders[type][format]()


def load_from_path(path: Path) -> GeneDataFrame:
    type, format = autodetect(path)
    if format in file_loaders[type]:
        return file_loaders[type][format]()(path)
    elif format in file_readers[type]:
        return join_any(read_from_path(path))
    raise LoaderNotFound(type)
