
from typing import Callable, Dict
from pathlib import Path

from .model import GeneDataFrame
from .utils import ConfigurableCallable
from .file_types import FileType, FileFormat
from .file_identify import autodetect
from .file_readers import file_readers, read_from_path

from . import nexus, tabfile

# This module has little use and should eventually be removed:
# The dataset will not always fit in memory and should be lazily loaded


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
        return tabfile.dataframe_from_path(path)


@file_loader(FileType.File, FileFormat.Nexus)
class NexusFileLoader(FileLoader):
    def call(self, path: Path) -> GeneDataFrame:
        return nexus.dataframe_from_path(path)


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
        return GeneDataFrame.from_stream(read_from_path(path))
    raise LoaderNotFound(type)
