#!/usr/bin/env python3

from typing import Callable, Dict
from enum import Enum, auto
from pathlib import Path
from re import fullmatch
from zipfile import ZipFile, is_zipfile
from zipp import Path as ZipPath  # BUGFIX: backport from Python 3.9.1

from .file_types import FileFormat


CallableTest = Callable[[Path], bool]


class TestType(Enum):
    File = auto()
    Directory = auto()
    Archive = auto()

    def __init__(self, value: int):
        self.tests: Dict[FileFormat, CallableTest] = dict()


def test(test: TestType, type: FileFormat) -> CallableTest:
    def decorator(func: CallableTest) -> CallableTest:
        test.tests[type] = func
        return func
    return decorator


@test(TestType.File, FileFormat.NexusFile)
def isNexusFile(path: Path) -> bool:
    with path.open() as file:
        return bool(file.read(6) == '#NEXUS')


@test(TestType.File, FileFormat.FastaFile)
def isFastaFile(path: Path) -> bool:
    with path.open() as file:
        return bool(file.read(1) == '>')


@test(TestType.File, FileFormat.PhylipFile)
def isPhylipFile(path: Path) -> bool:
    with path.open() as file:
        line = file.readline()
        return bool(fullmatch(r'\s*\d+\s+\d+\s*', line))


@test(TestType.File, FileFormat.TabFile)
def isTabFile(path: Path) -> bool:
    with path.open() as file:
        line = file.readline()
        return bool(fullmatch(r'([^\t]+\t)+[^\t]+', line))


@test(TestType.File, FileFormat.AliFile)
def isAliFile(path: Path) -> bool:
    with path.open() as file:
        for line in file:
            if line[0] in ['#', '\n']:
                continue
            return bool(line[0] == '>')


def _archiveTest(path: Path, test: CallableTest) -> bool:
    archive = ZipFile(path, 'r')
    for name in archive.namelist():
        path = ZipPath(archive, name)
        if not path.is_file():
            return False
        if not test(path):
            return False
    return True


@test(TestType.Archive, FileFormat.MultiFastaInput)
def isFastaArchive(path: Path) -> bool:
    return _archiveTest(path, isFastaFile)


@test(TestType.Archive, FileFormat.MultiPhylipInput)
def isPhylipArchive(path: Path) -> bool:
    return _archiveTest(path, isPhylipFile)


@test(TestType.Archive, FileFormat.MultiAliInput)
def isAliArchive(path: Path) -> bool:
    return _archiveTest(path, isAliFile)


# def _directoryTest(directory: Path, test: CallableTest) -> bool:
#     for path in directory.glob('*'):
#         if not path.is_file():
#             return False
#         if not test(path):
#             return False
#     return True
#
#
# @test(TestType.Directory, FileFormat.MultiFasta)
# def isFastaDirectory(path: Path) -> bool:
#     return _directoryTest(path, isFastaFile)


class UnknownFileFormat(Exception):
    def __init__(self, path: Path):
        self.path = path
        super().__init__(f'Unknown file type for {str(path)}')


class FileNotFound(Exception):
    def __init__(self, path: Path):
        self.path = path
        super().__init__(f'File not found: {str(path)}')


def autodetect(path: Path) -> FileFormat:
    """Attempt to automatically determine sequence file type"""
    if not path.exists():
        raise FileNotFound(path)
    tests = {}
    if path.is_dir():
        tests = TestType.Directory.tests
    elif path.is_file():
        if is_zipfile(path):
            tests = TestType.Archive.tests
        else:
            tests = TestType.File.tests
    for type, test in tests.items():
        if test(path):
            return type
    raise UnknownFileFormat(path)
