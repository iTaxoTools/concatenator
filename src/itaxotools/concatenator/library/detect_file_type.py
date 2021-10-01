#!/usr/bin/env python3

from typing import Callable, Dict
from enum import Enum, auto
from pathlib import Path
from re import fullmatch
from zipfile import ZipFile, is_zipfile
from zipp import Path as ZipPath  # BUGFIX: backport from Python 3.9.1

from .file_types import FileType


CallableTest = Callable[[Path], bool]


class TestType(Enum):
    File = auto()
    Directory = auto()
    Archive = auto()

    def __init__(self, value: int):
        self.tests: Dict[FileType, CallableTest] = dict()


def test(test: TestType, type: FileType) -> CallableTest:
    def decorator(func: CallableTest) -> CallableTest:
        test.tests[type] = func
        return func
    return decorator


@test(TestType.File, FileType.NexusFile)
def isNexusFile(path: Path) -> bool:
    with path.open() as file:
        return bool(file.read(6) == '#NEXUS')


@test(TestType.File, FileType.FastaFile)
def isFastaFile(path: Path) -> bool:
    with path.open() as file:
        return bool(file.read(1) == '>')


@test(TestType.File, FileType.PhylipFile)
def isPhylipFile(path: Path) -> bool:
    with path.open() as file:
        line = file.readline()
        return bool(fullmatch(r'\s*\d+\s+\d+\s*', line))


@test(TestType.File, FileType.TabFile)
def isTabFile(path: Path) -> bool:
    with path.open() as file:
        line = file.readline()
        return bool(fullmatch(r'([^\t]+\t)+[^\t]+', line))


# Can only be part of an archive
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


@test(TestType.Archive, FileType.MultiFastaInput)
def isFastaArchive(path: Path) -> bool:
    return _archiveTest(path, isFastaFile)


@test(TestType.Archive, FileType.MultiPhylipInput)
def isPhylipArchive(path: Path) -> bool:
    return _archiveTest(path, isPhylipFile)


@test(TestType.Archive, FileType.MultiAliInput)
def isAliArchive(path: Path) -> bool:
    return _archiveTest(path, isAliFile)


class UnknownFileType(Exception):
    def __init__(self, path: Path):
        self.path = path
        super().__init__(f'Unknown file type for {str(path)}')


class FileNotFound(Exception):
    def __init__(self, path: Path):
        self.path = path
        super().__init__(f'File not found: {str(path)}')


def autodetect(path: Path) -> FileType:
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
    raise UnknownFileType(path)
