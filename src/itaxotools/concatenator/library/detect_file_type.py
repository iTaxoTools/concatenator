#!/usr/bin/env python3

from typing import Callable, Dict
from enum import Enum, auto
from pathlib import Path
from re import fullmatch

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
def isNexus(path: Path) -> bool:
    with open(path) as file:
        return bool(file.read(6) == '#NEXUS')


@test(TestType.File, FileType.FastaFile)
def isFasta(path: Path) -> bool:
    with open(path) as file:
        return bool(file.read(1) == '>')


@test(TestType.File, FileType.PhylipFile)
def isPhylip(path: Path) -> bool:
    with open(path) as file:
        line = file.readline()
        return bool(fullmatch(r'\W*\d+\W+\d+\W*', line))


@test(TestType.File, FileType.TabFile)
def isTabFile(path: Path) -> bool:
    with open(path) as file:
        line = file.readline()
        return bool(fullmatch(r'.+\t.+\W*', line))


class UnknownFileType(Exception):
    def __init__(self, path):
        self.path = path
        super().__init__(f'Unknown file type for {str(path)}')


class FileNotFound(Exception):
    def __init__(self, path):
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
        tests = TestType.File.tests
    for type, test in tests.items():
        if test(path):
            return type
    raise UnknownFileType(path)
