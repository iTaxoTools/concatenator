#!/usr/bin/env python3

from typing import Callable, Dict, Tuple, Iterator
from pathlib import Path
from re import fullmatch
from zipfile import is_zipfile

from .file_types import FileFormat, FileType
from .file_utils import iterateZipArchive, iterateDirectory


CallableTest = Callable[[Path], bool]
CallableTestDecorator = Callable[[CallableTest], CallableTest]

tests: Dict[FileType, Dict[FileFormat, CallableTest]] = {
    type: dict() for type in FileType}


def test(type: FileType, format: FileFormat) -> CallableTestDecorator:
    def decorator(func: CallableTest) -> CallableTest:
        tests[type][format] = func
        return func
    return decorator


@test(FileType.File, FileFormat.Nexus)
def isNexusFile(path: Path) -> bool:
    with path.open() as file:
        return bool(file.read(6) == '#NEXUS')


@test(FileType.File, FileFormat.Fasta)
def isFastaFile(path: Path) -> bool:
    with path.open() as file:
        return bool(file.read(1) == '>')


@test(FileType.File, FileFormat.Phylip)
def isPhylipFile(path: Path) -> bool:
    with path.open() as file:
        line = file.readline()
        return bool(fullmatch(r'\s*\d+\s+\d+\s*', line))


@test(FileType.File, FileFormat.Tab)
def isTabFile(path: Path) -> bool:
    with path.open() as file:
        line = file.readline()
        return bool(fullmatch(r'([^\t]+\t)+[^\t]+', line))


@test(FileType.File, FileFormat.Ali)
def isAliFile(path: Path) -> bool:
    # Also catches Fasta format, check this last
    with path.open() as file:
        for line in file:
            if line[0] in ['#', '\n']:
                continue
            return bool(line[0] == '>')
    return False


def _containerTest(
    parts: Iterator[Path],
    test: CallableTest,
) -> bool:
    for part in parts:
        if not part.is_file():
            return False
        if not test(part):
            return False
    return True


def _register_multifile_test(
    format: FileFormat,
    tester: CallableTest
) -> None:

    @test(FileType.Directory, format)
    def _testMultifileDir(path: Path) -> bool:
        return _containerTest(iterateDirectory(path), tester)

    @test(FileType.ZipArchive, format)
    def _testMultifileZip(path: Path) -> bool:
        return _containerTest(iterateZipArchive(path), tester)


for format, tester in {
    FileFormat.Fasta: isFastaFile,
    FileFormat.Phylip: isPhylipFile,
    FileFormat.Ali: isAliFile,
}.items():
    _register_multifile_test(format, tester)


isAliZip = tests[FileType.ZipArchive][FileFormat.Ali]
isFastaZip = tests[FileType.ZipArchive][FileFormat.Fasta]
isPhylipZip = tests[FileType.ZipArchive][FileFormat.Phylip]


class UnknownFileType(Exception):
    def __init__(self, path: Path):
        self.path = path
        super().__init__(f'Unknown FileType for {str(path)}')


class UnknownFileFormat(Exception):
    def __init__(self, path: Path):
        self.path = path
        super().__init__(f'Unknown FileFormat for {str(path)}')


class FileNotFound(Exception):
    def __init__(self, path: Path):
        self.path = path
        super().__init__(f'File not found: {str(path)}')


def autodetect_type(path: Path) -> FileType:
    if path.is_dir():
        return FileType.Directory
    elif path.is_file():
        if is_zipfile(path):
            return FileType.ZipArchive
        else:
            return FileType.File
    raise UnknownFileType(path)


def autodetect_format(path: Path, type: FileType) -> FileFormat:
    for format in tests[type]:
        if tests[type][format](path):
            return format
    raise UnknownFileFormat(path)


def autodetect(path: Path) -> Tuple[FileType, FileFormat]:
    """Attempt to automatically determine sequence file type"""
    if not path.exists():
        raise FileNotFound(path)
    type = autodetect_type(path)
    format = autodetect_format(path, type)
    return (type, format)
