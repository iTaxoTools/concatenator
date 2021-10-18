#!/usr/bin/env python3

from pathlib import Path

from .concatenator import main  # noqa

from .library.file_types import FileType, FileFormat
from .library.file_identify import autodetect
from .library.file_readers import read_from_path, get_reader
from .library.file_writers import write_to_path, get_writer, get_extension
from .library.file_loaders import load_from_path
from .library import operators

def convert(
    source: Path,
    dest: Path,
    name: str,
    to_type: FileType,
    to_format: FileFormat,
) -> None:
    in_type, in_format = autodetect(source)
    print(f'Input: {source.name}: {str(in_type)}, {str(in_format)}')
    to_path = dest / (name + get_extension(to_type, to_format))
    write_to_path(read_from_path(source), to_path, to_type, to_format)
    print(f'Output: {to_path.name}: {str(to_type)}, {str(to_format)}')
