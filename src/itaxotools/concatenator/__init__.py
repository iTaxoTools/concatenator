#!/usr/bin/env python3

from pathlib import Path

from .concatenator import main  # noqa

from .library.file_types import FileType, FileFormat
from .library.file_identify import autodetect
from .library.file_readers import read_from_path
from .library.file_loaders import load_from_path
from .library.file_writers import write_from_stream, format_file_name


def convert(
    source: Path,
    dest: Path,
    name: str,
    to_type: FileType,
    to_format: FileFormat,
) -> None:
    in_type, in_format = autodetect(source)
    print(f'Input: {source.name}: {str(in_type)}, {str(in_format)}')
    to_path = dest / format_file_name(name, to_type, to_format)
    write_from_stream(read_from_path(source), to_path, to_type, to_format)
    print(f'Output: {to_path.name}: {str(to_type)}, {str(to_format)}')
