#!/usr/bin/env python3

from itaxotools.concatenator.library.file_types import FileType, FileFormat
from itaxotools.concatenator.library.file_readers import read_from_path
from itaxotools.concatenator.library.detect_file_type import autodetect
from itaxotools.concatenator.library.file_writers import write_from_iterator, format_file_name

from pathlib import Path

import sys


def convert(
    source: Path,
    dest: Path,
    name: str,
    to_type: FileType,
    to_format: FileFormat,
) -> None:
    in_type, in_format = autodetect(source)
    print(f'Input: {source.name}: {in_type.description}, {in_format.description}')
    to_path = dest / format_file_name(name, to_type, to_format)
    write_from_iterator(read_from_path(source), to_path, to_type, to_format)
    print(f'Output: {to_path.name}: {to_type.description}, {to_format.description}')

if __name__ == '__main__':
    source = Path(sys.argv[1])
    dest = Path(sys.argv[2])
    name = str(sys.argv[3])
    type = getattr(FileType, sys.argv[4])
    format = getattr(FileFormat, sys.argv[5])
    convert(source, dest, name, type, format)
