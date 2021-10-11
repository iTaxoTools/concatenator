#!/usr/bin/env python3

from itaxotools.concatenator import convert
from pathlib import Path

import sys


if __name__ == '__main__':
    source = Path(sys.argv[1])
    dest = Path(sys.argv[2])
    name = str(sys.argv[3])
    type = getattr(FileType, sys.argv[4])
    format = getattr(FileFormat, sys.argv[5])
    convert(source, dest, name, type, format)
