#!/usr/bin/env python3

from itaxotools.concatenator import FileType, FileFormat, convert
from pathlib import Path

import sys


if __name__ == '__main__':
    if len(sys.argv) < 6:
        print(' Usage: convert.py SOURCE DEST NAME TYPE FORMAT')
        print('    ex: convert.py examples/concat.tab examples new File Fasta')
        print('  TYPE: FROM: File Directory ZipArchive')
        print('FORMAT: FROM: Ali Fasta Phylip Nexus Tab')
        exit()
    source = Path(sys.argv[1])
    dest = Path(sys.argv[2])
    name = str(sys.argv[3])
    type = getattr(FileType, sys.argv[4])
    format = getattr(FileFormat, sys.argv[5])
    convert(source, dest, name, type, format)
