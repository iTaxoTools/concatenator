#!/usr/bin/env python3

from .concatenator import main

from .library.file_types import FileType
from .library.detect_file_type import autodetect
from .library.file_iterators import iterator_from_path
from .library.file_readers import dataframe_from_path
