#!/usr/bin/env python3

from typing import Optional, Callable, TextIO, Any, Union, List, Tuple, BinaryIO
from enum import Enum
import os
import shutil

import library.concat as concat
import library.tab_to_nexus as tab_to_nexus
import library.nexus_to_tab as nexus_to_tab
import library.dnaconvert as dnaconvert
import library.tab_to_multifile as tab_to_multifile
from library.file_utils import make_binary, make_binary_in


class FileType(Enum):
    """Types of files that can be processed by operations"""

    TabFile = ("Tab file", ".tab", False)
    NexusFile = ("NEXUS file", ".nex", False)
    ConcatTabFile = ("Concatenated Tab file", ".tab", False)
    ConcatFasta = ("Concatenated FASTA file", ".fas", False)
    ConcatPhylip = ("Concatenated Phylip file", ".phy", False)
    MultiFasta = ("Multifile FASTA archive", ".zip", True)
    MultiPhylip = ("Multifile Phylip archive", ".zip", True)
    MultiAli = ("Multifile Ali archive", ".zip", True)

    def __init__(self, description: str, extension: str, timestamp: bool):
        self.description = description
        self.extension = extension
        self.timestamp = timestamp


Parameter = Union[None, int, str]


def dnaconvert_wrapper(format: str) -> Callable[[TextIO, TextIO], None]:
    def wrap(infile: TextIO, outfile: TextIO) -> None:
        dnaconvert.process(infile, outfile, format)

    return wrap


class Operation(Enum):
    """Operations that the program can perform on files"""

    TabToNexus = (FileType.TabFile, "Tabfile to NEXUS", None, None)
    NexusToTab = (FileType.NexusFile, "NEXUS to Tabfile", None, None)
    ConcatTab = (FileType.TabFile, "Concatenate", None, None)
    DnaConvert = (
        FileType.ConcatTabFile,
        "DNAconvert\n(concatenated)",
        "Convert into {} with DNAconvert",
        ["fasta", "phylip"],
    )
    TabToMFasta = (FileType.TabFile, "Tabfile to multifile Fasta", None, None)
    TabToMPhylip = (FileType.TabFile, "Tabfile to multifile Phylip", None, None)
    TabToMAli = (FileType.TabFile, "Tabfile to multifile Ali", None, None)

    def __init__(
        self,
        input_type: FileType,
        button_text: str,
        description: Optional[str],
        parameter_type: Any,
    ):
        # parameter_type should be one of `None`, `int`, `str` or a list of strings
        self.input_type = input_type
        self.button_text = button_text
        if description:
            self.description = description
        else:
            self.description = button_text
        self.parameter_type = parameter_type

    def output_type(self, parameter: Parameter) -> FileType:
        if self == Operation.TabToNexus:
            return FileType.NexusFile
        elif self == Operation.NexusToTab:
            return FileType.TabFile
        elif self == Operation.ConcatTab:
            return FileType.ConcatTabFile
        elif self == Operation.DnaConvert:
            if parameter == "fasta":
                return FileType.ConcatFasta
            elif parameter == "phylip":
                return FileType.ConcatPhylip
            else:
                assert False
        elif self == Operation.TabToMFasta:
            return FileType.MultiFasta
        elif self == Operation.TabToMPhylip:
            return FileType.MultiPhylip
        elif self == Operation.TabToMAli:
            return FileType.MultiAli
        else:
            assert False

    def apply(self, parameter: Parameter) -> Callable[[BinaryIO, BinaryIO], None]:
        if self == Operation.TabToNexus:
            return make_binary(tab_to_nexus.process)
        elif self == Operation.NexusToTab:
            return make_binary(nexus_to_tab.process)
        elif self == Operation.ConcatTab:
            return make_binary(concat.process)
        elif self == Operation.DnaConvert:
            if type(parameter) is str:
                return make_binary(dnaconvert_wrapper(parameter))
            else:
                assert False
        elif self == Operation.TabToMFasta:
            return make_binary_in(tab_to_multifile.process_fasta)
        elif self == Operation.TabToMPhylip:
            return make_binary_in(tab_to_multifile.process_phylip)
        elif self == Operation.TabToMAli:
            return make_binary_in(tab_to_multifile.process_ali)
        else:
            assert False


def run_pipeline(
    infile: BinaryIO, outfile: BinaryIO, pipeline: List[Tuple[Operation, Parameter]]
) -> None:
    if not pipeline:
        shutil.copyfileobj(infile, outfile)
        return
    read_file = infile
    file_pairs = []
    for _ in range(len(pipeline) - 1):
        (r, w) = os.pipe()
        write_file = open(w, mode="wb")
        file_pairs.append((read_file, write_file))
        read_file = open(r, mode="rb")
    file_pairs.append((read_file, outfile))
    for (infile, outfile), (operation, parameter) in zip(file_pairs, pipeline):
        operation.apply(parameter)(infile, outfile)
        infile.close()
        outfile.close()
    pass
