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
import library.multifile_to_tab as multifile_to_tab
from library.file_utils import make_binary, make_binary_in, make_binary_out


class FileType(Enum):
    """Types of files that can be processed by operations"""

    TabFile = ("Tab file", ".tab", False)
    NexusFile = ("NEXUS file", ".nex", False)
    ConcatTabFile = ("Concatenated Tab file", ".tab", False)
    ConcatFasta = ("Concatenated FASTA file", ".fas", False)
    ConcatPhylip = ("Concatenated Phylip file", ".phy", False)
    MultiFastaOutput = ("Multifile FASTA output archive", ".zip", True)
    MultiPhylipOutput = ("Multifile Phylip output archive", ".zip", True)
    MultiAliOutput = ("Multifile Ali output archive", ".zip", True)
    MultiFastaInput = ("Multifile FASTA input archive", ".zip", True)
    MultiPhylipInput = ("Multifile Phylip input archive", ".zip", True)
    MultiAliInput = ("Multifile Ali input archive", ".zip", True)

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
    MFastaToTab = (FileType.MultiFastaInput, "Multifile Fasta to tabfile", None, None)
    MAliToTab = (FileType.MultiAliInput, "Multifile Ali to tabfile", None, None)
    MPhylipToTab = (
        FileType.MultiPhylipInput,
        "Multifile Phylip to tabfile",
        None,
        None,
    )

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
        type = {
            Operation.TabToNexus: FileType.NexusFile,
            Operation.NexusToTab: FileType.TabFile,
            Operation.ConcatTab: FileType.ConcatTabFile,
            Operation.DnaConvert: FileType.ConcatFasta
            if parameter == "fasta"
            else FileType.ConcatPhylip
            if parameter == "phylip"
            else None,
            Operation.TabToMFasta: FileType.MultiFastaOutput,
            Operation.TabToMPhylip: FileType.MultiPhylipOutput,
            Operation.TabToMAli: FileType.MultiAliOutput,
            Operation.MFastaToTab: FileType.TabFile,
            Operation.MAliToTab: FileType.TabFile,
            Operation.MPhylipToTab: FileType.TabFile,
        }.get(self)
        assert type is not None
        return type

    def apply(self, parameter: Parameter) -> Callable[[BinaryIO, BinaryIO], None]:
        operation: Optional[Callable[[BinaryIO, BinaryIO], None]] = {
            Operation.TabToNexus: make_binary(tab_to_nexus.process),
            Operation.NexusToTab: make_binary(nexus_to_tab.process),
            Operation.ConcatTab: make_binary(concat.process),
            Operation.DnaConvert: make_binary(dnaconvert_wrapper(parameter))
            if type(parameter) is str
            else None,
            Operation.TabToMFasta: make_binary_in(tab_to_multifile.process_fasta),
            Operation.TabToMPhylip: make_binary_in(tab_to_multifile.process_phylip),
            Operation.TabToMAli: make_binary_in(tab_to_multifile.process_ali),
            Operation.MFastaToTab: make_binary_out(multifile_to_tab.process_fasta),
            Operation.MAliToTab: make_binary_out(multifile_to_tab.process_ali),
            Operation.MPhylipToTab: make_binary_out(multifile_to_tab.process_phylip),
        }.get(self)
        assert operation is not None
        return operation


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
