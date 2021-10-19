#!/usr/bin/env python3

from typing import Optional, Callable, TextIO, Any, Union, List, Tuple, BinaryIO
from enum import Enum
import os
import shutil

from . import concat
from . import tab_to_nexus
from . import nexus_to_tab
from . import dnaconvert
from . import tab_to_multifile
from . import multifile_to_tab
from . import tab_to_partition_finder
from . import nexus_to_partition_finder
from . import tab_to_codon_tab
from .file_utils import make_binary, make_binary_in, make_binary_out
from .file_types import FileFormat


Parameter = Union[None, int, str]


def dnaconvert_wrapper(format: str) -> Callable[[TextIO, TextIO], None]:
    def wrap(infile: TextIO, outfile: TextIO) -> None:
        dnaconvert.process(infile, outfile, format)

    return wrap


class Operation(Enum):
    """Operations that the program can perform on files"""

    TabToNexus = (FileFormat.TabFile, "Tabfile to NEXUS", None, None)
    NexusToTab = (FileFormat.NexusFile, "NEXUS to Tabfile", None, None)
    ConcatTab = (FileFormat.TabFile, "Concatenate", None, None)
    DnaConvert = (
        FileFormat.ConcatTabFile,
        "DNAconvert\n(concatenated)",
        "Convert into {} with DNAconvert",
        ["fasta", "phylip"],
    )
    TabToMFasta = (FileFormat.TabFile, "Tabfile to multifile Fasta", None, None)
    TabToMPhylip = (FileFormat.TabFile, "Tabfile to multifile Phylip", None, None)
    TabToMAli = (FileFormat.TabFile, "Tabfile to multifile Ali", None, None)
    MFastaToTab = (FileFormat.MultiFastaInput, "Multifile Fasta to tabfile", None, None)
    MAliToTab = (FileFormat.MultiAliInput, "Multifile Ali to tabfile", None, None)
    MPhylipToTab = (
        FileFormat.MultiPhylipInput,
        "Multifile Phylip to tabfile",
        None,
        None,
    )
    TabToPartition = (
        FileFormat.TabFile,
        "Tabfile to Partitionfinder archive",
        None,
        None,
    )
    NexusToPartition = (
        FileFormat.NexusFile,
        "NEXUS file to Partitionfinder archive",
        None,
        None,
    )
    TabToCodons = (FileFormat.TabFile, "Split codon positions in tabfile", None, None)

    def __init__(
        self,
        input_type: FileFormat,
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

    def output_type(self, parameter: Parameter) -> FileFormat:
        type = {
            Operation.TabToNexus: FileFormat.NexusFile,
            Operation.NexusToTab: FileFormat.TabFile,
            Operation.ConcatTab: FileFormat.ConcatTabFile,
            Operation.DnaConvert: FileFormat.ConcatFasta
            if parameter == "fasta"
            else FileFormat.ConcatPhylip
            if parameter == "phylip"
            else None,
            Operation.TabToMFasta: FileFormat.MultiFastaOutput,
            Operation.TabToMPhylip: FileFormat.MultiPhylipOutput,
            Operation.TabToMAli: FileFormat.MultiAliOutput,
            Operation.MFastaToTab: FileFormat.TabFile,
            Operation.MAliToTab: FileFormat.TabFile,
            Operation.MPhylipToTab: FileFormat.TabFile,
            Operation.TabToPartition: FileFormat.PartitionFinderOutput,
            Operation.NexusToPartition: FileFormat.PartitionFinderOutput,
            Operation.TabToCodons: FileFormat.CodonTab,
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
            Operation.TabToPartition: make_binary_in(tab_to_partition_finder.process),
            Operation.NexusToPartition: make_binary_in(
                nexus_to_partition_finder.process
            ),
            Operation.TabToCodons: make_binary(tab_to_codon_tab.process),
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
