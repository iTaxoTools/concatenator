#!/usr/bin/env python3

from enum import Enum, auto


class FileFormat(Enum):
    """Supported file formats supported by oprations"""

    TabFile = ("Tab file", ".tab", False)
    AliFile = ("Ali file", ".ali", False)
    NexusFile = ("NEXUS file", ".nex", False)
    FastaFile = ("Fasta file", ".fas", False)
    PhylipFile = ("Phylip file", ".phy", False)
    ConcatTabFile = ("Concatenated Tab file", ".tab", False)
    ConcatFasta = ("Concatenated FASTA file", ".fas", False)
    ConcatPhylip = ("Concatenated Phylip file", ".phy", False)
    MultiFastaDirectory = ("Multifile FASTA directory", "", True)
    MultiFastaArchive = ("Multifile FASTA archive", ".zip", True)
    MultiFastaOutput = ("Multifile FASTA output archive", ".zip", True)
    MultiPhylipOutput = ("Multifile Phylip output archive", ".zip", True)
    MultiAliOutput = ("Multifile Ali output archive", ".zip", True)
    MultiFastaInput = ("Multifile FASTA input archive", ".zip", True)
    MultiPhylipInput = ("Multifile Phylip input archive", ".zip", True)
    MultiAliInput = ("Multifile Ali input archive", ".zip", True)
    PartitionFinderOutput = ("Partitionfinder archive", ".zip", True)
    CodonTab = ("Codon positions tab file", ".tab", False)

    def __init__(self, description: str, extension: str, timestamp: bool):
        self.description = description
        self.extension = extension
        self.timestamp = timestamp


class FileType(Enum):
    """File type, describing how a file format is stored"""

    SingleFile = auto()
    MultiFileZip = auto()
    MultiFileDir = auto()
