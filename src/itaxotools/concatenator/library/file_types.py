#!/usr/bin/env python3

from enum import Enum


class FileType(Enum):
    """Types of files that can be processed by operations"""

    TabFile = ("Tab file", ".tab", False)
    NexusFile = ("NEXUS file", ".nex", False)
    FastaFile = ("Fasta file", ".fas", False)
    PhylipFile = ("Phylip file", ".phy", False)
    ConcatTabFile = ("Concatenated Tab file", ".tab", False)
    ConcatFasta = ("Concatenated FASTA file", ".fas", False)
    ConcatPhylip = ("Concatenated Phylip file", ".phy", False)
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
