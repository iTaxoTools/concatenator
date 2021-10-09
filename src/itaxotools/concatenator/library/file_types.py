#!/usr/bin/env python3

from enum import Enum
from typing import Optional


class FileFormat(Enum):
    """File formats supported by operations"""

    Tab = ("Tab-separated Table", ".tab", False)
    Nexus = ("Interleaved NEXUS", ".nex", False)
    Ali = ("Ali", ".ali", False)
    Fasta = ("Fasta", ".fas", False)
    Phylip = ("Phylip", ".phy", False)

    # Pending cleanup

    TabFile = ("Tab file", ".tab", False)
    AliFile = ("Ali file", ".ali", False)
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


class FileType(Enum):
    """Containers for file formats"""

    File = ('Single File', None)
    Directory = ('Directory', '')
    ZipArchive = ('Zip Archive', '.zip')

    def __init__(self, description: str, extension: Optional[str]):
        self.description = description
        self.extension = extension
