#!/usr/bin/env python3

from enum import Enum, auto, Union

class FileType(Enum):
    """Types of files that can be processed by operations"""
    TabFile = auto()
    NexusFile = auto()    
    ConcatTabFile = auto()

Parameter = Union[None, int, str]

class Operation(Enum):
    """Operations that the program can perform on files"""
    TabToNexus = (FileType.TabFile, "Tabfile to NEXUS", )
    NexusToTab = (FileType.NexusFile, "NEXUS to Tabfile")
    ConcatTab = (FileType.TabFile, "Concatenate")

    def __init__(self, input_type: FileType, description: str):
       self.input_type = input_type
       self.description = description
    
