#!/usr/bin/env python3

from typing import Optional, Callable, TextIO, Any, Union
from enum import Enum

import library.concat as concat
import library.tab_to_nexus as tab_to_nexus
import library.nexus_to_tab as nexus_to_tab

class FileType(Enum):
    """Types of files that can be processed by operations"""
    TabFile = ("Tab file", ".tab")
    NexusFile = ("NEXUS file", ".nex")
    ConcatTabFile = ("Concatenated Tab file", ".tab")

    def __init__(self, description: str, extension: str):
        self.description = description
        self.extension = extension

Parameter = Union[None, int, str]

class Operation(Enum):
    """Operations that the program can perform on files"""
    TabToNexus = (FileType.TabFile, "Tabfile to NEXUS", None, None)
    NexusToTab = (FileType.NexusFile, "NEXUS to Tabfile", None, None)
    ConcatTab = (FileType.TabFile, "Concatenate", None, None)

    def __init__(self, input_type: FileType, button_text: str, description: Optional[str], parameter_type: Any):
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
        else:
            assert False
    
    def apply(self, parameter: Parameter) -> Callable[[TextIO, TextIO], None]:
        if self == Operation.TabToNexus:
            return tab_to_nexus.process
        elif self == Operation.NexusToTab:
            return nexus_to_tab.process
        elif self == Operation.ConcatTab:
            return concat.process
        else:
            assert False
