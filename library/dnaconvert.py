#!/usr/bin/env python3

from typing import TextIO

import sys
import os
import tempfile
import subprocess

resource_path = getattr(sys, '_MEIPASS', sys.path[0])

try:
    dnaconvert_path = None
    with open(os.path.join(resource_path, "options.tab")) as options:
        for line in options:
            option, _, value = line.rstrip().partition("\t")
            if option == "DNAconvert":
                dnaconvert_path = value
except FileNotFoundError:
    dnaconvert_path = None

def process(infile: TextIO, outfile: TextIO, format: str) -> None:
    if not dnaconvert_path:
        raise ValueError("DNAconvert path is not specified")
    dna_input = tempfile.NamedTemporaryFile(mode="w")
    dna_output = tempfile.NamedTemporaryFile(mode="r")
    for line in infile:
        dna_input.write(line)
    try:
        subprocess.run(["dnaconvert_path", "--cmd", "--informat", "tab", "--outformat", format, dna_input.name, dna_output.name])
    except FileNotFoundError:
        raise ValueError(f"File {dnaconvert_path} doesn't exist")
    except PermissionError:
        raise ValueError(f"File {dnaconvert_path} is not executable")
    else:
        for line in dna_output:
            outfile.write(line)
    dna_input.close()
    dna_output.close()
    os.unlink(dna_input.name)
    os.unlink(dna_output.name)
