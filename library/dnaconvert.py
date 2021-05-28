#!/usr/bin/env python3

from typing import TextIO

import sys
import os
import tempfile
import subprocess

resource_path = getattr(sys, '_MEIPASS', sys.path[0])

try:
    dnaconvert_path = None
    with open(os.path.join(resource_path, "data", "options.tab")) as options:
        for line in options:
            option, _, value = line.rstrip().partition("\t")
            if option == "DNAconvert":
                dnaconvert_path = value
except FileNotFoundError:
    dnaconvert_path = None

def process(infile: TextIO, outfile: TextIO, format: str) -> None:
    if not dnaconvert_path:
        raise ValueError("DNAconvert path is not specified")
    _, dna_input = tempfile.mkstemp()
    _, dna_output = tempfile.mkstemp()
    dnaconvert_dir, _ = os.path.split(dnaconvert_path)
    with open(dna_input, mode="w") as dna_infile:
        for line in infile:
            dna_infile.write(line)
    try:
        subprocess.run([dnaconvert_path, "--cmd", "--informat", "tab", "--outformat", format, dna_input, dna_output], cwd=dnaconvert_dir)
    except FileNotFoundError as e:
        raise ValueError(f"File {dnaconvert_path} doesn't exist") from e
    except PermissionError as e:
        raise ValueError(f"File {dnaconvert_path} is not executable") from e
    else:
        with open(dna_output) as dna_outfile:
            for line in dna_outfile:
                outfile.write(line)
    os.unlink(dna_input)
    os.unlink(dna_output)
