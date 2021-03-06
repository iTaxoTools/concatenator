# Concatenator
Performs a sequence of transformations on an input file.<br>
[*You may also be interested in the new GUI*](https://github.com/iTaxoTools/ConcatenatorGui)

### Installation
Clone and install the latest version (requires Python 3.8 or later):
```
git clone https://github.com/iTaxoTools/concatenator.git
cd concatenator
pip install . -f packages.html
```

### Executables
Download and run the standalone executables without installing Python.</br>
[See the latest release here.](https://github.com/iTaxoTools/concatenator/releases/latest)

## Usage
Launch the GUI:
```
concatenator
```

Construct a sequence of operations by clicking the buttons on the left side
(`Operation` panel), then click `Execute` button.

### Supported operations
* `Tabfile to NEXUS`: Converts a Tab file to a NEXUS file.
  - Input: Tab file with separate gene columns.
  - Output: Nexus file.
* `NEXUS to Tabfile`: Converts a NEXUS file to a Tab file.
  - Input: Nexus file.
  - Output: Tab file with separate gene columns.
* `Concatenate`: Concatenates all the columns with the name containing 'sequence' in a Tab file.
  - Input: Tab file with separate gene columns.
  - Output: Concatenated Tab file.
* `DNAconvert`: Convert a concatenated Tab file into a chosen format.
  - Input: Concatenated Tab file.
  - Output: Concatenated file of the chosen format.
* `Tabfile to multifile Fasta`: Convert a Tab file to a set of Fasta files for each gene.
  - Input: Tab file with separate gene columns
  - Output: Output archive with Fasta files
* `Tabfile to multifile Phylip`: Convert a Tab file to a set of Phylip files for each gene.
  - Input: Tab file with separate gene columns
  - Output: Output archive with Phylip files
* `Tabfile to multifile Ali`: Convert a Tab file to a set of Ali files for each gene.
  - Input: Tab file with separate gene columns
  - Output: Output archive with Ali files
* `Multifile Fasta to tabfile`: Convert an archive of Fasta files to a Tab file
  - Input: Input archive with Fasta files
  - Output: Tab file with separate gene columns
* `Multifile Ali to tabfile`: Convert an archive of Ali files to a Tab file
  - Input: Input archive with Ali files
  - Output: Tab file with separate gene columns
* `Multifile Phylip to tabfile`: Convert an archive of Phylip files to a Tab file
  - Input: Input archive with Phylip files
  - Output: Tab file with separate gene columns
* `Tabfile to Partitionfinder archive`: Convert a Tab file into the Partitionfinder input format
  - Input: Tab file with separate gene columns
  - Output: Archived Partitionfinder input
* `NEXUS to Partitionfinder archive`: Convert a NEXUS file into the Partitionfinder input format
  - Input: Nexus file
  - Output: Archived Partitionfinder input
* `Split codon positions in tabfile`: Separates genes in a Tabfile into codon position charsets
  - Input: Tab file with separate gene columns
  - Output: Tab file with separate codon position charsets

### Python API

The API exposes the new backend. See `scripts/convert.py` for an example.

### Packaging

It is advised to use PyInstaller within a virtual environment:
```
pip install ".[dev]" -f packages.html
pyinstaller scripts/concatenator.spec
```
