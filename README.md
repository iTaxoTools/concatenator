# Concatenator
Performs a sequence of transformations on an input file.

# Usage
Construct a sequence of operations by clicking the buttons on the left side (`Operation` panel), then click `Execute` button.

# Supported operations
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

# Requirements
The file `data\options.tab` should contain a valid path to DNAconvert executable, otherwise `DNAconvert` operation will not work.
