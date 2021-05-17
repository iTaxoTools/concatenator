#!/usr/bin/env python

from typing import TextIO, Iterator
import tempfile
from collections import defaultdict

import pandas as pd

from library.utils import *

def write(genes: Iterator[pd.DataFrame], output: TextIO) -> None:
    buf = tempfile.TemporaryFile(mode="w+")
    charsets = {}
    ntax = 0

    for gene in genes:
        output_fragment = pd.DataFrame()
        output_fragment['seqid'] = into_seqids(gene.iloc[:, :-1].copy())
        output_fragment['sequence'] = gene.iloc[:, -1]

        gene_name = gene.columns[-1][len("sequence_"):]
        gene_len = len(gene.iat[0, -1])
        charsets[gene_name] = gene_len

        output_fragment.to_string(buf, header=False, index=False)
        print('\n', file=buf)

        ntax = len(output_fragment)

    buf.seek(0)

    nchar = sum(charsets.values())

    print('#NEXUS\n', file=output)
    print('begin data;\n', file=output)
    print('format datatype=DNA missing=N missing=? Gap=- Interleave=yes;\n', file=output)
    print(f'dimensions Nchar={nchar} Ntax={ntax}\n', file=output)
    print('matrix\n', file=output)

    for line in buf:
        output.write(line)
    print(';\nend;', file=output)

    output.write('\n\n\n')

    print('begin sets;\n', file=output) 

    gene_position = 1

    for gene_name, gene_len in charsets.items():
        print('charset ', gene_name, ' = ', gene_position, '-', gene_position + gene_len - 1, ';', sep='', file=output)
        gene_position += gene_len
    
    print('\nend;', file=output)
    buf.close()
