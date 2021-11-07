
from typing import List

from .model import PathLike
from .types import Charset


def write_nexus(path: PathLike, charsets: List[Charset]) -> None:

    with path.open('w') as file:
        file.write('#NEXUS\n\n')
        file.write('BEGIN SETS;\n\n')

        for charset in charsets:
            file.write(f'charset {str(charset)}\n')
        file.write('\n')

        for charset in (cs for cs in charsets if cs.frame):
            for codon in charset.codon_sets():
                file.write(f'charset {str(codon)}\n')
            file.write('\n')

        file.write('END;\n')
