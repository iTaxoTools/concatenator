
from typing import List

from .model import PathLike
from .types import Charset


def write_nexus(
    path: PathLike,
    charsets: List[Charset],
    include_full_markers_with_codons: bool = False,
) -> None:

    with path.open('w', encoding='utf-8', errors='surrogateescape') as out:
        out.write('#NEXUS\n\n')
        out.write('BEGIN SETS;\n\n')

        if include_full_markers_with_codons:
            for charset in charsets:
                out.write(f'charset {str(charset)}\n')
            out.write('\n')
            for charset in (cs for cs in charsets if cs.frame):
                for codon in charset.codon_sets():
                    out.write(f'charset {str(codon)}\n')
                out.write('\n')
        else:
            full_charsets = [cs for cs in charsets if not cs.frame]
            split_charsets = [cs for cs in charsets if cs.frame]

            if full_charsets:
                for charset in full_charsets:
                    out.write(f'charset {str(charset)}\n')
                out.write('\n')

            if split_charsets:
                for charset in split_charsets:
                    for codon in charset.codon_sets():
                        out.write(f'charset {str(codon)}\n')
                    out.write('\n')

        out.write('END;\n')
