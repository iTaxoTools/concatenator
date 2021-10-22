#!/usr/bin/env python3

from typing import Iterator, List, Tuple, Callable, Optional
import random
import string

import pandas as pd

from itaxotools.concatenator import FileType, FileFormat

_TABLE_LEN = 20

_METADATA_COLUMNS = ("species", "specimen-voucher", "locality")

_SEQ_LEN = (9, 17)


def only_seqid(filetype: FileType, format: FileFormat) -> bool:
    return ((filetype == FileType.File and format == FileFormat.Nexus)
            or format in [FileFormat.Ali, FileFormat.Phylip, FileFormat.Fasta])


def single_sequence(filetype: FileType, format: FileFormat) -> bool:
    return (filetype == FileType.File and
            format in [FileFormat.Ali, FileFormat.Phylip, FileFormat.Fasta])


def uniform_sequence_len(filetype: FileType, format: FileFormat) -> bool:
    return format in [FileFormat.Nexus, FileFormat.Phylip, FileFormat.Ali]


def _gen_seq(length: Optional[int] = None) -> str:
    if length:
        seq_len = length
    else:
        seq_len = random.randint(*_SEQ_LEN)
    return "".join(random.choices("ACGT", k=seq_len))


def _gen_str(length: int) -> Callable[[], str]:
    return lambda: "".join(random.choices(string.ascii_lowercase, k=length))


def _gen_columns(generators: List[Tuple[str, Callable[[], str]]]) -> pd.DataFrame:
    columns = pd.Index(map(lambda p: p[0], generators))
    row_gen = ([gen() for _, gen in generators] for _ in range(_TABLE_LEN))
    return pd.DataFrame(row_gen, columns=columns)


def gen_table(filetype: FileType, format: FileFormat) -> pd.DataFrame:
    gen_seqid = (str(i) for i in range(_TABLE_LEN))
    if only_seqid(filetype, format):
        metadata_generators = [("seqid", gen_seqid.__next__)]
    else:
        metadata_generators = ([("seqid", gen_seqid.__next__)] +
                               [(name, _gen_str(8)) for name in _METADATA_COLUMNS])
    if uniform_sequence_len(filetype, format):
        seq_len = random.randint(*_SEQ_LEN)

        def gen_seq() -> str:
            return _gen_seq(seq_len)
    else:
        def gen_seq() -> str:
            return _gen_seq()
    if single_sequence(filetype, format):
        seq_generators = [("sequence_" + _gen_str(3)(), gen_seq)]
    else:
        gene_num = random.randint(3, 17)
        seq_generators = [("sequence_" + _gen_str(3)(), gen_seq)
                          for _ in range(gene_num)]
    return _gen_columns(metadata_generators +
                        seq_generators)


def stream_table(table: pd.DataFrame) -> Iterator[pd.Series]:
    seq_columns = [column for column in table.columns
                   if column != "seqid" and column not in _METADATA_COLUMNS]
    index_columns = [column for column in table.columns if column not in seq_columns]
    index = pd.MultiIndex.from_frame(table[index_columns])
    for column in seq_columns:
        yield pd.Series(table[column].to_list(), index=index, name=column)
