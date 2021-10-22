#!/usr/bin/env python3

from typing import Iterator, List, Tuple, Callable
import random
import string

import pandas as pd

from itaxotools.concatenator import FileType, FileFormat

_TABLE_LEN = 20

_METADATA_COLUMNS = ("species", "specimen-voucher", "locality")

_SEQ_LEN = (9, 17)


def _gen_seq() -> str:
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
    gene_num = random.randint(3, 11)
    seq_generators = [("sequence_" + _gen_str(3)(), _gen_seq) for _ in range(gene_num)]
    return _gen_columns([("seqid", gen_seqid.__next__)] +
                        [(name, _gen_str(8)) for name in _METADATA_COLUMNS] +
                        seq_generators)


def stream_table(table: pd.DataFrame) -> Iterator[pd.Series]:
    seq_columns = [column for column in table.columns
                   if column != "seqid" and column not in _METADATA_COLUMNS]
    index = pd.MultiIndex.from_frame(table[["seqid"] + list(_METADATA_COLUMNS)])
    for column in seq_columns:
        yield pd.Series(table[column].to_list(), index=index, name=column)
