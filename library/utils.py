#!/usr/bin/env python

from typing import Optional

import pandas as pd


def into_seqids(table: pd.DataFrame) -> pd.Series:
    """
    Concatenates a collection of string columns into a column of seqids.

    Removes unallowed characters.
    """
    for column in table.columns:
        # Remove unallowed character in the beginning
        table[column] = table[column].str.replace(r"^[^A-Za-z09]+", "", n=1)
        # Remove unallowed character in the end
        table[column] = table[column].str.replace(r"[^A-Za-z0-9]+$", "", n=1)
        # Replace substrings of unallowed characters in the middle with _
        table[column] = table[column].str.replace(r"[^A-Za-z0-9]+", "_")
    return table.iloc[:, 0].str.cat(table.iloc[:, 1:], sep="_")


def max_length_if_not_same(column: pd.Series) -> Optional[int]:
    """
    Returns length of the longest string in the ``column``, if the strings have different length
    """
    max_length = column.str.len().max()
    min_length = column.str.len().min()
    if max_length > min_length:
        return max_length
    else:
        return None


def make_equal_length(
    column: pd.Series, max_length: Optional[int] = None, fillchar: str = "N"
) -> pd.Series:
    """
    Pads strings in the ``column`` to ``max_length``.

    If ``max_length`` is None, uses the length of the longest string
    """
    if not max_length:
        max_length = column.str.len().max()
    if max_length:
        return column.str.ljust(max_length, fillchar)
    else:
        return column
