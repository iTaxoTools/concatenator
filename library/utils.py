#!/usr/bin/env python

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
