#!/usr/bin/env python

import pandas as pd

def into_seqids(table: pd.DataFrame) -> pd.Series:
    for column in table.columns:
        table[column] = table[column].str.replace(r'^[^A-Za-z09]+', '', n=1)
        table[column] = table[column].str.replace(r'[^A-Za-z0-9]+$', '', n=1)
        table[column] = table[column].str.replace(r'[^A-Za-z0-9]+', '_')
    return table.iloc[:, 0].str.cat(table.iloc[:,1:], sep='_')

