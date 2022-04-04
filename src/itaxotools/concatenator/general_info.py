#!/usr/bin/env python3
from __future__ import annotations

from enum import Enum, auto

import pandas as pd
import numpy as np


class InfoColumns(Enum):
    Taxon = auto()
    Gene = auto()
    SeqCount = auto()
    NucleotideCount = auto()
    MissingCount = auto()
    SeqLenMin = auto()
    SeqLenMax = auto()
    GCCount = auto()


class GeneralInfo:
    def __init__(self, dataframe: pd.DataFrame):
        assert set(dataframe.index.names) == {InfoColumns.Taxon, InfoColumns.Gene}
        assert set(dataframe.columns) == set(InfoColumns) - {
            InfoColumns.Taxon,
            InfoColumns.Gene,
        }
        self.dataframe = dataframe

    def __add__(self, other: GeneralInfo) -> GeneralInfo:
        result = pd.DataFrame(index=self.dataframe.index.union(other.dataframe.index))
        for add_column in (
            InfoColumns.SeqCount,
            InfoColumns.NucleotideCount,
            InfoColumns.MissingCount,
            InfoColumns.SeqLenMin,
            InfoColumns.SeqLenMax,
            InfoColumns.GCCount,
        ):
            result[add_column] = self.dataframe[add_column].add(
                other.dataframe[add_column], fill_value=0
            )
        result[InfoColumns.SeqLenMin] = self.dataframe[InfoColumns.SeqLenMin].combine(
            other.dataframe[InfoColumns.SeqLenMin], min, fill_value=np.inf
        )
        result[InfoColumns.SeqLenMax] = self.dataframe[InfoColumns.SeqLenMax].combine(
            other.dataframe[InfoColumns.SeqLenMax], max, fill_value=0
        )
        return GeneralInfo(result)
