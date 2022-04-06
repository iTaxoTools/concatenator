#!/usr/bin/env python3
from __future__ import annotations

from typing import Iterator
from enum import Enum, auto
from dataclasses import dataclass

import pandas as pd
import numpy as np

from .file_types import FileFormat


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

    def total_data(self) -> pd.Series:
        result = pd.Series(
            index=[
                "Number of taxa",
                "Number of genes (markers)",
                "% missing data",
                "GC content of sequences",
                "Total number of nucleotides in concatenated data set",
                "Minimum number of nucleotides per taxon",
                "Maximum number of nucleotides per taxon",
                "Average number of nucleotides per taxon",
                "Minimum number of markers per taxon",
                "Maximum number of markers per taxon",
                "Average number of markers per taxon",
                "Minimum number of taxa per markers",
                "Maximum number of taxa per marker",
                "Average number of taxa per marker",
            ]
        )
        dataframe = self.dataframe.reset_index()
        result["Number of taxa"] = dataframe[InfoColumns.Taxon].nunique()
        result["Number of genes (markers)"] = dataframe[InfoColumns.Gene].nunique()
        result["% missing data"] = (
            dataframe[InfoColumns.MissingCount].sum()
            / dataframe[InfoColumns.NucleotideCount].sum()
            * 100
        )
        result["GC content of sequences"] = (
            dataframe[InfoColumns.GCCount].sum()
            / dataframe[InfoColumns.NucleotideCount].sum()
        )
        result["Total number of nucleotides in concatenated data set"] = dataframe[
            InfoColumns.NucleotideCount
        ].sum()
        result["Minimum number of nucleotides per taxon"] = (
            dataframe.groupby(InfoColumns.Taxon)[InfoColumns.NucleotideCount]
            .sum()
            .min()
        )
        result["Maximum number of nucleotides per taxon"] = (
            dataframe.groupby(InfoColumns.Taxon)[InfoColumns.NucleotideCount]
            .sum()
            .max()
        )
        result["Average number of nucleotides per taxon"] = (
            result["Total number of nucleotides in concatenated data set"]
            / result["Number of taxa"]
        )
        markers_per_taxon = dataframe.groupby(InfoColumns.Taxon)[
            InfoColumns.Gene
        ].count()
        result["Minimum number of markers per taxon"] = markers_per_taxon.min()
        result["Maximum number of markers per taxon"] = markers_per_taxon.max()
        result["Average number of markers per taxon"] = markers_per_taxon.mean()
        taxa_per_marker = dataframe.groupby(InfoColumns.Gene)[InfoColumns.Taxon].count()
        result["Minimum number of taxa per markers"] = taxa_per_marker.min()
        result["Maximum number of taxa per marker"] = taxa_per_marker.max()
        result["Average number of taxa per marker"] = taxa_per_marker.mean()
        return result

    @staticmethod
    def by_input_file(tables: Iterator[FileGeneralInfo]) -> pd.DataFrame:
        rows = []
        for file_info in tables:
            row = file_info.table.dataframe.reset_index().agg(
                {
                    InfoColumns.Taxon: "nunique",
                    InfoColumns.Gene: "nunique",
                    InfoColumns.MissingCount: "sum",
                    InfoColumns.NucleotideCount: "sum",
                    InfoColumns.SeqLenMin: "min",
                    InfoColumns.SeqLenMax: "max",
                    InfoColumns.GCCount: "sum",
                }
            )
            row["input file name"] = file_info.filename
            row["input file format"] = file_info.file_format
            rows.append(row)
        result = pd.DataFrame(rows)
        result[InfoColumns.MissingCount] = (
            result[InfoColumns.MissingCount] / result[InfoColumns.NucleotideCount] * 100
        )
        result[InfoColumns.GCCount] = (
            result[InfoColumns.GCCount] / result[InfoColumns.NucleotideCount]
        )
        result.drop(columns=InfoColumns.NucleotideCount, inplace=True)
        result.rename(
            columns={
                InfoColumns.Taxon: "number of taxa",
                InfoColumns.Gene: "number of genes (markers)",
                InfoColumns.MissingCount: "% missing data",
                InfoColumns.SeqLenMin: "sequence length minimum",
                InfoColumns.SeqLenMax: "sequence length maximum",
                InfoColumns.GCCount: "GC content of sequences",
            }
        )
        return result


@dataclass
class FileGeneralInfo:
    filename: str
    file_format: FileFormat
    table: GeneralInfo