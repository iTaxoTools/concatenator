#!/usr/bin/env python3
from __future__ import annotations

from typing import Iterator, Tuple
from enum import Enum, auto
from dataclasses import dataclass

import pandas as pd
import networkx as nx

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

    @classmethod
    def empty(cls) -> GeneralInfo:
        return cls(
            pd.DataFrame(
                index=pd.MultiIndex.from_tuples(
                    [], names=[InfoColumns.Taxon, InfoColumns.Gene]
                ),
                columns=[
                    column
                    for column in InfoColumns
                    if column not in {InfoColumns.Taxon, InfoColumns.Gene}
                ],
            )
        )

    def __add__(self, other: GeneralInfo) -> GeneralInfo:
        left = self.dataframe.reset_index()
        right = other.dataframe.reset_index()
        result = pd.concat([left, right], ignore_index=True)
        result = result.groupby([InfoColumns.Gene, InfoColumns.Taxon]).agg(
            {
                InfoColumns.SeqCount: "sum",
                InfoColumns.NucleotideCount: "sum",
                InfoColumns.MissingCount: "sum",
                InfoColumns.SeqLenMin: "min",
                InfoColumns.SeqLenMax: "max",
                InfoColumns.GCCount: "sum",
            }
        )
        return GeneralInfo(result)

    def total_data(self) -> pd.Series:
        """
        Returns a pandas Series with index:
                Number of markers
                Number of samples
                Total number of nucleotides
                Missing nucleotides (%)
                GC content (%)
                Minimum number of nucleotides per sample
                Maximum number of nucleotides per sample
                Average number of nucleotides per sample
                Minimum number of markers per sample
                Maximum number of markers per sample
                Average number of markers per sample
                Minimum number of samples per marker
                Maximum number of samples per marker
                Average number of samples per marker
        """
        result = pd.Series(
            index=[
                "Number of markers",
                "Number of samples",
                "Total number of nucleotides",
                "Missing nucleotides (%)",
                "GC content (%)",
                "Minimum number of nucleotides per sample",
                "Maximum number of nucleotides per sample",
                "Average number of nucleotides per sample",
                "Minimum number of markers per sample",
                "Maximum number of markers per sample",
                "Average number of markers per sample",
                "Minimum number of samples per marker",
                "Maximum number of samples per marker",
                "Average number of samples per marker",
            ],
            dtype="object",
        )
        dataframe = self.dataframe.reset_index()
        result["Number of samples"] = dataframe[InfoColumns.Taxon].nunique()
        result["Number of markers"] = dataframe[InfoColumns.Gene].nunique()
        result["Missing nucleotides (%)"] = (
            dataframe[InfoColumns.MissingCount].sum()
            / (
                dataframe[InfoColumns.NucleotideCount].sum()
                + dataframe[InfoColumns.MissingCount].sum()
            )
            * 100
        )
        result["GC content (%)"] = (
            dataframe[InfoColumns.GCCount].sum()
            / (dataframe[InfoColumns.NucleotideCount].sum())
            * 100
        )
        result["Total number of nucleotides"] = dataframe[
            InfoColumns.NucleotideCount
        ].sum()
        result["Minimum number of nucleotides per sample"] = (
            dataframe.groupby(InfoColumns.Taxon)[InfoColumns.NucleotideCount]
            .sum()
            .min()
        )
        result["Maximum number of nucleotides per sample"] = (
            dataframe.groupby(InfoColumns.Taxon)[InfoColumns.NucleotideCount]
            .sum()
            .max()
        )
        result["Average number of nucleotides per sample"] = (
            result["Total number of nucleotides"]
            / result["Number of samples"]
        )
        markers_per_taxon = dataframe.groupby(InfoColumns.Taxon)[
            InfoColumns.Gene
        ].count()
        result["Minimum number of markers per sample"] = markers_per_taxon.min()
        result["Maximum number of markers per sample"] = markers_per_taxon.max()
        result["Average number of markers per sample"] = markers_per_taxon.mean()
        taxa_per_marker = dataframe.groupby(InfoColumns.Gene)[InfoColumns.Taxon].count()
        result["Minimum number of samples per marker"] = taxa_per_marker.min()
        result["Maximum number of samples per marker"] = taxa_per_marker.max()
        result["Average number of samples per marker"] = taxa_per_marker.mean()
        return result

    @staticmethod
    def by_input_file(tables: Iterator[FileGeneralInfo]) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with columns:
            Input file name
            Input file format
            Number of markers
            Number of samples
            Sequence length minimum
            Sequence length maximum
            Missing nucleotides (%)
            GC content (%)
        """
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
            row["Input file name"] = file_info.filename
            row["Input file format"] = file_info.file_format
            rows.append(row)
        result = pd.DataFrame(rows)
        result[InfoColumns.GCCount] = (
            result[InfoColumns.GCCount] / (result[InfoColumns.NucleotideCount]) * 100
        )
        result[InfoColumns.MissingCount] = (
            result[InfoColumns.MissingCount]
            / (result[InfoColumns.NucleotideCount] + result[InfoColumns.MissingCount])
            * 100
        )
        result.drop(columns=InfoColumns.NucleotideCount, inplace=True)
        result.rename(
            columns={
                InfoColumns.Taxon: "Number of samples",
                InfoColumns.Gene: "Number of markers",
                InfoColumns.SeqLenMin: "Sequence length minimum",
                InfoColumns.SeqLenMax: "Sequence length maximum",
                InfoColumns.MissingCount: "Missing nucleotides (%)",
                InfoColumns.GCCount: "GC content (%)",
            },
            inplace=True,
        )
        result.set_index("Input file name", inplace=True)
        return result[
            [
                "Input file format",
                "Number of markers",
                "Number of samples",
                "Sequence length minimum",
                "Sequence length maximum",
                "Missing nucleotides (%)",
                "GC content (%)",
            ]
        ]

    def by_taxon(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with index "taxon name" and columns:
            Markers with data
            Total number of nucleotides
            Sequence length minimum
            Sequence length maximum
            Missing nucleotides (%)
            Missing markers (%)
            GC content (%)
        """
        dataframe = self.dataframe.reset_index()
        result = dataframe.groupby(InfoColumns.Taxon).agg(
            **{
                "Markers with data": pd.NamedAgg(
                    column=InfoColumns.Gene, aggfunc="nunique"
                ),
                "Total number of nucleotides": pd.NamedAgg(
                    column=InfoColumns.NucleotideCount, aggfunc="sum"
                ),
                "Sequence length minimum": pd.NamedAgg(
                    column=InfoColumns.SeqLenMin, aggfunc="min"
                ),
                "Sequence length maximum": pd.NamedAgg(
                    column=InfoColumns.SeqLenMax, aggfunc="max"
                ),
                "GC content (%)": pd.NamedAgg(
                    column=InfoColumns.GCCount, aggfunc="sum"
                ),
            }
        )
        for_gene_max_length = dataframe[[InfoColumns.Gene]].copy()
        for_gene_max_length["length"] = (
            dataframe[InfoColumns.NucleotideCount] + dataframe[InfoColumns.MissingCount]
        )
        total_length = float(
            for_gene_max_length.groupby(InfoColumns.Gene)["length"].max().sum()
        )
        result["Missing markers (%)"] = (
            1
            - result["Markers with data"]
            / dataframe[InfoColumns.Gene].nunique()
        ) * 100
        result["Missing nucleotides (%)"] = (
            1 - result["Total number of nucleotides"] / total_length
        ) * 100
        result["GC content (%)"] *= 100 / result["Total number of nucleotides"]
        result["Missing nucleotides (%)"] = result["Missing nucleotides (%)"].astype("float")
        result["GC content (%)"] = result["GC content (%)"].astype("float")
        result.index.name = "taxon name"

        return result

    def by_gene(self, gene_info: GeneInfo) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with index "gene name" and columns:
            Number of samples with data
            Total number of nucleotides
            Missing nucleotides (%)
            Missing samples (%)
            GC content (%)
            Re-aligned by MAFFT
            Padded to compensate for unequal sequence lengths
            Padded to start with 1st codon position
        """
        dataframe = self.dataframe.reset_index()
        result = dataframe.groupby(InfoColumns.Gene).agg(
            **{
                "Number of samples with data": pd.NamedAgg(
                    column=InfoColumns.Taxon, aggfunc="nunique"
                ),
                "Total number of nucleotides": pd.NamedAgg(
                    column=InfoColumns.NucleotideCount, aggfunc="sum"
                ),
                "Missing nucleotides (%)": pd.NamedAgg(
                    column=InfoColumns.MissingCount, aggfunc="sum"
                ),
                "GC content (%)": pd.NamedAgg(
                    column=InfoColumns.GCCount, aggfunc="sum"
                ),
            }
        )
        result["Missing samples (%)"] = (
            1
            - result["Number of samples with data"]
            / dataframe[InfoColumns.Taxon].nunique()
        ) * 100
        result["Missing nucleotides (%)"] *= 100 / (
            result["Total number of nucleotides"]
            + result["Missing nucleotides (%)"]
        )
        result["GC content (%)"] *= (
            100 / result["Total number of nucleotides"]
        )
        result = result.join(gene_info.dataframe, how="left")
        for yes_no_column in GeneInfoColumns:
            result[yes_no_column] = result[yes_no_column].map(
                {True: "yes", False: "no"}
            )
        result.rename(
            columns={
                GeneInfoColumns.MafftRealigned: "Re-aligned by MAFFT",
                GeneInfoColumns.PaddedLength: "Padded to compensate for unequal sequence lengths",
                GeneInfoColumns.PaddedCodonPosition: "Padded to start with 1st codon position",
            },
            inplace=True,
        )
        result["Missing nucleotides (%)"] = result["Missing nucleotides (%)"].astype("float")
        result["GC content (%)"] = result["GC content (%)"].astype("float")
        result.index.name = "gene name"

        return result

    def disjoint_taxon_groups(self) -> Iterator[set]:
        """
        Yields sets of taxon, such taxons in different sets have no genes in common
        """
        left_taxa = self.dataframe.index.to_frame(index=False)
        right_taxa = left_taxa.copy()
        left_taxa.rename(columns={InfoColumns.Taxon: "left"}, inplace=True)
        right_taxa.rename(columns={InfoColumns.Taxon: "right"}, inplace=True)
        connected_taxa = left_taxa.merge(
            right_taxa, on=InfoColumns.Gene, suffixes=[None, None]
        )[["left", "right"]]
        graph = nx.from_pandas_edgelist(
            connected_taxa, source="left", target="right", edge_attr=None
        )
        for component in nx.connected_components(graph):
            yield component

    def unconnected_taxons(self) -> Iterator[Tuple[str, str]]:
        """
        Yields pairs of taxon names that have no genes in common
        """
        left_taxa = self.dataframe.index.to_frame(index=False)
        right_taxa = left_taxa.copy()
        left_taxa.rename(columns={InfoColumns.Taxon: "left"}, inplace=True)
        right_taxa.rename(columns={InfoColumns.Taxon: "right"}, inplace=True)
        connected_taxa = left_taxa.merge(
            right_taxa, on=InfoColumns.Gene, suffixes=[None, None]
        )[["left", "right"]]
        connected = pd.MultiIndex.from_frame(connected_taxa)
        all_pairs = pd.MultiIndex.from_product(
            [connected_taxa["left"].unique(), connected_taxa["left"].unique()]
        )
        unconnected = all_pairs.difference(connected)
        for taxon1, taxon2 in unconnected:
            if taxon1 < taxon2:
                yield (taxon1, taxon2)


@dataclass
class FileGeneralInfo:
    filename: str
    file_format: FileFormat
    table: GeneralInfo


class GeneInfoColumns(Enum):
    MafftRealigned = auto()  # Re-aligned by MAFFT
    PaddedLength = auto()  # Padded to compensate for unequal sequence lengths
    PaddedCodonPosition = auto()  # Padded to start with 1st codon position


class GeneInfo:
    """
    Contains a pandas DataFrame with index `InfoColumns.Gene` and columns:
        GeneInfoColumns.MafftRealigned (Re-aligned by MAFFT)
        GeneInfoColumns.PaddedLength (Padded to compensate for unequal sequence lengths)
        GeneInfoColumns.PaddedCodonPosition (Padded to start with 1st codon position)
    """

    def __init__(self, dataframe: pd.DataFrame):
        assert dataframe.index.name == InfoColumns.Gene
        assert set(dataframe.columns) == set(GeneInfoColumns)
        self.dataframe = dataframe
