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
                Number of taxa
                Number of genes (markers)
                % missing data
                GC content of sequences
                Total number of nucleotides in concatenated data set
                Minimum number of nucleotides per taxon
                Maximum number of nucleotides per taxon
                Average number of nucleotides per taxon
                Minimum number of markers per taxon
                Maximum number of markers per taxon
                Average number of markers per taxon
                Minimum number of taxa per markers
                Maximum number of taxa per marker
                Average number of taxa per marker
        """
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
            ],
            dtype="object",
        )
        dataframe = self.dataframe.reset_index()
        result["Number of taxa"] = dataframe[InfoColumns.Taxon].nunique()
        result["Number of genes (markers)"] = dataframe[InfoColumns.Gene].nunique()
        result["% missing data"] = (
            dataframe[InfoColumns.MissingCount].sum()
            / (
                dataframe[InfoColumns.NucleotideCount].sum()
                + dataframe[InfoColumns.MissingCount].sum()
            )
            * 100
        )
        result["GC content of sequences"] = (
            dataframe[InfoColumns.GCCount].sum()
            / (dataframe[InfoColumns.NucleotideCount].sum())
            * 100
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
        """
        Returns a pandas DataFrame with columns:
            input file name
            input file format
            number of samples
            number of markers
            sequence length minimum
            sequence length maximum
            % missing nucleotides
            % GC content
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
            row["input file name"] = file_info.filename
            row["input file format"] = file_info.file_format
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
                InfoColumns.Taxon: "number of samples",
                InfoColumns.Gene: "number of markers",
                InfoColumns.SeqLenMin: "sequence length minimum",
                InfoColumns.SeqLenMax: "sequence length maximum",
                InfoColumns.MissingCount: "% missing nucleotides",
                InfoColumns.GCCount: "% GC content",
            },
            inplace=True,
        )
        return result[
            [
                "input file name",
                "input file format",
                "number of samples",
                "number of markers",
                "sequence length minimum",
                "sequence length maximum",
                "% missing nucleotides",
                "% GC content",
            ]
        ]

    def by_taxon(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with index "taxon name" and columns:
            number of markers with data
            total number of nucleotides
            % of missing data (nucleotides)
            % of missing data (markers)
            sequence length minimum
            sequence length maximum
            GC content of sequences
        """
        dataframe = self.dataframe.reset_index()
        result = dataframe.groupby(InfoColumns.Taxon).agg(
            **{
                "number of markers with data": pd.NamedAgg(
                    column=InfoColumns.Gene, aggfunc="nunique"
                ),
                "total number of nucleotides": pd.NamedAgg(
                    column=InfoColumns.NucleotideCount, aggfunc="sum"
                ),
                "sequence length minimum": pd.NamedAgg(
                    column=InfoColumns.SeqLenMin, aggfunc="min"
                ),
                "sequence length maximum": pd.NamedAgg(
                    column=InfoColumns.SeqLenMax, aggfunc="max"
                ),
                "GC content of sequences": pd.NamedAgg(
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
        result["% of missing data (markers)"] = (
            1
            - result["number of markers with data"]
            / dataframe[InfoColumns.Gene].nunique()
        ) * 100
        result["% of missing data (nucleotides)"] = (
            1 - result["total number of nucleotides"] / total_length
        ) * 100
        result["GC content of sequences"] *= 100 / result["total number of nucleotides"]
        result["% of missing data (nucleotides)"] = result["% of missing data (nucleotides)"].astype("float")
        result["GC content of sequences"] = result["GC content of sequences"].astype("float")
        result.index.name = "taxon name"

        return result

    def by_gene(self, gene_info: GeneInfo) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with index "gene name" and columns:
            number of taxa with data
            total number of nucleotide in alignment
            % of missing data (nucleotides)
            % of missing data (taxa)
            GC content of sequences
            re-aligned by Mafft yes/no
            padded to compensate for unequal sequence lengths yes/no
            padded to start with 1st codon position yes/no
        """
        dataframe = self.dataframe.reset_index()
        result = dataframe.groupby(InfoColumns.Gene).agg(
            **{
                "number of taxa with data": pd.NamedAgg(
                    column=InfoColumns.Taxon, aggfunc="nunique"
                ),
                "total number of nucleotides in alignment": pd.NamedAgg(
                    column=InfoColumns.NucleotideCount, aggfunc="sum"
                ),
                "% of missing data (nucleotides)": pd.NamedAgg(
                    column=InfoColumns.MissingCount, aggfunc="sum"
                ),
                "GC content of sequences": pd.NamedAgg(
                    column=InfoColumns.GCCount, aggfunc="sum"
                ),
            }
        )
        result["% of missing data (taxa)"] = (
            1
            - result["number of taxa with data"]
            / dataframe[InfoColumns.Taxon].nunique()
        ) * 100
        result["% of missing data (nucleotides)"] *= 100 / (
            result["total number of nucleotides in alignment"]
            + result["% of missing data (nucleotides)"]
        )
        result["GC content of sequences"] *= (
            100 / result["total number of nucleotides in alignment"]
        )
        result = result.join(gene_info.dataframe, how="left")
        for yes_no_column in GeneInfoColumns:
            result[yes_no_column] = result[yes_no_column].map(
                {True: "yes", False: "no"}
            )
        result.rename(
            columns={
                GeneInfoColumns.MafftRealigned: "re-aligned by Mafft yes/no",
                GeneInfoColumns.PaddedLength: "padded to compensate for unequal sequence lengths yes/no",
                GeneInfoColumns.PaddedCodonPosition: "padded to start with 1st codon position yes/no",
            },
            inplace=True,
        )
        result["% of missing data (nucleotides)"] = result["% of missing data (nucleotides)"].astype("float")
        result["GC content of sequences"] = result["GC content of sequences"].astype("float")
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
    MafftRealigned = auto()  # re-aligned by Mafft yes/no
    PaddedLength = auto()  # padded to compensate for unequal sequence lengths yes/no
    PaddedCodonPosition = auto()  # padded to start with 1st codon position yes/no


class GeneInfo:
    """
    Contains a pandas DataFrame with index `InfoColumns.Gene` and columns:
        GeneInfoColumns.MafftRealigned (re-aligned by Mafft yes/no)
        GeneInfoColumns.PaddedLength (padded to compensate for unequal sequence lengths yes/no)
        GeneInfoColumns.PaddedCodonPosition (padded to start with 1st codon position yes/no)
    """

    def __init__(self, dataframe: pd.DataFrame):
        assert dataframe.index.name == InfoColumns.Gene
        assert set(dataframe.columns) == set(GeneInfoColumns)
        self.dataframe = dataframe
