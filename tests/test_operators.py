from pathlib import Path
from typing import Iterator
import random

import pytest
import pandas as pd
import numpy as np

from itaxotools.concatenator.library.model import GeneSeries, GeneStream, Operator
from itaxotools.concatenator.library.types import TextCase, Charset
from itaxotools.concatenator.library.codons import GeneticCode, ReadingFrame
from itaxotools.concatenator.library.operators import (
    OpSanitizeGeneNames,
    OpSanitizeSpeciesNames,
    OpSequenceCase,
    OpTranslateMissing,
    OpTranslateGap,
    OpSpreadsheetCompatibility,
    OpFilterGenes,
    OpStencilGenes,
    OpBlock,
    OpReverseComplement,
    OpReverseNegativeReadingFrames,
    OpPadReadingFrames,
    OpDropEmpty,
    OpExtractCharsets,
    OpUpdateMetadata,
    OpMakeUniform,
    OpGeneralInfo,
    OpIndexMerge,
    OpGeneralInfoPerFile,
    OpGeneralInfoPerGene,
    OpTagSet,
    OpTagDelete,
    OpDropIfAllEmpty,
)
from itaxotools.concatenator.library.file_readers import read_from_path
from itaxotools.concatenator.library.general_info import (
    FileGeneralInfo,
    GeneralInfo,
    InfoColumns,
    GeneInfo,
    GeneInfoColumns,
)
from itaxotools.concatenator.library.file_types import FileFormat

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


def assert_gene_meta_equal(gene1, gene2):
    assert gene1.missing == gene2.missing
    assert gene1.gap == gene2.gap
    assert gene1.genetic_code == gene2.genetic_code
    assert gene1.reading_frame == gene2.reading_frame
    assert gene1.codon_names == gene2.codon_names


def test_tag_set_delete():
    series = pd.Series(
        {"seq1": "ACGT"},
        name="gene",
    )
    series.index.name = "seqid"
    gene = GeneSeries(series)
    added = OpTagSet("test", 42)(gene)
    assert added.tags.test == 42
    removed = OpTagDelete("test")(added)
    assert "test" not in removed.tags


@pytest.fixture
def gene_unsanitized() -> GeneSeries:
    series = pd.Series(
        {
            "ΔšëqΔ¹Δ": "GCAGTATAA",
        },
        name="ΔgëneΔ²Δ",
    )
    series.index.name = "seqid"
    return GeneSeries(series)


def test_sanitize_genes(gene_unsanitized):
    gene = gene_unsanitized
    altered = OpSanitizeGeneNames()(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.name == "gene_2"


def test_sanitize_species(gene_unsanitized):
    gene = gene_unsanitized
    altered = OpSanitizeSpeciesNames()(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.index[0] == "seq_1"


@pytest.fixture
def gene_case_mixed() -> GeneSeries:
    series = pd.Series(
        {
            "seq1": "gcaGTATAA",
        },
        name="gene",
    )
    series.index.name = "seqid"
    return GeneSeries(series)


def test_case_unchanged(gene_case_mixed):
    gene = gene_case_mixed
    altered = OpSequenceCase(TextCase.Unchanged)(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.loc["seq1"] == "gcaGTATAA"


def test_case_upper(gene_case_mixed):
    gene = gene_case_mixed
    altered = OpSequenceCase(TextCase.Upper)(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.loc["seq1"] == "GCAGTATAA"


def test_case_lower(gene_case_mixed):
    gene = gene_case_mixed
    altered = OpSequenceCase(TextCase.Lower)(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.loc["seq1"] == "gcagtataa"


def test_reverse_complement(gene_case_mixed):
    gene = gene_case_mixed
    altered = OpReverseComplement()(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.loc["seq1"] == "TTATACtgc"


@pytest.fixture
def gene_missing_gap() -> GeneSeries:
    series = pd.Series(
        {
            "seq1": "--GTnN?-TAA",
        },
        name="gene",
    )
    series.index.name = "seqid"
    return GeneSeries(series, missing="N?", gap="-")


def test_replace_missing(gene_missing_gap):
    gene = gene_missing_gap
    altered = OpTranslateMissing("?")(gene)
    assert altered.missing == "?"
    assert altered.series.loc["seq1"] == "--GTn??-TAA"


def test_replace_missing(gene_missing_gap):
    gene = gene_missing_gap
    altered = OpTranslateGap("*")(gene)
    assert altered.gap == "*"
    assert altered.series.loc["seq1"] == "**GTnN?*TAA"


def test_spreadsheet_compatibility(gene_missing_gap):
    gene = gene_missing_gap
    altered = OpSpreadsheetCompatibility()(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.loc["seq1"] == "N-GTnN?-TAA"


@pytest.fixture
def stream_simple() -> GeneStream:
    series1 = pd.Series(
        {
            "seq1": "ATCGCCTAA",
        },
        name="gene1",
    )
    series1.index.name = "seqid"
    gene1 = GeneSeries(series1)
    series2 = pd.Series(
        {
            "seq1": "GCCTAA",
        },
        name="gene2",
    )
    series2.index.name = "seqid"
    gene2 = GeneSeries(series2, reading_frame=ReadingFrame(-2))
    return GeneStream(iter([gene1, gene2]))


def test_filter(stream_simple):
    stream = stream_simple
    altered = stream.pipe(OpFilterGenes({"gene1"}))
    assert len(altered) == 1
    assert altered["gene1"]
    assert not altered["gene2"]


def test_stencil(stream_simple):
    stream = stream_simple
    altered = stream.pipe(OpStencilGenes(OpBlock(), {"gene1"}))
    assert len(altered) == 1
    assert not altered["gene1"]
    assert altered["gene2"]


def test_charsets(stream_simple):
    stream = stream_simple
    gene1 = stream["gene1"]
    gene2 = stream["gene2"]
    op_charsets = OpExtractCharsets()
    altered = stream.pipe(op_charsets)
    genes = list(altered)
    assert len(genes) == 2
    assert_gene_meta_equal(gene1, genes[0])
    assert_gene_meta_equal(gene2, genes[1])
    charsets = op_charsets.charsets
    assert charsets[0].name == "gene1"
    assert charsets[1].name == "gene2"
    assert charsets[0].frame == ReadingFrame(0)
    assert charsets[1].frame == ReadingFrame(-2)
    assert charsets[0].position == 1
    assert charsets[1].position == 10
    assert charsets[0].length == 9
    assert charsets[1].length == 6
    assert charsets[0].position_end == 9
    assert charsets[1].position_end == 15


@pytest.fixture
def stream_reading_frames() -> GeneStream:
    series1 = pd.Series(
        {
            "seq1": "ATCGCCTAA",
        },
        name="gene1",
    )
    series1.index.name = "seqid"
    gene1 = GeneSeries(series1, reading_frame=ReadingFrame(1))
    series2 = pd.Series(
        {
            "seq1": "GCCTAA",
        },
        name="gene2",
    )
    series2.index.name = "seqid"
    gene2 = GeneSeries(series2, reading_frame=ReadingFrame(-2), missing="n")
    series3 = pd.Series(
        {
            "seq1": "TAA",
        },
        name="gene3",
    )
    series3.index.name = "seqid"
    gene3 = GeneSeries(series3, reading_frame=ReadingFrame(3), missing="Nn?")
    return GeneStream(iter([gene1, gene2, gene3]))


def test_negative_reading_frames(stream_reading_frames):
    stream = stream_reading_frames
    altered = stream.pipe(OpReverseNegativeReadingFrames())
    assert len(altered) == 3
    assert altered["gene1"].reading_frame == ReadingFrame(1)
    assert altered["gene2"].reading_frame == ReadingFrame(2)
    assert altered["gene3"].reading_frame == ReadingFrame(3)
    assert altered["gene2"].series.loc["seq1"] == "TTAGGC"


def test_pad_reading_frames(stream_reading_frames):
    stream = stream_reading_frames
    altered = stream.pipe(OpPadReadingFrames())
    assert len(altered) == 3
    assert altered["gene1"].reading_frame == ReadingFrame(1)
    assert altered["gene2"].reading_frame == ReadingFrame(-1)
    assert altered["gene3"].reading_frame == ReadingFrame(1)
    assert altered["gene1"].series.loc["seq1"] == "ATCGCCTAA"
    assert altered["gene2"].series.loc["seq1"] == "GCCTAAnn"
    assert altered["gene3"].series.loc["seq1"] == "?TAA"


@pytest.fixture
def gene_meta_simple() -> GeneSeries:
    series = pd.Series(
        {
            "seq1": "TAA",
        },
        name="gene",
    )
    series.index.name = "seqid"
    return GeneSeries(series)


def test_replace_missing(gene_meta_simple):
    gene = gene_meta_simple
    metas = {
        "gene": {
            "missing": "X",
            "gap": "Y",
            "reading_frame": ReadingFrame(3),
            "genetic_code": GeneticCode(11),
            "names": ("g1", "g2", "g3"),
        }
    }
    altered = OpUpdateMetadata(metas)(gene)
    assert altered.missing == "X"
    assert altered.gap == "Y"
    assert altered.reading_frame == ReadingFrame(3)
    assert altered.genetic_code == GeneticCode(11)
    assert altered.names == ("g1", "g2", "g3")
    assert altered.series.loc["seq1"] == "TAA"


def test_make_uniform():
    series = pd.Series(
        {
            "seq1": "TAA",
            "seq2": "GCATAA",
        },
        name="gene",
    )
    series.index.name = "seqid"

    gene = GeneSeries(series)
    altered = OpMakeUniform("")(gene)
    assert_gene_meta_equal(gene, altered)
    assert altered.series.loc["seq1"] == "TAA"

    gene = GeneSeries(series, gap="*")
    altered = OpMakeUniform("-")(gene)
    assert_gene_meta_equal(gene, altered)
    assert altered.series.loc["seq1"] == "TAA---"

    gene = GeneSeries(series, reading_frame=ReadingFrame(-1))
    altered = OpMakeUniform("?")(gene)
    assert_gene_meta_equal(gene, altered)
    assert altered.series.loc["seq1"] == "???TAA"


def test_drop_empty():
    series = pd.Series(
        {
            "seq1": "---NNNTAA",
            "seq2": "N--?-***-nn-",
        },
        name="gene",
    )
    series.index.name = "seqid"

    gene = GeneSeries(series, missing="Nn?", gap="*-")
    altered = OpDropEmpty()(gene)
    assert_gene_meta_equal(gene, altered)
    assert len(altered.series) == 1
    assert altered.series.loc["seq1"] == "---NNNTAA"


def test_drop_empty_all():
    series = pd.Series(
        {
            "seq1": "----",
        },
        name="gene",
    )
    series.index.name = "seqid"

    gene = GeneSeries(series, missing="-")
    altered = OpDropEmpty()(gene)
    assert altered is None


def test_index_merge():
    series = pd.Series(
        {
            ("seq1", "voucher1", "locality1"): "ACGT",
        },
        name="gene",
    )
    series.index.names = ("seqid", "voucher", "locality")

    gene = GeneSeries(series)
    altered = OpIndexMerge(index="taxon", glue="+")(gene)
    assert isinstance(altered.series.index, pd.Index)
    assert altered.series.index.name == "taxon"
    assert len(altered.series.index) == 1
    assert altered.series.index[0] == "seq1+voucher1+locality1"


def test_general_info():
    TOTAL_DATA = pd.read_pickle(TEST_DATA_DIR / "test_general_info" / "total_data.pkl")
    BY_TAXON = pd.read_pickle(TEST_DATA_DIR / "test_general_info" / "by_taxon.pkl")
    genestream = read_from_path(TEST_DATA_DIR / "sequences.tab")
    operator = OpGeneralInfo()
    for _ in genestream.pipe(operator):
        pass
    table = operator.table
    print(table.dataframe.to_string())
    print(table.total_data().to_string())
    assert table.total_data().equals(TOTAL_DATA)
    print(table.by_taxon().to_string())
    assert table.by_taxon().equals(BY_TAXON)
    for taxa in table.disjoint_taxon_groups():
        print(taxa)
        assert taxa == {
            "sample4_Mantella crocea_FGZC346_Ranomafana",
            "sample3_Mantella aurantiaca_FGZC345_Andasibe",
            "sample1_Mantella aurantiaca_ZCMV1234_Andasibe",
            "sample2_Mantella aurantiaca_ZCMV9876_Andasibe",
            "sample5_Mantella crocea_MNHN1991_Ranomafana",
        }


def test_general_info_per_file():
    op = OpGeneralInfoPerFile()
    stream = read_from_path(TEST_DATA_DIR / "simple.tab")
    piped = stream.pipe(op)
    for _ in piped:
        pass
    stream = read_from_path(TEST_DATA_DIR / "sequences.tab")
    piped = stream.pipe(op)
    for _ in piped:
        pass

    table = op.get_info()
    assert table.iloc[0]["input file name"] == "simple.tab"
    assert table.iloc[0]["input file format"] == FileFormat.Tab
    assert table.iloc[0]["number of samples"] == 1
    assert table.iloc[0]["number of markers"] == 2
    assert table.iloc[0]["sequence length minimum"] == 4
    assert table.iloc[0]["sequence length maximum"] == 8
    assert table.iloc[0]["% missing nucleotides"] == 40.0
    assert table.iloc[0]["% GC content"] == 50.0

    assert table.iloc[1]["input file name"] == "sequences.tab"
    assert table.iloc[1]["input file format"] == FileFormat.Tab
    assert table.iloc[1]["number of samples"] == 5
    assert table.iloc[1]["number of markers"] == 3
    assert table.iloc[1]["sequence length minimum"] == 8
    assert table.iloc[1]["sequence length maximum"] == 10
    assert table.iloc[1]["% missing nucleotides"] == 0.0
    assert (table.iloc[1]["% GC content"] - 41 / 123 * 100) < 0.0001


def test_general_info_per_gene():
    stream = read_from_path(TEST_DATA_DIR / "per_gene.tab")
    op_general = OpGeneralInfo()
    op_per_gene = OpGeneralInfoPerGene()
    op_stencil_mafft = OpStencilGenes(OpTagSet("MafftRealigned", True), ["16S"])
    op_stencil_len = OpStencilGenes(OpTagSet("PaddedLength", True), ["16S"])
    op_stencil_codon = OpStencilGenes(OpTagSet("PaddedCodonPosition", True), ["cytb"])
    piped = (
        stream.pipe(op_stencil_mafft)
        .pipe(op_stencil_len)
        .pipe(op_stencil_codon)
        .pipe(op_general)
        .pipe(op_per_gene)
    )
    for _ in piped:
        pass
    table = op_per_gene.get_info(op_general.table)
    print(table.to_string())
    assert table.loc["16S"]["number of taxa with data"] == 3
    assert table.loc["16S"]["total number of nucleotides in alignment"] == 12
    assert table.loc["16S"]["% of missing data (nucleotides)"] == 50.0
    assert table.loc["16S"]["GC content of sequences"] == 50.0
    assert table.loc["16S"]["% of missing data (taxa)"] == 40.0
    assert table.loc["16S"]["re-aligned by Mafft yes/no"] == "yes"
    assert (
        table.loc["16S"]["padded to compensate for unequal sequence lengths yes/no"]
        == "yes"
    )
    assert table.loc["16S"]["padded to start with 1st codon position yes/no"] == "no"
    assert table.loc["cytb"]["number of taxa with data"] == 5
    assert table.loc["cytb"]["total number of nucleotides in alignment"] == 40
    assert (
        abs(table.loc["cytb"]["% of missing data (nucleotides)"] - 4 / 12 * 100)
        < 0.0001
    )
    assert table.loc["cytb"]["GC content of sequences"] == 50.0
    assert table.loc["cytb"]["% of missing data (taxa)"] == 0.0
    assert table.loc["cytb"]["re-aligned by Mafft yes/no"] == "no"
    assert (
        table.loc["cytb"]["padded to compensate for unequal sequence lengths yes/no"]
        == "no"
    )
    assert table.loc["cytb"]["padded to start with 1st codon position yes/no"] == "yes"


def test_general_info_disjoint():
    stream = read_from_path(TEST_DATA_DIR / "disjoint.tab")
    operator = OpGeneralInfo()
    for gene in stream.pipe(operator):
        print(gene.series)
    groups = list(operator.table.disjoint_taxon_groups())
    assert len(groups) == 2


def test_unconnected_taxons():
    stream = read_from_path(TEST_DATA_DIR / "long_with_missing2.tsv")
    operator = OpGeneralInfo()
    for _ in stream.pipe(operator):
        pass
    unconnected_pairs = list(operator.table.unconnected_taxons())
    assert len(unconnected_pairs) == 12


def test_drop_if_all_empty():
    series = pd.Series(
        {
            "seq1": "----",
            "seq2": "----",
            "seq3": "----",
        },
        name="gene",
    )
    series.index.name = "seqid"

    gene = GeneSeries(series, missing="-")
    altered = OpDropIfAllEmpty()(gene)
    assert altered is None


def test_drop_if_all_empty_kept():
    series = pd.Series(
        {
            "seq1": "ACGT",
            "seq2": "----",
            "seq3": "----",
        },
        name="gene",
    )
    series.index.name = "seqid"

    gene = GeneSeries(series, missing="-")
    altered = OpDropIfAllEmpty()(gene)
    assert len(altered.series) == 3


def test_long_missing_general_info():
    BY_TAXON = pd.read_pickle(TEST_DATA_DIR / "long_with_missing2_tsv.pkl")
    genestream = read_from_path(TEST_DATA_DIR / "long_with_missing2.tsv")
    operator = OpGeneralInfo()
    for _ in genestream.pipe(operator):
        pass
    table = operator.table
    assert table.by_taxon().equals(BY_TAXON)
