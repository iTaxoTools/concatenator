from pathlib import Path
from typing import Iterator
import random

import pytest
import pandas as pd
import numpy as np

from itaxotools.concatenator.library.model import GeneSeries, GeneStream
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


def assert_gene_meta_equal(gene1, gene2):
    assert gene1.missing == gene2.missing
    assert gene1.gap == gene2.gap
    assert gene1.genetic_code == gene2.genetic_code
    assert gene1.reading_frame == gene2.reading_frame
    assert gene1.codon_names == gene2.codon_names


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
    genestream = read_from_path(Path(__file__).with_name("sequences.tab"))
    operator = OpGeneralInfo()
    for _ in genestream.pipe(operator):
        pass
    table = operator.table
    print(table.dataframe.to_string())
    print(table.total_data().to_string())
    print(table.by_taxon().to_string())
    for taxa in table.disjoint_taxon_groups():
        print(taxa)
    genestream = read_from_path(Path(__file__).with_name("sequences.tab"))
    filenames = ["foo", "bar", "baz"]
    file_formats = [FileFormat.Tab, FileFormat.Nexus, FileFormat.Ali]

    def stream_with_files() -> Iterator[FileGeneralInfo]:
        for filename, file_format, gene in zip(filenames, file_formats, genestream):
            op = OpGeneralInfo()
            op(gene)
            yield FileGeneralInfo(filename, file_format, op.table)

    print(GeneralInfo.by_input_file(stream_with_files()).to_string())
    genestream = read_from_path(Path(__file__).with_name("sequences.tab"))
    gene_index = pd.Index([gene.name for gene in genestream])
    gene_index.name = InfoColumns.Gene
    gene_info = GeneInfo(
        pd.DataFrame(
            {
                GeneInfoColumns.MafftRealigned: random.choices(
                    [True, False], k=len(gene_index)
                ),
                GeneInfoColumns.PaddedLength: random.choices(
                    [True, False], k=len(gene_index)
                ),
                GeneInfoColumns.PaddedCodonPosition: random.choices(
                    [True, False], k=len(gene_index)
                ),
            },
            index=gene_index,
        )
    )
    print(table.by_gene(gene_info).to_string())
