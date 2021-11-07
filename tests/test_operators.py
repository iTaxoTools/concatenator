
import pytest
import pandas as pd

from itaxotools.concatenator.library.model import GeneSeries, GeneStream
from itaxotools.concatenator.library.types import TextCase
from itaxotools.concatenator.library.codons import GeneticCode, ReadingFrame
from itaxotools.concatenator.library.operators import (
    OpSanitizeGeneNames, OpSanitizeSpeciesNames, OpSequenceCase,
    OpTranslateMissing, OpTranslateGap, OpSpreadsheetCompatibility,
    OpFilterGenes, OpStencilGenes, OpBlock, OpReverseComplement)


def assert_gene_meta_equal(gene1, gene2):
    assert gene1.missing == gene2.missing
    assert gene1.gap == gene2.gap
    assert gene1.genetic_code == gene2.genetic_code
    assert gene1.reading_frame == gene2.reading_frame
    assert gene1.codon_names == gene2.codon_names


@pytest.fixture
def gene_unsanitized() -> GeneSeries:
    series = pd.Series({
            'ΔšëqΔ¹Δ': 'GCAGTATAA',
        }, name='ΔgëneΔ²Δ')
    series.index.name = 'seqid'
    return GeneSeries(series)


def test_sanitize_genes(gene_unsanitized):
    gene = gene_unsanitized
    altered = OpSanitizeGeneNames()(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.name == 'gene_2'


def test_sanitize_species(gene_unsanitized):
    gene = gene_unsanitized
    altered = OpSanitizeSpeciesNames()(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.index[0] == 'seq_1'


@pytest.fixture
def gene_case_mixed() -> GeneSeries:
    series = pd.Series({
            'seq1': 'gcaGTATAA',
        }, name='gene')
    series.index.name = 'seqid'
    return GeneSeries(series)


def test_case_unchanged(gene_case_mixed):
    gene = gene_case_mixed
    altered = OpSequenceCase(TextCase.Unchanged)(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.iloc[0] == 'gcaGTATAA'


def test_case_upper(gene_case_mixed):
    gene = gene_case_mixed
    altered = OpSequenceCase(TextCase.Upper)(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.iloc[0] == 'GCAGTATAA'


def test_case_lower(gene_case_mixed):
    gene = gene_case_mixed
    altered = OpSequenceCase(TextCase.Lower)(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.iloc[0] == 'gcagtataa'


def test_reverse_complement(gene_case_mixed):
    gene = gene_case_mixed
    altered = OpReverseComplement()(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.iloc[0] == 'TTATACtgc'


@pytest.fixture
def gene_missing_gap() -> GeneSeries:
    series = pd.Series({
            'seq1': '--GTAN?-TAA',
        }, name='gene')
    series.index.name = 'seqid'
    return GeneSeries(series, missing='N?', gap='-')


def test_replace_missing(gene_missing_gap):
    gene = gene_missing_gap
    altered = OpTranslateMissing('?')(gene)
    assert altered.missing == '?'
    assert altered.series.iloc[0] == '--GTA??-TAA'


def test_replace_missing(gene_missing_gap):
    gene = gene_missing_gap
    altered = OpTranslateGap('*')(gene)
    assert altered.gap == '*'
    assert altered.series.iloc[0] == '**GTAN?*TAA'


def test_spreadsheet_compatibility(gene_missing_gap):
    gene = gene_missing_gap
    altered = OpSpreadsheetCompatibility()(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.series.iloc[0] == 'N-GTAN?-TAA'


@pytest.fixture
def stream_simple() -> GeneStream:
    series1 = pd.Series({
            'seq1': 'ATCGCCTAA',
        }, name='gene1')
    series1.index.name = 'seqid'
    gene1 = GeneSeries(series1)
    series2 = pd.Series({
            'seq1': 'GCCTAA',
        }, name='gene2')
    series2.index.name = 'seqid'
    gene2 = GeneSeries(series2)
    return GeneStream(iter([gene1, gene2]))


def test_filter(stream_simple):
    stream = stream_simple
    altered = stream.pipe(OpFilterGenes({'gene1'}))
    assert len(altered) == 1
    assert altered['gene1']
    assert not altered['gene2']


def test_stencil(stream_simple):
    stream = stream_simple
    altered = stream.pipe(OpStencilGenes(OpBlock(), {'gene1'}))
    assert len(altered) == 1
    assert not altered['gene1']
    assert altered['gene2']
