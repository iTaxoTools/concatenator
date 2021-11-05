
import pytest
import pandas as pd

from itaxotools.concatenator.library.model import (
    GeneSeries, GeneDataFrame, BadGeneJoin)
from itaxotools.concatenator.library.codons import GeneticCode, ReadingFrame


def assert_gene_meta_equal(gene1, gene2):
    assert gene1.missing == gene2.missing
    assert gene1.gap == gene2.gap
    assert gene1.genetic_code == gene2.genetic_code
    assert gene1.reading_frame == gene2.reading_frame
    assert gene1.codon_names == gene2.codon_names


def test_gene_dataframe_simple():
    series1 = pd.Series({
            'seq1': 'GCAGTATAA',
        }, name='gene1')
    series1.index.name = 'seqid'
    gene1 = GeneSeries(series1, missing='?', reading_frame = ReadingFrame(0))

    series2 = pd.Series({
            'seq1': 'GGCTAA',
            'seq2': 'GCCTAA',
        }, name='gene2')
    series2.index.name = 'seqid'
    gene2 = GeneSeries(series2, missing='N', reading_frame = ReadingFrame(1))

    stream = [gene1, gene2]
    gdf = GeneDataFrame.from_stream(stream)
    assert gdf['gene1'].series.loc['seq1'] == series1.loc['seq1']
    assert pd.isnull(gdf['gene1'].series.loc['seq2'])
    assert gdf['gene2'].series.loc['seq1'] == series2.loc['seq1']
    assert gdf['gene2'].series.loc['seq2'] == series2.loc['seq2']
    assert_gene_meta_equal(gdf['gene1'], gene1)
    assert_gene_meta_equal(gdf['gene2'], gene2)


def test_gene_dataframe_existing():
    series1 = pd.Series({
            'seq1': 'GCAGTATAA',
        }, name='gene')
    series1.index.name = 'seqid'
    gene1 = GeneSeries(series1)

    series2 = pd.Series({
            'seq2': 'GCCTAA',
        }, name='gene')
    series2.index.name = 'seqid'
    gene2 = GeneSeries(series2)

    stream = [gene1, gene2]
    gdf = GeneDataFrame.from_stream(stream)
    assert gdf['gene'].series.loc['seq1'] == series1.loc['seq1']
    assert gdf['gene'].series.loc['seq2'] == series2.loc['seq2']
    assert_gene_meta_equal(gdf['gene'], gene1)
    assert_gene_meta_equal(gdf['gene'], gene2)


def test_gene_dataframe_incompatible():
    series1 = pd.Series({
            'seq1': 'GCAGTATAA',
        }, name='gene1')
    series1.index.name = 'seqid'
    gene1 = GeneSeries(series1, missing='?')

    series2 = pd.Series({
            'seq2': 'GCCTAA',
        }, name='gene1')
    series2.index.name = 'seqid'
    gene2 = GeneSeries(series2, missing='N')

    with pytest.raises(BadGeneJoin):
        stream = [gene1, gene2]
        gdf = GeneDataFrame.from_stream(stream)


def test_gene_dataframe_duplicate():
    series1 = pd.Series({
            'seq1': 'GCAGTATAA',
        }, name='gene1')
    series1.index.name = 'seqid'
    gene1 = GeneSeries(series1, missing='?')

    series2 = pd.Series({
            'seq1': 'GCCTAA',
        }, name='gene1')
    series2.index.name = 'seqid'
    gene2 = GeneSeries(series2, missing='?')

    with pytest.raises(BadGeneJoin):
        stream = [gene1, gene2]
        gdf = GeneDataFrame.from_stream(stream)


def test_reading_frame():
    assert not bool(ReadingFrame(0))
    assert ReadingFrame(1) == 1
    assert ReadingFrame(1).label == '+1'
    assert ReadingFrame(-2).label == '-2'


def test_genetic_code():
    assert not bool(GeneticCode(0))
    assert not GeneticCode(0).stops
    assert GeneticCode(1) == 1
    assert GeneticCode(1).stops
    assert GeneticCode(1).text == 'Standard'
    assert GeneticCode(0).text == 'Unknown'
