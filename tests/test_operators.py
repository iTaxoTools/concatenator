
import pytest
import pandas as pd

from itaxotools.concatenator.library.model import GeneSeries, GeneStream
from itaxotools.concatenator.library.codons import GeneticCode, ReadingFrame
from itaxotools.concatenator.library.operators import (
    OpSanitizeGeneNames, )


def assert_gene_meta_equal(gene1, gene2):
    assert gene1.missing == gene2.missing
    assert gene1.gap == gene2.gap
    assert gene1.genetic_code == gene2.genetic_code
    assert gene1.reading_frame == gene2.reading_frame
    assert gene1.codon_names == gene2.codon_names


def test_sanitize_genes():
    series = pd.Series({
            'seq1': 'GCAGTATAA',
        }, name='ΔgëneΔ²Δ')
    series.index.name = 'seqid'
    gene = GeneSeries(series)

    altered = OpSanitizeGeneNames()(gene)
    assert_gene_meta_equal(altered, gene)
    assert altered.name == 'gene_2'
