
from typing import Callable, Dict, List, Iterator, TypeVar, Set, Union
from functools import reduce

import pandas as pd

from .model import Operator, GeneSeries, GeneStream, GeneDataFrame
from .utils import (
    Translation, ConfigurableCallable, Param,
    OrderedSet, removeprefix,
    )


T = TypeVar('T')


# to be removed, replaced by GeneStream.pipe
def chain(funcs: List[Callable[[T], T]]) -> Callable[[T], T]:
    return reduce(lambda f, g: lambda x: f(g(x)), funcs)


class InvalidGeneSeries(Exception):
    def __init__(self, gene: GeneSeries, error: str):
        self.gene = gene
        super().__init__(error)


class OpCheckValid(Operator):
    def call(self, gene: GeneSeries) -> GeneSeries:
        if not isinstance(gene, GeneSeries):
            raise InvalidGeneSeries(gene, 'Not a GeneSeries!')
        if not isinstance(gene.series, pd.Series):
            raise InvalidGeneSeries(gene, 'Gene series is not pandas.Series!')
        if not gene.name:
            raise InvalidGeneSeries(gene, 'Missing gene data: "name"')
        if not gene.missing:
            raise InvalidGeneSeries(gene, 'Missing gene data: "missing"')
        if not gene.gap:
            raise InvalidGeneSeries(gene, 'Missing gene data: "gap"')
        return gene


class OpIndexMerge(Operator):
    index: str = 'seqid'
    glue: str = Param('_')

    def call(self, gene: GeneSeries) -> GeneSeries:
        gene = gene.copy()
        indices = gene.series.index.to_frame().fillna('', inplace=False)
        gene.series.index = pd.Index(
            indices.apply(self.glue.join, axis=1),
            name=self.index)
        return gene


class OpIndexFilter(Operator):
    index: str = 'seqid'

    def call(self, gene: GeneSeries) -> GeneSeries:
        gene = gene.copy()
        gene.series.index = gene.series.index.to_frame()[self.index]
        return gene


class OpDropEmpty(Operator):
    missing: str = Param('')

    def call(self, gene: GeneSeries) -> GeneSeries:
        gene = gene.copy()
        gene.series = gene.series.dropna(inplace=False)
        if self.missing:
            gene.series = gene.series[
                ~ gene.series.str.fullmatch(f'[{self.missing}]+')]
        gene.series = gene.series[gene.series.str.len() > 0]
        return gene


class OpPadRight(Operator):
    padding: str = Param('-')

    def call(self, gene: GeneSeries) -> GeneSeries:
        if not self.padding:
            return gene
        assert len(self.padding) == 1
        gene = gene.copy()
        gene.series = gene.series.fillna('', inplace=False)
        max_length = gene.series.str.len().max()
        gene.series = gene.series.str.ljust(max_length, self.padding)
        return gene


class OpTranslateSequences(Operator):
    translation: Translation = Param({})

    def call(self, gene: GeneSeries) -> GeneSeries:
        gene = gene.copy()
        gene.series = series.str.translate(self.translation)
        return gene


class OpFilterCharsets(Operator):
    filter: Set = Param(set())

    def call(self, gene: GeneSeries) -> GeneSeries:
        if gene.name not in self.filter:
            return None
        return gene


class OpTranslateGenes(Operator):
    translation: Dict = Param(dict())

    def call(self, gene: GeneSeries) -> GeneSeries:
        if gene.name in self.translation:
            if self.translation[gene.name] is None:
                return None
            gene = gene.copy()
            gene.name = self.translation[gene.name]
        return gene


class OpChainCharsets(Operator):
    allow_duplicates: bool = Param(False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._memory = set()

    def call(self, gene: GeneSeries) -> GeneSeries:
        if self.allow_duplicates:
            return gene
        if gene.name in gene._memory:
            return None
        self._memory.add(gene.name)
        return gene


class OpApplyToGene(Operator):
    func: Callable[[GeneSeries], GeneSeries] = Param(None)

    def call(self, gene: GeneSeries) -> GeneSeries:
        if self.func is None:
            return gene
        return self.func(gene)


class OpApplyToSeries(Operator):
    func: Callable[[pd.Series], pd.Series] = Param(None)

    def call(self, gene: GeneSeries) -> GeneSeries:
        if self.func is None:
            return gene
        gene = gene.copy()
        gene.series = self.func(gene.series)
        return gene


def join_any(stream: GeneStream) -> GeneDataFrame:
    # To be replaced by a better method, using a file buffer
    sentinel = '\u0000'
    all_keys = OrderedSet()

    def guarded(names: Iterator[str]) -> List[str]:
        return [name for name in names if name.startswith(sentinel)]

    def guard(names: Iterator[str]) -> List[str]:
        return [sentinel + name for name in names]

    def unguard(names: Iterator[str]) -> List[str]:
        return [removeprefix(name, sentinel) for name in names]

    def fold_keys(stream: GeneStream) -> Iterator[GeneSeries]:
        for gene in stream:
            gene = gene.copy()
            keys = guard(gene.series.index.names)
            gene.series.index = pd.MultiIndex.from_frame(
                gene.series.index.to_frame(), names=keys)
            all_keys.update(keys)
            gene.series = gene.series.reset_index(keys)
            yield gene

    genes = fold_keys(stream)
    gdf = GeneDataFrame.from_gene(next(genes))
    all = gdf.dataframe
    for gene in genes:
        merge_keys = all_keys & guarded(gene.series.columns)
        all = pd.merge(all, gene.series, how='outer', on=list(merge_keys))
        # We assume all gene series are compatible
        print(gene, gene.missing)
        gdf.missing = gene.missing
        gdf.gap = gene.gap
    print('end', gdf.missing)
    all.set_index(list(all_keys), inplace=True)
    all.index.names = unguard(all.index.names)
    gdf.dataframe = all
    return gdf
