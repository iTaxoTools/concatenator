
from typing import Callable, Dict, List, Iterator, Optional, Set

import pandas as pd

from itaxotools.DNAconvert.library.utils import sanitize

from .model import Operator, GeneSeries, GeneStream, GeneDataFrame
from .types import TextCase, Charset
from .utils import (
    Translation, Field, OrderedSet, removeprefix,
    reverse_complement, has_uniform_length)
from .codons import final_column_reading_frame, ReadingFrame


class InvalidGeneSeries(Exception):
    def __init__(self, gene: GeneSeries, error: str):
        self.gene = gene
        super().__init__(error)


class OpPass(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        return gene


class OpBlock(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        return None


class OpCheckDuplicates(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._memory = set()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name in self._memory:
            raise InvalidGeneSeries(gene, f'Duplicate gene: {repr(gene.name)}')
        self._memory.add(gene.name)
        return gene


class OpCheckValid(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.op_check_duplicates = OpCheckDuplicates()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if not isinstance(gene, GeneSeries):
            raise InvalidGeneSeries(gene, 'Not a GeneSeries!')
        if not isinstance(gene.series, pd.Series):
            raise InvalidGeneSeries(gene, 'Gene series is not pandas.Series!')
        if gene.series.index.duplicated().any():
            raise InvalidGeneSeries(gene, 'Duplicate indices')
        if not gene.name:
            raise InvalidGeneSeries(gene, 'Missing gene data: "name"')
        if not gene.missing:
            raise InvalidGeneSeries(gene, 'Missing gene data: "missing"')
        if not gene.gap:
            raise InvalidGeneSeries(gene, 'Missing gene data: "gap"')
        return self.op_check_duplicates(gene)


class OpIndexMerge(Operator):
    index: str = 'seqid'
    glue: str = Field('glue', value='_')

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        indices = gene.series.index.to_frame().fillna('', inplace=False)
        gene.series.index = pd.Index(
            indices.apply(self.glue.join, axis=1),
            name=self.index)
        return gene


class OpIndexFilter(Operator):
    index: str = 'seqid'

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series.index = gene.series.index.to_frame()[self.index]
        return gene


class OpDropEmpty(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series = gene.series.dropna(inplace=False)
        if gene.missing:
            gene.series = gene.series[
                ~ gene.series.str.fullmatch(f'[{gene.missing}]+')]
        gene.series = gene.series[gene.series.str.len() > 0]
        return gene


class OpPadRight(Operator):
    padding: str = Field('padding', value='-')

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if not self.padding:
            return gene
        assert len(self.padding) == 1
        gene = gene.copy()
        gene.series = gene.series.fillna('', inplace=False)
        max_length = gene.series.str.len().max()
        gene.series = gene.series.str.ljust(max_length, self.padding)
        return gene


class OpTranslateSequences(Operator):
    translation: Translation = Field('translation', value={})

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series = gene.series.str.translate(self.translation)
        return gene


class OpTranslateMissing(Operator):
    missing: str = Field('missing', value='?')

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if not self.missing or not gene.missing:
            return gene
        assert len(self.missing) == 1
        mapping = {char: self.missing for char in gene.missing}
        translation = str.maketrans(mapping)
        gene = OpTranslateSequences(translation)(gene)
        gene.missing = self.missing
        return gene


class OpTranslateGap(Operator):
    gap: str = Field('gap', value='-')

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if not self.gap or not gene.gap:
            return gene
        assert len(self.gap) == 1
        mapping = {char: self.gap for char in gene.gap}
        translation = str.maketrans(mapping)
        gene = OpTranslateSequences(translation)(gene)
        gene.gap = self.gap
        return gene


class OpFilterGenes(Operator):
    genes: Set = Field('genes', value=set())

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name not in self.genes:
            return None
        return gene


class OpStencilGenes(Operator):
    operator: Operator = Field('operator', value=OpPass())
    genes: Set = Field('genes', value=set())

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name in self.genes:
            return self.operator(gene)
        return gene


class OpTranslateGenes(Operator):
    translation: Dict = Field('translation', value=dict())

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name in self.translation:
            if self.translation[gene.name] is None:
                return None
            gene = gene.copy()
            gene.name = self.translation[gene.name]
        return gene


class OpChainGenes(Operator):
    allow_duplicates: bool = Field('allow_duplicates', value=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._memory = set()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if self.allow_duplicates:
            return gene
        if gene.name in self._memory:
            return None
        self._memory.add(gene.name)
        return gene


class OpDetectReadingFrame(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        final_reading_frame = final_column_reading_frame(
            gene.series, gene.genetic_code, gene.reading_frame)
        gene.reading_frame = final_reading_frame
        return gene


class OpSanitizeGeneNames(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.name = sanitize(gene.name)
        return gene


class OpSanitizeSpeciesNames(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        indices = gene.series.index.names
        data = gene.series.reset_index()
        data[indices] = data[indices].applymap(sanitize)
        gene.series = data.set_index(indices).iloc[:, 0]
        return gene


class OpSequenceCase(Operator):
    case: TextCase = Field('case', value=TextCase.Unchanged)

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if self.case.method is None:
            return gene
        gene = gene.copy()
        gene.series = gene.series.apply(self.case.method)
        return gene


class OpSpreadsheetCompatibility(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series = gene.series.str.replace('^-', 'N', regex=True)
        return gene


class OpReverseComplement(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series = gene.series.map(reverse_complement)
        return gene


class OpReverseNegativeReadingFrames(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.reading_frame >= 0:
            return gene
        gene = gene.copy()
        gene.reading_frame = ReadingFrame(-gene.reading_frame)
        return OpReverseComplement()(gene)


class OpPadReadingFrames(Operator):
    padding: str = Field('padding', value='')
    lengths = {2: 2, 3: 1}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert len(self.padding) in [0, 1]

    @staticmethod
    def pad_left(seq: str, pad: str) -> str:
        return f'{pad}{seq}'

    @staticmethod
    def pad_right(seq: str, pad: str) -> str:
        return f'{seq}{pad}'

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.reading_frame in [0, 1, -1]:
            return gene
        gene = gene.copy()
        func = self.pad_left if gene.reading_frame > 0 else self.pad_right
        reading_frame = abs(gene.reading_frame)
        if self.padding:
            padding = self.padding
        else:
            padding = '?' if '?' in gene.missing else gene.missing[0]
        padding = padding * self.lengths[reading_frame]
        gene.series = gene.series.apply(lambda seq: func(seq, padding))
        sign = 1 - 2 * int(gene.reading_frame < 0)
        gene.reading_frame = ReadingFrame(sign)
        return gene


class OpExtractCharsets(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.charsets = list()
        self.cursor = 1

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        assert has_uniform_length(gene.series)
        length = len(gene.series.iat[0])
        charset = Charset(
            gene.name, self.cursor, length,
            gene.reading_frame, gene.codon_names)
        self.charsets.append(charset)
        self.cursor += length
        return gene


class OpUpdateMetadata(Operator):
    metas: Dict[str, Dict] = Field('metas', value={})

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name not in self.metas:
            return gene
        gene = gene.copy()
        for k, v in self.metas[gene.name].items():
            setattr(gene, k, v)
        return gene


class OpApplyToGene(Operator):
    func: Callable[[GeneSeries], GeneSeries] = Field('func', value=None)

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if self.func is None:
            return gene
        return self.func(gene)


class OpApplyToSeries(Operator):
    func: Callable[[pd.Series], pd.Series] = Field('func', value=None)

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if self.func is None:
            return gene
        gene = gene.copy()
        gene.series = self.func(gene.series)
        return gene


# Pending removal, functionality to be merged into GeneDataFrame?
def _join_any(stream: GeneStream) -> GeneDataFrame:
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
