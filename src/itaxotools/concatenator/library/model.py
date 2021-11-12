
from __future__ import annotations
from typing import Iterator, Optional, Protocol
from itertools import tee
from copy import copy

import pandas as pd

from .utils import ConfigurableCallable, fill_empty
from .file_utils import PathLike
from .codons import GeneticCode, ReadingFrame


class BadGeneJoin(Exception):
    def __init__(self, existing: GeneSeries, adding: GeneSeries, error: str):
        self.existing = existing
        self.adding = adding
        super().__init__(error)


class GeneSeries:
    """Holds all sequences and metadata for a certain gene"""

    defaults = dict(
        codon_names = ('**_1st', '**_2nd', '**_3rd'),
        missing='Nn?',
        gap='-',
    )

    def __init__(
        self,
        series: Optional[pd.Series] = None,
        genetic_code: GeneticCode = GeneticCode(0),
        reading_frame: ReadingFrame = ReadingFrame(0),
        codon_names: Tuple[str, str, str] = defaults['codon_names'],
        missing: str = defaults['missing'],
        gap: str = defaults['gap'],
    ):
        self.series = series
        self.genetic_code = genetic_code
        self.reading_frame = reading_frame
        self.codon_names = codon_names
        self.missing = missing
        self.gap = gap

    def copy(self):
        other = copy(self)
        if self.series is not None:
            other.series = self.series.copy(deep=False)
        return other

    @property
    def name(self):
        if self.series is None:
            raise Exception('GeneSeries.series is None')
        return self.series.name

    @name.setter
    def name(self, value: str):
        if self.series is None:
            raise Exception('GeneSeries.series is None')
        self.series.name = value


class Operator(ConfigurableCallable):
    """The basis for manipulating GeneSeries"""

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        raise NotImplementedError()

    def iter(self, stream: GeneStream):
        for gene in stream:
            result = self(gene)
            if result is not None:
                yield result


class GeneStream:
    """An iterator over GeneSeries items with some extra functionalities"""

    def __init__(self, iterator: Iterator[GeneSeries]):
        self.iterator = iterator

    def __iter__(self):
        return self.iterator

    def __next__(self):
        return next(self.iterator)

    def __len__(self) -> int:
        # Used for testing, should be avoided otherwise
        self.iterator, iter = tee(self.iterator, 2)
        return sum(1 for _ in iter)

    def __getitem__(self, name: str) -> Optional[GeneSeries]:
        # Used for testing, should be avoided otherwise
        self.iterator, iter = tee(self.iterator, 2)
        matches = [gene for gene in iter if gene.name == name]
        if not matches:
            return None
        if len(matches) > 1:
            raise KeyError(f'Duplicate gene: {name}')
        return matches[0]

    def pipe(self, op: Operator) -> GeneStream:
        return GeneStream((op.iter(self)))


class GeneDataFrame:
    """Holds the sequences and metadata for multiple genes"""

    # this is generally a bad idea and should be fully replaced by GeneStream

    def __init__(self, dataframe: Optional[pd.DataFrame] = None, **kwargs):
        self._genes: Dict[str, GeneSeries] = dict()
        if dataframe is None:
            self.dataframe = pd.DataFrame(dtype=str)
        else:
            self.dataframe = dataframe
            for name in dataframe:
                series = pd.Series(dtype=str, name=name)
                self._genes[name] = GeneSeries(series, **kwargs)

    def __getitem__(self, name: str) -> GeneSeries:
        if name not in self._genes:
            raise Exception(f'Gene {repr(name)} not found')
        gene = self._genes[name].copy()
        gene.series = self.dataframe[name]
        return gene

    @classmethod
    def from_stream(cls, stream: GeneStream, filler = '') -> GeneDataFrame:
        obj = cls()
        for gene in stream:
            obj.add_gene(gene)
        if filler:
            assert len(filler) == 1
            for col in obj.dataframe:
                obj.dataframe[col] = fill_empty(obj.dataframe[col], filler)
        return obj

    def to_stream(self) -> GeneStream:
        return GeneStream(self[name] for name in self.dataframe)

    @property
    def genes(self):
        return list(self._genes.keys())

    def _join(self, series: pd.Series) -> pd.DataFrame:
        # This may drop rows when joining multiindex, should probably reset
        # unique indices and join the common ones, then set them back.
        # In case of no common index, should concatenate instead?
        df = self.dataframe.join(series, how='outer')
        df = df.fillna('')
        self.dataframe = df

    def _concat(self, gene: GeneSeries) -> None:
        existing = self[gene.name]
        if not all([
            existing.genetic_code == gene.genetic_code,
            existing.reading_frame == gene.reading_frame,
            existing.codon_names == gene.codon_names,
            existing.missing == gene.missing,
            existing.gap == gene.gap,
        ]):
            raise BadGeneJoin(
                existing, gene,
                f'Incompatible metadata for duplicate gene {repr(gene.name)}')
        both = pd.concat([existing.series, gene.series])
        if both.index.duplicated().any():
            raise BadGeneJoin(
                existing, gene,
                f'Duplicate index for duplicate gene {repr(gene.name)}')
        df = self.dataframe.drop(gene.name, axis=1)
        self.dataframe = df.join(both, how='outer')

    def add_gene(self, gene: GeneSeries) -> None:
        if gene.name in self.genes:
            self._concat(gene)
            return
        gene = gene.copy()
        if self.dataframe.empty:
            self.dataframe = pd.DataFrame(gene.series)
        else:
            self._join(gene.series)
        self._genes[gene.name] = gene
        gene.series = None


class GeneIO(Protocol):
    """
    Simply reads/writes the GeneSeries. Data processing is handled
    by the calling FileReader/FileWriter.
    """

    def gene_from_path(self, path: PathLike, **kwargs) -> GeneSeries:
        ...

    def gene_to_path(self, gene: GeneSeries, path: PathLike, **kwargs) -> None:
        ...


class StreamIO(Protocol):
    """
    Simply reads/writes the GeneStream. Data processing is handled
    by the calling FileReader/FileWriter.
    """

    def stream_from_path(self, path: PathLike, **kwargs) -> GeneStream:
        ...

    def stream_to_path(self, stream: GeneStream, path: PathLike,
                       **kwargs) -> None:
        ...
