
from __future__ import annotations
from typing import Callable, Iterator, Optional, Protocol, TextIO
from copy import copy

import pandas as pd

from .utils import ConfigurableCallable
from .file_utils import PathLike


class GeneDataFrame:
    """Holds the sequences and metadata for multiple genes"""

    # this is generally a bad idea and should be fully replaced by GeneStream
    defaults = dict(
        missing = '?N',
        gap = '-',
    )

    def __init__(
        self,
        dataframe: pd.DataFrame,
        missing: str = defaults['missing'],
        gap: str = defaults['gap'],
    ):
        self.dataframe = dataframe
        self.missing = missing
        self.gap = gap

    def stream(self):
        return GeneStream.from_dataframe(
            self.dataframe, missing=self.missing, gap=self.gap)

    @classmethod
    def from_gene(cls, gene: GeneSeries):
        return cls(
            pd.DataFrame(gene.series),
            missing=gene.missing,
            gap=gene.gap)


class GeneSeries:
    """Holds all sequences and metadata for a certain gene"""

    defaults = dict(
        # codons = ('**_1st', '**_2nd', '**_3rd'),
        missing = '?N',
        gap = '-',
    )

    def __init__(
        self,
        series: pd.Series,
        # genetic_code: GeneticCode = GeneticCode(0),
        # reading_frame: ReadingFrame = ReadingFrame(0),
        # codons: Tuple[str, str, str] = self.defaults['codons'],
        missing: str = defaults['missing'],
        gap: str = defaults['gap'],
    ):
        self.series = series
        # self.genetic_code = genetic_code
        # self.reading_frame = reading_frame
        # self.codons = codons
        self.missing = missing
        self.gap = gap

    def copy(self):
        other = copy(self)
        other.series = self.series.copy(deep=False)
        return other

    @property
    def name(self):
        return self.series.name

    @name.setter
    def name(self, value: str):
        self.series.name = value


class Operator(ConfigurableCallable):
    """The basis for manipulating GeneSeries"""

    def call(self, series: GeneSeries) -> Optional[GeneSeries]:
        raise NotImplementedError()

    def iter(self, stream: GeneStream):
        for series in stream:
            result = self(series)
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

    def pipe(self, op: Operator) -> GeneStream:
        return GeneStream((op.iter(self)))

    @classmethod
    def from_dataframe(cls, df: pd.Dataframe, **kwargs):
        return cls((GeneSeries(df[col], **kwargs) for col in df))


class GeneIO(Protocol):
    """
    Simply reads/writes the GeneSeries. Data processing is handled
    by the calling FileReader/FileWriter.
    """

    def gene_from_path(self, path: PathLike) -> GeneSeries:
        ...

    def gene_to_path(self, gene: GeneSeries, path: PathLike) -> None:
        ...


class StreamIO(Protocol):
    """
    Simply reads/writes the GeneStream. Data processing is handled
    by the calling FileReader/FileWriter.
    """

    def stream_from_path(self, path: PathLike, **kwargs) -> GeneStream:
        ...

    def stream_to_path(self, gene: GeneStream, path: PathLike) -> None:
        ...
