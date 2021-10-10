
from typing import Callable, Union

import pandas as pd

from .utils import Filter, Stream, Translation

from . import SPECIES


IndexedData = Union[pd.DataFrame, pd.Series]


class Operator:
    def __init__(self, func: Callable[[pd.Series], pd.Series]):
        self.func = func

    def __call__(self, series: pd.Series) -> pd.Series:
        return self.func(series)

    @property
    def to_filter(self) -> Filter:
        def _filter(stream: Stream) -> Stream:
            for series in stream:
                yield self(series)
        return _filter


def operator(func: Callable[[pd.Series], pd.Series]) -> Operator:
    return Operator(func)


class InvalidSeries(Exception):
    def __init__(self, series: pd.Series, error: str):
        self.series = series
        super().__init__(error)


@operator
def check_valid(series: pd.Series) -> pd.Series:
    if not isinstance(series, pd.Series):
        raise InvalidSeries(series, 'Not a pandas.Series!')
    if not series.name:
        raise InvalidSeries(series, 'Missing series name')
    if not isinstance(series.index, pd.MultiIndex):
        raise InvalidSeries(series, 'Series must have MultiIndex')
    if SPECIES not in series.index.names:
        raise InvalidSeries(series, f'Missing "{SPECIES}" index')
    if series.index.names[0] != SPECIES:
        raise InvalidSeries(series, f'"{SPECIES}" index must be first')
    return series


@operator
def index_to_multi(series: IndexedData) -> IndexedData:
    if not isinstance(series.index, pd.MultiIndex):
        series.index = pd.MultiIndex.from_arrays(
            [series.index], names=[series.index.name])
    return series


@operator
def species_to_front(series: IndexedData) -> IndexedData:
    ordered = list(series.index.names)
    ordered.remove(SPECIES)
    ordered = [SPECIES] + ordered
    return series.reorder_levels(ordered)


def translate(translation: Translation) -> Filter:
    @operator
    def _translate(series: pd.Series) -> pd.Series:
        return series.str.translate(translation)
    return _translate
