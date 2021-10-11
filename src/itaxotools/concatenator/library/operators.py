
from typing import Callable, List, Iterator, Union

import pandas as pd

from .utils import Filter, Stream, Translation, removeprefix, make_equal_length

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


def index_merge(glue: str = '_') -> Operator:
    @operator
    def _index_merge(series: IndexedData) -> IndexedData:
        series.index = series.index.to_frame().apply(glue.join, axis=1)
        return series
    return _index_merge


@operator
def index_species_only(series: IndexedData) -> IndexedData:
    series.index = series.index.to_frame()[SPECIES]
    return series


@operator
def drop_empty(series: IndexedData) -> IndexedData:
    series.dropna(inplace=True)
    return series[series.str.len() > 0]


def pad_right(fillchar: str = '-') -> Operator:
    @operator
    def _pad_right(series: pd.Series) -> pd.Series:
        return make_equal_length(series.dropna(), fillchar=fillchar)
    return _pad_right


def translate(translation: Translation) -> Operator:
    @operator
    def _translate(series: pd.Series) -> pd.Series:
        return series.str.translate(translation)
    return _translate


def join_any(stream: Stream) -> pd.DataFrame:
    """Outer join for any MultiIndex"""
    sentinel = '\u0000'
    all_keys = set()

    def guarded(names: Iterator[str]) -> List[str]:
        return [name for name in names if name.startswith(sentinel)]

    def guard(names: Iterator[str]) -> List[str]:
        return [sentinel + name for name in names]

    def unguard(names: Iterator[str]) -> List[str]:
        return [removeprefix(name, sentinel) for name in names]

    def fold_keys(stream: Stream) -> Iterator[pd.DataFrame]:
        for series in stream:
            series.index.names = guard(series.index.names)
            keys = series.index.names
            all_keys.update(keys)
            yield series.reset_index(keys)

    species = fold_keys(stream)
    all = pd.DataFrame(next(species))
    for series in species:
        merge_keys = set(guarded(series.columns)) & all_keys
        all = pd.merge(all, series, how='outer', on=list(merge_keys))
    all.set_index(list(all_keys), inplace=True)
    all.index.names = unguard(all.index.names)
    return species_to_front(index_to_multi(all))
