
from typing import Callable, Dict, Iterator, List
from functools import reduce

import pandas as pd

from .utils import removeprefix


Stream = Iterator[pd.Series]
Filter = Callable[[Stream], Stream]

# Such as returned by str.maketrans
Translation = Dict[int, int]


class InvalidSeries(Exception):
    def __init__(self, series: pd.Series, error: str):
        self.series = series
        print(series)
        super().__init__(error)


def check_valid(stream: Stream) -> Stream:
    for series in stream:
        if not isinstance(series, pd.Series):
            raise InvalidSeries(series, 'Not a pandas.Series!')
        if not series.name:
            raise InvalidSeries(series, 'Missing series name')
        if not isinstance(series.index, pd.MultiIndex):
            raise InvalidSeries(series, 'Series must have MultiIndex')
        if series.index.names[0] != 'species':
            raise InvalidSeries(series, 'Missing "species" index')
        yield series


def index_to_multi(stream: Stream) -> Stream:
    for series in stream:
        if not isinstance(series.index, pd.MultiIndex):
            series.index = pd.MultiIndex.from_arrays(
                [series.index], names=[series.index.name])
        yield series


def _join_species(stream: Stream) -> Stream:
    """Join by species"""
    sentinel = '\u0000'
    all_extras = set()

    def guard(names: Iterator[str]) -> List[str]:
        return [sentinel + name for name in names]

    def unguard(names: Iterator[str]) -> List[str]:
        return [removeprefix(name, sentinel) for name in names]

    def by_species(stream: Stream) -> Stream:
        for series in stream:
            series.index.names = guard(series.index.names)
            extras = series.index.names[1:]  # only keep species for index
            all_extras.update(extras)
            yield series.reset_index(extras)

    all = pd.concat(by_species(stream), axis=1)
    all = all.set_index(list(all_extras), append=True)
    all.index.names = unguard(all.index.names)
    for col in all:
        yield all[col].fillna('')


def join(stream: Stream) -> Stream:
    for series in index_to_multi(_join_species(stream)):
        yield series


def translate(translation: Translation) -> Filter:
    def _translate(stream: Stream) -> Stream:
        for series in stream:
            yield series.str.translate(translation)
    return _translate


def chain(filters: List[Filter]) -> Filter:
    return reduce(lambda f, g: lambda x: f(g(x)), filters)
