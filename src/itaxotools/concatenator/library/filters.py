
from typing import Callable, Dict, Iterator, List
from functools import reduce

import pandas as pd

from .utils import Stream, Filter, removeprefix
from .operators import index_to_multi


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
    for series in index_to_multi.to_filter(_join_species(stream)):
        yield series


def chain(filters: List[Filter]) -> Filter:
    return reduce(lambda f, g: lambda x: f(g(x)), filters)
