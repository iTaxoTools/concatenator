
from typing import Iterator, List
from functools import reduce

import pandas as pd

from .utils import Stream, Filter, removeprefix
from .operators import index_to_multi, species_to_front


def chain(filters: List[Filter]) -> Filter:
    return reduce(lambda f, g: lambda x: f(g(x)), filters)


def join_any(stream: Stream) -> Stream:
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
    for col in all:
        yield all[col]


def join(stream: Stream) -> Stream:
    filters = chain([species_to_front.to_filter, index_to_multi.to_filter])
    for series in filters(join_any(stream)):
        yield series
