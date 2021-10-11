
from typing import Callable, List, Iterator, TypeVar, Union
from functools import reduce

import pandas as pd

from .utils import Filter, Stream, Translation, removeprefix, make_equal_length

from . import SPECIES


IndexedData = Union[pd.DataFrame, pd.Series]
T = TypeVar('T')


def chain(funcs: List[Callable[[T], T]]) -> Callable[[T], T]:
    return reduce(lambda f, g: lambda x: f(g(x)), funcs)


class _Operator_meta(type):
    def __new__(cls, name, bases, classdict):
        result = type.__new__(cls, name, bases, dict(classdict))
        if bases:
            result._members = [x for x in classdict
                if (not x.startswith('_') and not callable(classdict[x]))]
        return result


class Operator(metaclass=_Operator_meta):
    def __init__(self, *args, **kwargs):
        members = self._members.copy()
        for arg in args:
            if not members:
                raise TypeError(f'Too many arguments')
            setattr(self, members.pop(0), arg)
        for kwarg in kwargs:
            if kwarg not in members:
                raise TypeError(f'Unexpected keyword: {kwarg}')
            setattr(self, kwarg, kwargs[kwarg])

    def __call__(self, series: pd.Series) -> pd.Series:
        return self.op(series)

    def op(self, series: pd.Series) -> pd.Series:
        raise NotImplementedError()

    @property
    def to_filter(self) -> Filter:
        def _filter(stream: Stream) -> Stream:
            for series in stream:
                yield self(series)
        return _filter


class InvalidSeries(Exception):
    def __init__(self, series: pd.Series, error: str):
        self.series = series
        super().__init__(error)


class OpCheckValid(Operator):
    def op(self, series: pd.Series) -> pd.Series:
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


class OpIndexToMulti(Operator):
    def op(self, series: pd.Series) -> pd.Series:
        if not isinstance(series.index, pd.MultiIndex):
            series.index = pd.MultiIndex.from_arrays(
                [series.index], names=[series.index.name])
        return series


class OpSpeciesToFront(Operator):
    def op(self, series: pd.Series) -> pd.Series:
        ordered = list(series.index.names)
        ordered.remove(SPECIES)
        ordered = [SPECIES] + ordered
        return series.reorder_levels(ordered)


class OpIndexMerge(Operator):
    glue: str = '_'
    def op(self, series: pd.Series) -> pd.Series:
        series.index = series.index.to_frame().apply(self.glue.join, axis=1)
        return series


class OpIndexSpeciesOnly(Operator):
    def op(self, series: pd.Series) -> pd.Series:
        series.index = series.index.to_frame()[SPECIES]
        return series


class OpDropEmpty(Operator):
    def op(self, series: pd.Series) -> pd.Series:
        series.dropna(inplace=True)
        return series[series.str.len() > 0]


class OpPadRight(Operator):
    fill: str = '-'
    def op(self, series: pd.Series) -> pd.Series:
        return make_equal_length(series.dropna(), fillchar=self.fill)


class OpTranslate(Operator):
    translation: Translation = {}
    def op(self, series: pd.Series) -> pd.Series:
        return series.str.translate(self.translation)


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
    return OpSpeciesToFront()(OpIndexToMulti()(all))


def join_any_to_stream(stream: Stream) -> Stream:
    for series in join_any(stream):
        yield series