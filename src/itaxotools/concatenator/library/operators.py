
from typing import Callable, Dict, List, Iterator, TypeVar, Set, Union
from functools import reduce

import pandas as pd

from .utils import (
    Stream, Filter, MultiFilter, Translation, ConfigurableCallable, Param,
    removeprefix, make_equal_length,
    )

from . import SPECIES


IndexedData = Union[pd.DataFrame, pd.Series]
T = TypeVar('T')


def chain(funcs: List[Callable[[T], T]]) -> Callable[[T], T]:
    return reduce(lambda f, g: lambda x: f(g(x)), funcs)


class Operator(ConfigurableCallable):
    @property
    def to_filter(self) -> Filter:
        def _filter(stream: Stream) -> Stream:
            for series in stream:
                result = self(series)
                if result is not None:
                    yield result
        return _filter

    @property
    def multi_filter(self) -> MultiFilter:
        def _filter(streams: Iterator[Stream]) -> Stream:
            for stream in streams:
                for series in stream:
                    result = self(series)
                    if result is not None:
                        yield result
        return _filter


class InvalidSeries(Exception):
    def __init__(self, series: pd.Series, error: str):
        self.series = series
        super().__init__(error)


class OpCheckValid(Operator):
    def call(self, series: pd.Series) -> pd.Series:
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
    def call(self, series: pd.Series) -> pd.Series:
        if not isinstance(series.index, pd.MultiIndex):
            series = series.copy(deep=False)
            series.index = pd.MultiIndex.from_arrays(
                [series.index], names=[series.index.name])
        return series


class OpSpeciesToFront(Operator):
    def call(self, series: pd.Series) -> pd.Series:
        ordered = list(series.index.names)
        if SPECIES not in ordered:
            return series
        ordered.remove(SPECIES)
        ordered = [SPECIES] + ordered
        series = series.copy(deep=False)
        return series.reorder_levels(ordered)


class OpIndexMerge(Operator):
    glue: str = Param('_')

    def call(self, series: pd.Series) -> pd.Series:
        series = series.copy(deep=False)
        indices = series.index.to_frame().fillna('', inplace=False)
        series.index = pd.Index(
            indices.apply(self.glue.join, axis=1),
            name='merged')
        return series


class OpIndexSpeciesOnly(Operator):
    def call(self, series: pd.Series) -> pd.Series:
        series = series.copy(deep=False)
        series.index = series.index.to_frame()[SPECIES]
        return series


class OpDropEmpty(Operator):
    missing: str = Param('')

    def call(self, series: pd.Series) -> pd.Series:
        series = series.dropna(inplace=False)
        if self.missing:
            series = series[~ series.str.fullmatch(f'[{self.missing}]+')]
        return series[series.str.len() > 0]


class OpPadRight(Operator):
    padding: str = Param('-')

    def call(self, series: pd.Series) -> pd.Series:
        if not self.padding:
            return series
        series = series.fillna('', inplace=False)
        return make_equal_length(series, fillchar=self.padding)


class OpTranslateSequences(Operator):
    translation: Translation = Param({})

    def call(self, series: pd.Series) -> pd.Series:
        return series.str.translate(self.translation)


class OpFilterCharsets(Operator):
    filter: Set = Param(set())

    def call(self, series: pd.Series) -> pd.Series:
        if series.name not in self.filter:
            return None
        return series


class OpTranslateCharsets(Operator):
    translation: Dict = Param(dict())

    def call(self, series: pd.Series) -> pd.Series:
        if series.name in self.translation:
            if self.translation[series.name] is None:
                return None
            series = series.copy(deep=False)
            series.name = self.translation[series.name]
        return series


class OpChainCharsets(Operator):
    allow_duplicates: bool = Param(False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._memory = set()

    def call(self, series: pd.Series) -> pd.Series:
        if self.allow_duplicates:
            return series
        if series.name in self._memory:
            return None
        self._memory.add(series.name)
        return series


class OpApply(Operator):
    func: Callable[[pd.Series], pd.Series] = Param(None)

    def call(self, series: pd.Series) -> pd.Series:
        if self.func is None:
            return series
        return self.func(series)


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
            keys = guard(list(series.index.names))
            series.index = pd.MultiIndex.from_frame(
                series.index.to_frame(), names=keys)
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
