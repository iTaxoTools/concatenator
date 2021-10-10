
from typing import Callable, Dict, Iterator, List
from functools import reduce

import pandas as pd


Stream = Iterator[pd.Series]
Filter = Callable[[Stream], Stream]

# Such as returned by str.maketrans
Translation = Dict[int, int]


def translate(translation: Translation) -> Filter:
    def _translate(stream: Stream) -> Stream:
        for series in stream:
            yield series.str.translate(translation)
    return _translate


def join(stream: Stream) -> Stream:
    all = pd.concat(stream, axis=1)
    for col in all:
        yield all[col].fillna('')


def chain(filters: List[Filter]) -> Filter:
    return reduce(lambda f, g: lambda x: f(g(x)), filters)
