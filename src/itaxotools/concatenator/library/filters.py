
from typing import List
from functools import reduce

from .utils import Stream, Filter
from .operators import join_any


def chain(filters: List[Filter]) -> Filter:
    return reduce(lambda f, g: lambda x: f(g(x)), filters)


def join_any_to_stream(stream: Stream) -> Stream:
    for series in join_any(stream):
        yield series
