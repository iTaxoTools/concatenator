
from typing import Callable, List, TypeVar
from functools import reduce

from .utils import Stream
from .operators import join_any


T = TypeVar('T')


def chain(funcs: List[Callable[[T], T]]) -> Callable[[T], T]:
    return reduce(lambda f, g: lambda x: f(g(x)), funcs)


def join_any_to_stream(stream: Stream) -> Stream:
    for series in join_any(stream):
        yield series
