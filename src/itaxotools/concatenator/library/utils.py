#!/usr/bin/env python

from typing import Callable, Dict, Iterator, Optional, Any

import pandas as pd


Stream = Iterator[pd.Series]
Filter = Callable[[Stream], Stream]
MultiFilter = Callable[[Iterator[Stream]], Stream]

# Such as returned by str.maketrans
Translation = Dict[int, int]


class Param:
    def __init__(self, default: Any = None):
        self.default = default


class _ConfigurableCallable_meta(type):
    def __new__(cls, name, bases, classdict):
        result = type.__new__(cls, name, bases, classdict)
        result._params_ = [
            x for x in classdict if isinstance(classdict[x], Param)]
        for param in result._params_:
            setattr(result, param, classdict[param].default)
        return result


class ConfigurableCallable(metaclass=_ConfigurableCallable_meta):
    def __init__(self, *args, **kwargs):
        members = self._params_.copy()
        for arg in args:
            if not members:
                raise TypeError('Too many arguments')
            setattr(self, members.pop(0), arg)
        for kwarg in kwargs:
            if kwarg not in members:
                raise TypeError(f'Unexpected keyword: {kwarg}')
            setattr(self, kwarg, kwargs[kwarg])

    def __call__(self, *args, **kwargs) -> Any:
        return self.call(*args, **kwargs)

    def call(self, *args, **kwargs) -> Any:
        raise NotImplementedError()


# For Python 3.8 compatibility
def removeprefix(text: str, prefix: str) -> str:
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def into_seqids(table: pd.DataFrame) -> pd.Series:
    """
    Concatenates a collection of string columns into a column of seqids.

    Removes unallowed characters.
    """
    for column in table.columns:
        # Remove unallowed character in the beginning
        table[column] = table[column].str.replace(r"^[^A-Za-z0-9]+", "", n=1, regex=True)
        # Remove unallowed character in the end
        table[column] = table[column].str.replace(r"[^A-Za-z0-9]+$", "", n=1, regex=True)
        # Replace substrings of unallowed characters in the middle with _
        table[column] = table[column].str.replace(r"[^A-Za-z0-9]+", "_", regex=True)
    return table.iloc[:, 0].str.cat(table.iloc[:, 1:], sep="_")


def max_length_if_not_same(column: pd.Series) -> Optional[int]:
    """
    Returns length of the longest string in the ``column``,
    if the strings have different length
    """
    max_length = column.str.len().max()
    min_length = column.str.len().min()
    if max_length > min_length:
        return max_length
    else:
        return None


def make_equal_length(
    column: pd.Series, max_length: Optional[int] = None, fillchar: str = "N"
) -> pd.Series:
    """
    Pads strings in the ``column`` to ``max_length``.

    If ``max_length`` is None, uses the length of the longest string
    """
    if not max_length:
        max_length = column.str.len().max()
    if max_length:
        return column.str.ljust(max_length, fillchar)
    else:
        return column


def has_uniform_length(series: pd.Series) -> bool:
    max_length = series.str.len().max()
    min_length = series.str.len().min()
    return (max_length == min_length)
