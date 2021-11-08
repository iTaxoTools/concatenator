#!/usr/bin/env python

from typing import Callable, Dict, Iterator, Optional, Any
from itertools import chain
from enum import Enum

import pandas as pd

from itaxotools.common.param.core import Field, Group


# Such as returned by str.maketrans
Translation = Dict[int, int]


class _ConfigurableCallable_meta(type):
    def __new__(cls, name, bases, classdict):
        new_params = {
            param: classdict[param] for param in classdict
            if isinstance(classdict[param], Field)}
        for param in new_params:
            classdict.pop(param)
        obj = super().__new__(cls, name, bases, classdict)
        inherited_params = dict()
        if hasattr(obj, '_class_params_'):
            inherited_params = obj._class_params_.copy()
        obj._class_params_ = dict(**inherited_params, **new_params)
        return obj


class ConfigurableCallable(metaclass=_ConfigurableCallable_meta):
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def update(self, *args, **kwargs):
        self._params_ = {
            param: self._class_params_[param].copy()
            for param in self._class_params_}
        keys = list(self._params_.keys())
        for arg in args:
            if not keys:
                raise TypeError('Too many arguments')
            self._params_[keys.pop(0)].value = arg
        for kwarg in kwargs:
            if kwarg not in keys:
                raise TypeError(f'Unexpected keyword: {kwarg}')
            self._params_[kwarg].value = kwargs[kwarg]
        return self

    @property
    def params(self) -> Group:
        return Group(key='root', children=[
            param for param in self._params_.values()])

    def __getattr__(self, attr):
        if attr in self._params_:
            return self._params_[attr].value

    def __call__(self, *args, **kwargs) -> Any:
        return self.call(*args, **kwargs)

    def call(self, *args, **kwargs) -> Any:
        raise NotImplementedError()


class OrderedSet(dict):
    def __init__(self, iterator: Iterator = {}):
        super().__init__()
        self.update(iterator)

    def __and__(self, other):
        return OrderedSet(key for key in self if key in other)

    def __or__(self, other):
        return OrderedSet(chain((key for key in self), (key for key in other)))

    def update(self, iterator: Iterator):
        super().update({key: None for key in iterator})

    def add(self, item):
        super().update({item: None})


# For Python 3.8 compatibility
def removeprefix(text: str, prefix: str) -> str:
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def reverse_complement(nucleotides: str) -> str:
    reversed = nucleotides[::-1]
    translation = str.maketrans('ACGTacgt', 'TGCAtgca')
    return reversed.translate(translation)


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


def fill_empty(column: pd.Series, filler: str) -> pd.Series:
    column = column.fillna('')
    max_length = column.str.len().max()
    column = column.str.ljust(max_length, filler)
    return column


def has_uniform_length(series: pd.Series) -> bool:
    max_length = series.str.len().max()
    min_length = series.str.len().min()
    return (max_length == min_length)
