#!/usr/bin/env python

from typing import TextIO, List, TypeVar

T = TypeVar('T')


def pop_many(l: List[T], indices: List[int]) -> List[T]:
    """
    Pops the elements of `l` at `indices` and returns the list of them.

    Sorts `indices` in reverse order
    """
    indices.sort()
    indices.reverse()
    result = [l.pop(i) for i in indices]
    result.reverse()
    return result


def process(input: TextIO, output: TextIO) -> None:
    # read the header
    header = input.readline().rstrip().split("\t")

    # find the columns with description (i.e. not sequences) values
    description_indices = [i for i, column in enumerate(
        header) if "sequence" not in column.lower()]
    description_columns = pop_many(header, description_indices)

    # print the header for the output
    print(*description_columns, "sequence", sep="\t", file=output)
    for line in input:
        # separate description and sequences
        values = line.split("\t")
        description_values = pop_many(values, description_indices)

        # print the description
        print(*description_values, sep="\t", end="\t", file=output)
        # print the sequence
        print(*values, sep="", file=output)
