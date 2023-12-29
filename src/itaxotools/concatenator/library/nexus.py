#!/usr/bin/env python

from typing import TextIO, Iterator, ClassVar, Set, List, Optional, Tuple, Dict
from enum import Enum
import tempfile
import re
import logging
from itertools import chain

import pandas as pd

from .model import GeneStream, GeneDataFrame, PathLike
from .types import Justification, Charset
from .utils import *


def write(genes: Iterator[pd.DataFrame], output: TextIO) -> None:
    buf = tempfile.TemporaryFile(mode="w+")
    charsets = {}
    ntax = 0

    for gene in genes:
        gene.iloc[:, -1] = make_equal_length(gene.iloc[:, -1])
        output_fragment = pd.DataFrame()
        output_fragment["seqid"] = into_seqids(gene.iloc[:, :-1].copy())
        output_fragment["sequence"] = gene.iloc[:, -1]

        gene_name = gene.columns[-1][len("sequence_") :]
        gene_len = len(gene.iat[0, -1])
        charsets[gene_name] = gene_len

        output_fragment.to_string(buf, header=False, index=False)
        print("\n", file=buf)

        ntax = len(output_fragment)

    buf.seek(0)

    nchar = sum(charsets.values())

    print("#NEXUS\n", file=output)
    print("begin data;\n", file=output)
    print(
        "format datatype=DNA missing=N missing=? Gap=- Interleave=yes;\n", file=output
    )
    print(f"dimensions Nchar={nchar} Ntax={ntax};\n", file=output)
    print("matrix\n", file=output)

    for line in buf:
        output.write(line)
    print(";\nend;", file=output)

    output.write("\n\n\n")

    print("begin sets;\n", file=output)

    gene_position = 1

    for gene_name, gene_len in charsets.items():
        print(
            "charset ",
            gene_name,
            " = ",
            gene_position,
            "-",
            gene_position + gene_len - 1,
            ";",
            sep="",
            file=output,
        )
        gene_position += gene_len

    print("\nend;", file=output)
    buf.close()


def nexus_writer(
    stream: GeneStream,
    out: TextIO,
    justification: Justification = Justification.Left,
    separator: str = " ",
    include_full_markers_with_codons: bool = False,
) -> None:
    buffer = tempfile.TemporaryFile(
        mode="w+", encoding="utf-8", errors="surrogateescape"
    )
    charsets = list()
    missings = OrderedSet()
    gaps = OrderedSet()
    ntax = 0
    cursor = 1

    for gene in stream:
        series = gene.series
        assert has_uniform_length(series)

        length = len(series.iat[0])
        charsets.append(
            Charset(gene.name, cursor, length, gene.reading_frame, gene.codon_names)
        )
        cursor += length
        index_len = series.index.str.len().max()
        for index, sequence in series.items():
            buffer.write(
                (f"{justification.apply(index, index_len)}" f"{separator}{sequence}\n")
            )
        buffer.write("\n")
        ntax = len(series)
        missings |= OrderedSet(gene.missing.upper())
        gaps |= OrderedSet(gene.gap.upper())

    out.write("#NEXUS\n\n")
    out.write("BEGIN DATA;\n\n")
    out.write(f"Dimensions Nchar={cursor-1} Ntax={ntax};\n")
    out.write("Format Datatype=DNA")
    for missing in sorted(missings, reverse=True):
        out.write(f" Missing={missing}")
    for gap in sorted(gaps, reverse=True):
        out.write(f" Gap={gap}")
    out.write(" Interleave=yes;\n")
    out.write("Matrix\n\n")

    buffer.seek(0)
    for line in buffer:
        out.write(line)
    out.write(";\nEND;\n\n\n")
    out.write("BEGIN SETS;\n\n")
    buffer.close()

    if include_full_markers_with_codons:
        for charset in charsets:
            out.write(f'charset {str(charset)}\n')
        out.write('\n')
        for charset in (cs for cs in charsets if cs.frame):
            for codon in charset.codon_sets():
                out.write(f'charset {str(codon)}\n')
            out.write('\n')
    else:
        full_charsets = [cs for cs in charsets if not cs.frame]
        split_charsets = [cs for cs in charsets if cs.frame]

        if full_charsets:
            for charset in full_charsets:
                out.write(f'charset {str(charset)}\n')
            out.write('\n')

        if split_charsets:
            for charset in split_charsets:
                for codon in charset.codon_sets():
                    out.write(f'charset {str(codon)}\n')
                out.write('\n')

    out.write("END;\n")


def stream_to_path(stream: GeneStream, path: PathLike, *args, **kwargs) -> None:
    with path.open("w", encoding="utf-8", errors="surrogateescape") as file:
        nexus_writer(stream, file, *args, **kwargs)


def read(input: TextIO) -> pd.DataFrame:
    commands = NexusCommands(input)
    reader = NexusReader()
    for command, args in commands:
        reader.execute(command, args)
    return reader.return_table()


def dataframe_from_path(path: PathLike) -> GeneDataFrame:
    with path.open(encoding="utf-8", errors="surrogateescape") as file:
        commands = NexusCommands(file)
        reader = NexusReader()
        for command, args in commands:
            reader.execute(command, args)
        data = reader.return_table()
    data.set_index("seqid", inplace=True)
    missing = "".join(reader.missings)
    gap = "".join(reader.gaps)
    gdf = GeneDataFrame(data, missing=missing, gap=gap)
    return gdf


def stream_from_path(path: PathLike) -> GeneStream:
    gdf = dataframe_from_path(path)
    return gdf.to_stream()


class Tokenizer:
    """
    Token iterator for the NEXUS format

    Emits the stream of words and punctuation
    """

    punctuation: ClassVar[Set[str]] = set("=;\n")

    def __init__(self, file: TextIO):
        """
        Iterate over the token in 'file'
        """
        # check that the file is in NEXUS format
        magic_word = file.read(6)
        if magic_word != "#NEXUS":
            raise ValueError("The input file is not a nexus file")
        # contains the underlying file
        self.file = file
        # the currently read token
        self.token: List[str] = []
        # the currently read line
        self.line = ""
        # the reading position in the line
        self.line_pos = 0

    def peek_char(self) -> Optional[str]:
        """
        Returns the next char, without advancing the position.
        Returns None, if EOF is reached
        """
        try:
            c = self.line[self.line_pos]
            return c
        except IndexError:
            # line_pos is self.line.len()
            # it's equivalent to 0 in the next line
            self.line = self.file.readline()
            if self.line == "":
                # EOF is reached
                return None
            self.line_pos = 0
            c = self.line[0]
            return c

    def get_char(self) -> Optional[str]:
        """
        Emits a char, advancing the reading
        Returns None, if EOF is reached
        """
        c = self.peek_char()
        self.line_pos += 1
        return c

    def replace_token(self, token: List[str]) -> str:
        """
        Remember the next token and return the str representation of the current one
        """
        self.token, token = token, self.token
        return "".join(token)

    def skip_comment(self) -> None:
        """
        Advance the iterator past the end of a comment
        """
        while True:
            c = self.get_char()
            if c is None:
                # EOF is reached, should not happen in a well-formed file
                raise ValueError("Nexus: EOF inside a comment")
            elif c == "[":
                # comment inside a comment
                self.skip_comment
            elif c == "]":
                # end of the comment
                break

    def read_quoted(self) -> List[str]:
        """
        Reads a quoted string as one token
        """
        s = []
        while True:
            c = self.get_char()
            if c is None:
                # EOF is reached, should not happen in a well-formed file
                raise ValueError("Nexus: EOF inside a quoted value")
            elif c == "'":
                # possible end of the quoted value
                if self.peek_char == "'":
                    # '' is ', and not the end
                    s += ["'"]
                else:
                    # the end of the line
                    return s
            else:
                # update the line
                s += [c]

    def __iter__(self) -> "Tokenizer":
        """Tokenizer is an Iterator"""
        return self

    def __next__(self) -> str:
        if self.token:
            # the is a previously saved token
            return self.replace_token([])
        while True:
            c = self.get_char()
            if c is None:
                # EOF => return the last token
                if self.token:
                    token = self.replace_token([])
                    return token
                else:
                    raise StopIteration
            elif c in Tokenizer.punctuation:
                # punctuation is a token by itself => save it into the token
                if self.token:
                    token = self.replace_token([c])
                    return token
                return c
            elif c == "[":
                # a comment => skip it
                self.skip_comment()
            elif c == "'":
                # a quoted value => read it and save into the token
                token = self.replace_token(self.read_quoted())
                if token:
                    return token
            elif c.isspace():
                # whitespace => return the token, if it's the first whitespace
                if self.token:
                    token = self.replace_token([])
                    return token
            else:
                # otherwise => update the token
                self.token.append(c)

    @staticmethod
    def print_tokens(path: str) -> None:
        """
        Testing method to check that the tokens are read correctly
        """
        with open(path) as file:
            for token in Tokenizer(file):
                print(repr(token))


class NexusCommands:
    """
    Iterator that emits NEXUS command as a tuple of the command name and the arguments' iterator
    """

    def __init__(self, file: TextIO):
        """
        Iterate over the commands in 'file'
        """
        self.tokenizer = Tokenizer(file)

    def __iter__(self) -> "NexusCommands":
        """
        NexusCommands is an Iterator
        """
        return self

    def __next__(self) -> Tuple[str, Iterator[str]]:
        # next token is the command name
        command = '\n'
        while command == '\n':
            command = next(self.tokenizer).casefold()

        def arguments() -> Iterator[str]:
            # emit tokens until the ';'
            while True:
                try:
                    arg = next(self.tokenizer)
                except StopIteration:
                    # EOF is reached; should not happen in a well-formed file
                    raise ValueError("Nexus: EOF inside a command")
                if arg == ";":
                    break
                else:
                    yield arg

        return command, arguments()

    @staticmethod
    def print_commands(path: str) -> None:
        """
        Testing method to check that commands are read correctly
        """
        with open(path) as file:
            for command, args in NexusCommands(file):
                print(repr(command), [repr(arg) for arg in args])


NexusState = Enum("NexusState", ["Data", "Sets"])


class NexusReader:
    """
    A virtual machine that executes NEXUS command and constructs a pandas table
    """

    nexus_state: ClassVar[Dict[str, NexusState]] = dict(
        data=NexusState.Data, sets=NexusState.Sets
    )

    def __init__(self) -> None:
        self.table = pd.DataFrame()
        # maps starts of charsets to their names, should only contain gene charsets
        self.charsets: Dict[int, str] = {}
        self.state: Optional[NexusState] = None
        self.todo: Set[NexusState] = {NexusState.Data, NexusState.Sets}
        self.ntax: Optional[int] = None
        self.read_matrix = False
        self.missings = set()
        self.gaps = set()

    @staticmethod
    def iter_concrete_args(args: Iterator[str]) -> str:
        while True:
            try:
                arg = next(args)
                if arg != '\n':
                    yield arg
            except StopIteration:
                return

    @staticmethod
    def iter_all_args(args: Iterator[str]) -> str:
        yield from args

    @staticmethod
    def iter_until_newline(args: Iterator[str]) -> str:
        while True:
            try:
                arg = next(args)
                if arg == '\n':
                    return
                yield arg
            except StopIteration:
                return

    def begin_block(self, args: Iterator[str]) -> None:
        """
        Sets the state for the next block
        """
        args = self.iter_concrete_args(args)
        try:
            arg = next(args).casefold()
            state = NexusReader.nexus_state[arg]
            if state in self.todo:
                self.state = state
        except KeyError:
            pass
        except StopIteration:
            pass

    def complete_block(self) -> None:
        """
        Marks current block as completed
        """
        if self.state and self.state in self.todo:
            if self.state == NexusState.Data and self.table.empty:
                return
            self.todo.remove(self.state)

    def execute(self, command: str, args: Iterator[str]) -> None:
        """
        Executes a NEXUS command and updates the internal state

        Returns an iterator over sequences in a matrix, if the command is 'matrix'
        """
        # for each command execute the corresponding method
        if command == "begin":
            self.begin_block(args)
        elif command == "format":
            self.configure_format(args)
        elif command == "dimensions":
            self.read_dimensions(args)
        elif command == "end" or command == "endblock":
            self.complete_block()
        elif command == "matrix":
            self.sequences(args)
        elif command == "charset":
            self.add_charset(args)
        else:
            logging.warning("Unknown NEXUS command {command}")
        # skip unused arguments
        for _ in args:
            pass

    def configure_format(self, args: Iterator[str]) -> None:
        """
        Configure the internal state in responce to the 'format' command
        """
        args = self.iter_concrete_args(args)
        # If the 'datatype' is DNA, RNA, Nucleotide or Protein, prepare for reading
        for arg in args:
            if arg.casefold() == "datatype":
                if next(args) != "=":
                    continue
                if re.search(r"DNA|RNA|Nucleotide", next(args), flags=re.IGNORECASE):
                    self.read_matrix = True
            elif arg.casefold() == "missing":
                if next(args) != "=":
                    continue
                missing = next(args)
                self.missings |= {missing.upper() + missing.lower()}
            elif arg.casefold() == "gap":
                if next(args) != "=":
                    continue
                gap = next(args)
                self.gaps |= {gap.upper() + gap.lower()}
            elif arg.casefold() == "interleave":
                self.interleave = True

    def read_dimensions(self, args: Iterator[str]) -> None:
        args = self.iter_concrete_args(args)
        for arg in args:
            if arg.casefold() == "ntax":
                if next(args) != "=":
                    continue
                try:
                    self.ntax = int(next(args))
                except ValueError:
                    pass

    def sequences(self, args: Iterator[str]) -> None:
        """
        Reads the matrix into self.table
        """
        args = self.iter_all_args(args)
        if self.read_matrix and self.ntax and self.state == NexusState.Data:

            def get_lines(args):
                while True:
                    try:
                        peek = next(args)
                        if peek == '\n':
                            continue
                    except StopIteration:
                        return
                    yield self.iter_until_newline(chain([peek], args))

            def get_bunches(lines):
                while True:
                    try:
                        peek = iter(list(next(lines)))
                        rest = (next(lines) for _ in range(self.ntax - 1))
                        yield chain([peek], rest)
                    except StopIteration:
                        return

            lines = get_lines(args)
            bunches = get_bunches(lines)

            for column_count, bunch in enumerate(bunches):
                column = {}
                for line in bunch:
                    seqid = next(line)
                    gene = ''.join(line)
                    column[seqid] = gene
                self.table[column_count] = pd.Series(column)

        self.read_matrix = False
        self.ntax = None

    def add_charset(self, args: Iterator[str]) -> None:
        """
        Adds charset name to self.columns
        """
        args = self.iter_concrete_args(args)
        if self.state == NexusState.Sets:
            args_tuple = tuple(args)
            if len(args_tuple) != 3 or args_tuple[1] != "=":
                raise ValueError(f"Invalid arguments for 'charset': {args_tuple}")
            name, _, position = args_tuple
            if "\\" in position:
                return
            start_s, hyphen, _ = position.partition("-")
            if not (hyphen and start_s.isdigit()):
                raise ValueError(f"Invalid charset position: {repr(position)}")
            self.charsets[int(start_s)] = name

    def columns(self) -> List[str]:
        return ["seqid"] + [
            self.charsets[position] for position in sorted(self.charsets.keys())
        ]

    def return_table(self) -> pd.DataFrame:
        self.table.reset_index(inplace=True)
        columns = self.columns()
        if len(self.table.columns) == len(columns):
            self.table.columns = columns
            return self.table
        else:
            logging.warning(
                "Several blocks of sequences in interleaved format detected, "
                "but character sets missing or wrongly specified. "
                "All sequences will be read in as a single gene."
            )
            concat_table = self.table.iloc[:, 0:1].copy()
            if len(self.table.columns) >= 2:
                concat_table["sequence_gene"] = self.table.iloc[:, 1].str.cat(
                    self.table.iloc[:, 2:]
                )
            return concat_table
