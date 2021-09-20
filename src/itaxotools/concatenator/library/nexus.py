#!/usr/bin/env python

from enum import Enum
from typing import TextIO, Iterator, ClassVar, Set, List, Optional, Tuple, Dict
import tempfile
import re

import pandas as pd

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

        gene_name = gene.columns[-1][len("sequence_"):]
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


def read(input: TextIO) -> pd.DataFrame:
    commands = NexusCommands(input)
    reader = NexusReader()
    for command, args in commands:
        reader.execute(command, args)
    return reader.return_table()


class Tokenizer:
    """
    Token iterator for the NEXUS format

    Emits the stream of words and punctuation
    """

    punctuation: ClassVar[Set[str]] = set("=;")

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
                    "".join(self.token)
                else:
                    raise StopIteration
            elif c in Tokenizer.punctuation:
                # punctuation is a token by itself => save it into the token
                token = self.replace_token([c])
                if token:
                    return token
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
        self.columns = ["seqid"]
        self.state: Optional[NexusState] = None
        self.todo: Set[NexusState] = {NexusState.Data, NexusState.Sets}
        self.ntax: Optional[int] = None
        self.read_matrix = False

    def begin_block(self, args: Iterator[str]) -> None:
        """
        Sets the state for the next block
        """
        try:
            arg = next(args)
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

    def configure_format(self, args: Iterator[str]) -> None:
        """
        Configure the internal state in responce to the 'format' command
        """
        # If the 'datatype' is DNA, RNA, Nucleotide or Protein, prepare for reading
        for arg in args:
            if arg.casefold() == "datatype":
                if next(args) != "=":
                    continue
                if re.search(r"DNA|RNA|Nucleotide", next(args), flags=re.IGNORECASE):
                    self.read_matrix = True
            elif arg.casefold() == "interleave":
                self.interleave = True

    def read_dimensions(self, args: Iterator[str]) -> None:
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
        if self.read_matrix and self.ntax and self.state == NexusState.Data:
            column_count = 0
            while True:
                try:
                    column = {}
                    for _ in range(self.ntax):
                        seqid = next(args)
                        gene = next(args)
                        column[seqid] = gene
                except StopIteration:
                    break
                self.table[column_count] = pd.Series(column)
                column_count += 1
        self.read_matrix == False
        self.ntax == None

    def add_charset(self, args: Iterator[str]) -> None:
        """
        Adds charset name to self.columns
        """
        if self.state == NexusState.Sets:
            try:
                self.columns.append("sequence_" + next(args))
            except StopIteration:
                self.columns.append("")

    def return_table(self) -> pd.DataFrame:
        self.table.reset_index(inplace=True)
        self.table.columns = self.columns
        return self.table