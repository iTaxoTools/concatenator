#!/usr/bin/env python3

import os
from pathlib import Path
from datetime import datetime
from typing import Callable, Tuple, List
import logging

import tkinter as tk
import tkinter.ttk as ttk
import tkinter.filedialog as tkfiledialog
import tkinter.messagebox as tkmessagebox

from .operations import Operation, Parameter, run_pipeline
from .file_identify import autodetect
from .file_types import FileFormat


class Operator(ttk.Frame):
    def __init__(
        self,
        *args,
        operation: Operation,
        operations: List[Tuple[Operation, Parameter]],
        command=Callable[[], None],
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.operations = operations
        self.operation = operation
        self.parameter = tk.StringVar()
        self.command = command
        ttk.Button(
            self, text=self.operation.button_text, command=self.push_operation
        ).grid(row=0, column=0, sticky="w")
        if self.operation.parameter_type == int or self.operation.parameter_type == str:
            ttk.Entry(self, textvariable=self.parameter).grid(
                row=0, column=1, sticky="w"
            )
        elif isinstance(self.operation.parameter_type, list):
            ttk.Combobox(
                self,
                textvariable=self.parameter,
                values=self.operation.parameter_type,
                state="readonly",
            ).grid(row=0, column=1, sticky="w")
            self.parameter.set(self.operation.parameter_type[0])
        elif self.operation.parameter_type == None:
            pass
        else:
            assert False

    def push_operation(self) -> None:
        if self.operation.parameter_type == int:
            try:
                parameter = int(self.parameter.get())
            except ValueError:
                parameter = 0
        else:
            parameter = self.parameter.get()
        if self.operations:
            last_op, last_param = self.operations[-1]
            last_output_type = last_op.output_type(last_param)
            if last_op.output_type(last_param) != self.operation.input_type:
                self.show_type_error(
                    last_op.description.format(last_param), last_output_type
                )
                return
        self.operations.append((self.operation, parameter))
        self.command()

    def show_type_error(self, last_op: str, last_output_type: FileFormat) -> None:
        message = "\n".join(
            [
                "Can't add operation:",
                self.operation.description.format(self.parameter.get()),
                "",
                "The last operation:",
                last_op,
                "outputs:",
                last_output_type.description,
                "But the new operation:",
                self.operation.description.format(self.parameter.get()),
                "requires:",
                self.operation.input_type.description,
            ]
        )
        tkmessagebox.showerror("Error", message)


class TkWarnLogger(logging.Handler):
    """docstring for TkWarnLogger."""

    def __init__(self, level=logging.NOTSET):
        logging.Handler.__init__(self, level)
        self.addFilter(lambda record: record.levelno == logging.WARNING)

    def emit(self, record: logging.LogRecord) -> None:
        tkmessagebox.showwarning("Warning", record.getMessage())
        print(record.pathname, record.lineno, sep=": ")
        print("Warning:", record.getMessage(), "\n")


class ConcatGUI(ttk.Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.input_file = tk.StringVar()
        self.output_dir = tk.StringVar()

        self.operation = tk.StringVar(value="concat")

        self.operations: List[Tuple[Operation, Parameter]] = []
        self.operations_display = tk.StringVar()

        self.make_button_frame()
        self.make_entries_frame()
        self.make_operations_frame()

        self.rowconfigure(2, weight=1)
        self.columnconfigure(0, weight=1)

        self.grid(row=0, column=0, sticky="nsew")
        logger = logging.getLogger()
        logger.addHandler(TkWarnLogger())
        logger.setLevel(logging.WARNING)

    def make_button_frame(self) -> None:
        button_frame = ttk.Frame(self)

        current_column = 0

        ttk.Button(
            button_frame, text="Browse input file", command=self.browse_input
        ).grid(row=0, column=current_column)
        current_column += 1

        ttk.Button(
            button_frame, text="Browse output directory", command=self.browse_output
        ).grid(row=0, column=current_column)
        current_column += 1

        ttk.Button(button_frame, text="Execute", command=self.run).grid(
            row=0, column=current_column
        )
        current_column += 1

        button_frame.columnconfigure(current_column, weight=1)

        button_frame.grid(row=0, column=0, sticky="nsew")

    def make_entries_frame(self) -> None:
        entries_frame = ttk.Frame(self)

        current_row = 0

        ttk.Label(entries_frame, text="Input file").grid(
            row=current_row, column=0, sticky="w"
        )
        ttk.Entry(entries_frame, textvariable=self.input_file).grid(
            row=current_row, column=1, sticky="we"
        )
        current_row += 1

        ttk.Label(entries_frame, text="Output directory").grid(
            row=current_row, column=0, sticky="w"
        )
        ttk.Entry(entries_frame, textvariable=self.output_dir).grid(
            row=current_row, column=1, sticky="we"
        )
        current_row += 1

        entries_frame.rowconfigure(current_row, weight=1)
        entries_frame.columnconfigure(1, weight=1)

        entries_frame.grid(row=1, column=0, sticky="nsew")

    def make_operations_frame(self) -> None:
        operations_frame = ttk.Frame(self)
        list_frame = ttk.LabelFrame(operations_frame, text="Operations")

        current_row = 0

        for operation in Operation:
            Operator(
                list_frame,
                operation=operation,
                operations=self.operations,
                command=self.display_operations,
            ).grid(row=current_row, column=0, sticky="w")
            current_row += 1

        ttk.Button(list_frame, text="Undo", command=self.pop_operation).grid(
            row=current_row, column=0, sticky="w"
        )
        current_row += 1

        list_frame.rowconfigure(current_row, weight=1)
        list_frame.grid(row=0, column=0, sticky="nsew")

        ttk.Label(operations_frame, textvariable=self.operations_display).grid(
            row=0, column=1, sticky="nw"
        )
        operations_frame.rowconfigure(0, weight=1)
        operations_frame.columnconfigure(1, weight=1)

        operations_frame.grid(row=2, column=0, sticky="nsew")

    def pop_operation(self) -> None:
        self.operations.pop()
        self.display_operations()

    def browse_input(self) -> None:
        newpath = tkfiledialog.askopenfilename()
        if newpath:
            try:
                autodetect(Path(newpath))
            except Exception as e:
                tkmessagebox.showerror("Error", str(e))
                raise
            self.input_file.set(os.path.abspath(newpath))

    def browse_output(self) -> None:
        newpath = tkfiledialog.askdirectory()
        if newpath:
            self.output_dir.set(os.path.abspath(newpath))

    def generate_output_suffix(self) -> str:
        if self.operations:
            last_op, last_param = self.operations[-1]
            output_type = last_op.output_type(last_param)
            ext = output_type.extension
            if output_type.timestamp:
                return "_" + datetime.utcnow().strftime("%Y%m%dT%H%M%S") + ext
            else:
                return ext
        else:
            _, ext = os.path.splitext(self.input_file.get())
            return ext

    def run(self) -> None:
        filename, _ = os.path.splitext(os.path.basename(self.input_file.get()))
        suffix = self.generate_output_suffix()
        output_file = os.path.join(self.output_dir.get(), filename + suffix)
        try:
            with open(self.input_file.get(), mode="rb") as infile, open(
                output_file, mode="wb"
            ) as outfile:
                run_pipeline(infile, outfile, self.operations)
            tkmessagebox.showinfo("Done", "Processing is complete")
        except FileNotFoundError as e:
            tkmessagebox.showerror("Error", f"File {e.filename} not found")
            raise
        except Exception as e:
            tkmessagebox.showerror("Error", str(e))
            raise

    def display_operations(self):
        self.operations_display.set(
            "\n".join(
                operation.description.format(parameter)
                for operation, parameter in self.operations
            )
        )
