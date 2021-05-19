#!/usr/bin/env python

import os
from typing import Dict, Callable, TextIO, Tuple

import tkinter as tk
import tkinter.ttk as ttk
import tkinter.filedialog as tkfiledialog
import tkinter.messagebox as tkmessagebox

import library.concat as concat
import library.tab_to_nexus as tab_to_nexus



class ConcatGUI(ttk.Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.input_file = tk.StringVar()
        self.output_dir = tk.StringVar()

        self.operation = tk.StringVar(value="concat")

        self.operations: Dict[str, Tuple[Callable[[TextIO, TextIO], None], str]] = {
                "concat": (concat.process, ".tab"),
                "tab_to_nexus": (tab_to_nexus.process, ".nex"),
                }

        self.make_button_frame()
        self.make_entries_frame()
        self.make_operations_frame()

        self.rowconfigure(2, weight=1)
        self.columnconfigure(0, weight=1)

        self.grid(row=0, column=0, sticky="nsew")

    def make_button_frame(self) -> None:
        button_frame = ttk.Frame(self)

        current_column = 0

        ttk.Button(button_frame, text="Browse input file",
                   command=self.browse_input).grid(row=0, column=current_column)
        current_column += 1

        ttk.Button(button_frame, text="Browse output directory",
                   command=self.browse_output).grid(row=0, column=current_column)
        current_column += 1

        ttk.Button(button_frame, text="Execute", command=self.run).grid(
            row=0, column=current_column)
        current_column += 1

        button_frame.columnconfigure(current_column, weight=1)

        button_frame.grid(row=0, column=0, sticky='nsew')

    def make_entries_frame(self) -> None:
        entries_frame = ttk.Frame(self)

        current_row = 0

        ttk.Label(entries_frame, text="Input file").grid(
            row=current_row, column=0, sticky='w')
        ttk.Entry(entries_frame, textvariable=self.input_file).grid(
            row=current_row, column=1, sticky='we')
        current_row += 1

        ttk.Label(entries_frame, text="Output directory").grid(
            row=current_row, column=0, sticky='w')
        ttk.Entry(entries_frame, textvariable=self.output_dir).grid(
            row=current_row, column=1, sticky='we')
        current_row += 1

        entries_frame.rowconfigure(current_row, weight=1)
        entries_frame.columnconfigure(1, weight=1)

        entries_frame.grid(row=1, column=0, sticky="nsew")

    def make_operations_frame(self) -> None:
        operations_frame = ttk.Frame(self)
        list_frame = ttk.LabelFrame(operations_frame, text="Operations") 

        current_row=0 

        ttk.Radiobutton(list_frame, text="Concatenate Tabfile", variable=self.operation, value="concat").grid(
                row=current_row, column=0, sticky="w")
        current_row += 1

        ttk.Radiobutton(list_frame, text="Tabfile to NEXUS", variable=self.operation, value="tab_to_nexus").grid(
                row=current_row, column=0, sticky="w")
        current_row += 1

        list_frame.rowconfigure(current_row, weight=1)
        list_frame.grid(row=0, column=0, sticky="nsew")

        operations_frame.rowconfigure(0, weight=1)
        operations_frame.columnconfigure(1, weight=1)

        operations_frame.grid(row=2, column=0, sticky="nsew")

    def browse_input(self) -> None:
        newpath = tkfiledialog.askopenfilename()
        if newpath:
            self.input_file.set(os.path.abspath(newpath))

    def browse_output(self) -> None:
        newpath = tkfiledialog.askdirectory()
        if newpath:
            self.output_dir.set(os.path.abspath(newpath))

    def run(self) -> None:
        filename, _ = os.path.splitext(os.path.basename(self.input_file.get()))
        process, out_ext = self.operations[self.operation.get()]
        output_file = os.path.join(self.output_dir.get(), filename + out_ext)
        try:
            with open(self.input_file.get(), errors="replace") as infile, open(output_file, mode="w") as outfile:
                process(infile, outfile)
            tkmessagebox.showinfo("Done", "Processing is complete")
        except FileNotFoundError as e:
            tkmessagebox.showerror("Error", f"File {e.filename} not found")
            raise
        except Exception as e:
            tkmessagebox.showerror("Error", str(e))
            raise
