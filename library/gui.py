#!/usr/bin/env python

import os

import tkinter as tk
import tkinter.ttk as ttk
import tkinter.filedialog as tkfiledialog
import tkinter.messagebox as tkmessagebox

from library.process import process


class ConcatGUI(ttk.Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.input_file = tk.StringVar()
        self.output_dir = tk.StringVar()

        self.make_button_frame()
        self.make_entries_frame()

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

        ttk.Button(button_frame, text="Concatenate", command=self.run).grid(
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

    def browse_input(self) -> None:
        newpath = tkfiledialog.askopenfilename()
        if newpath:
            self.input_file.set(os.path.abspath(newpath))

    def browse_output(self) -> None:
        newpath = tkfiledialog.askdirectory()
        if newpath:
            self.output_dir.set(os.path.abspath(newpath))

    def run(self) -> None:
        filename = os.path.basename(self.input_file.get())
        output_file = os.path.join(self.output_dir.get(), filename)
        try:
            with open(self.input_file.get(), errors="replace") as infile, open(output_file, mode="w") as outfile:
                process(infile, outfile)
            tkmessagebox.showinfo("Done", "Concatenation is complete")
        except FileNotFoundError as e:
            tkmessagebox.showerror("Error", f"File {e.filename} not found")
            raise
        except Exception as e:
            tkmessagebox.showerror("Error", str(e))
            raise
