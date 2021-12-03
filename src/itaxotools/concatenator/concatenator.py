#!/usr/bin/env python

import tkinter as tk

from .library.concat_gui import ConcatGUI


def gui_main() -> None:
    root = tk.Tk()

    def close_window():
        root.destroy()
        root.quit()

    root.title("Concatenator")

    root.protocol("WM_DELETE_WINDOW", close_window)
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    ConcatGUI(root)

    root.mainloop()
    root.quit()


def main() -> None:
    gui_main()


if __name__ == "__main__":
    main()
