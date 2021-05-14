#!/usr/bin/env python

import tkinter as tk

from library.gui import ConcatGUI


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


if __name__ == "__main__":
    gui_main()
