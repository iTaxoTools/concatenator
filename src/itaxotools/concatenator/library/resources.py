#!/usr/bin/env python3

from typing import Any
from pathlib import Path
import sys
import os


def package_path_pyinstaller(package: Any) -> Path:
    if isinstance(package, str):
        parts = package.split('.')
    elif isinstance(package, type(sys)):
        parts = package.__name__.split('.')
    else:
        return None
    path = Path(sys._MEIPASS)
    for part in parts:
        path /= part
    return path


def package_path_file(package: Any) -> Path:
    if isinstance(package, str):
        file = sys.modules[package].__file__
    elif isinstance(package, type(sys)):
        file = package.__file__
    else:
        return None
    path = Path(os.path.dirname(file))
    return path


def package_path_importlib(package: Any) -> Path:
    return importlib.resources.files(package)


try:
    import importlib.resources.files
    package_path = package_path_importlib
except ModuleNotFoundError:
    if hasattr(sys, '_MEIPASS'):
        package_path = package_path_pyinstaller
    else:
        package_path = package_path_file

_resource_path = package_path('itaxotools.concatenator') / 'resources'


def get_resource(path: Any) -> str:
    return str(_resource_path / path)
