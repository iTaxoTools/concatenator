#!/usr/bin/env python3

from typing import Any
import importlib.resources

_resource_path = importlib.resources.files('itaxotools.concatenator') / "resources"


def get_resource(path: Any) -> str:
    return str(_resource_path / path)
