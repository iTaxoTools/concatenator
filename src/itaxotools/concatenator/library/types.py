
from typing import Callable, Optional
from enum import Enum


class Justification(Enum):
    Left = 'Left', str.ljust
    Right = 'Right', str.rjust
    Center = 'Center', str.center
    NoJust = 'None', None

    def __init__(self, description: str, method: Optional[Callable]):
        self.description = description
        self.method = method

    def apply(self, text: str, *args, **kwargs):
        if not self.method:
            return text
        return self.method(text, *args, **kwargs)
