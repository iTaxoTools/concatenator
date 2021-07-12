#!/usr/bin/env python3

from collections import defaultdict
from typing import Dict, Set, Optional, List
import json
import random
from library.gencode_detection import detect_reading_frame

bases = "TCAG"

codons = [c1 + c2 + c3 for c1 in bases for c2 in bases for c3 in bases]


def load_tables() -> Optional[Dict[int, Set[int]]]:
    """
    Reads data/stop_codons.json and returns mapping from translation table to indices of stop codons
    """
    try:
        stop_codons: Dict[str, List[int]] = json.load(open("data/stop_codons.json"))
    except json.JSONDecodeError:
        return None
    result = defaultdict(set)
    for stop_codon, tables in stop_codons.items():
        codon_i = codons.index(stop_codon)
        for table in tables:
            result[table].add(codon_i)
    return result


def random_codon(tables: Dict[int, Set[int]], table: int) -> int:
    """
    Returns a random codon, which is not a stop codon for the table, according to tables
    """
    while True:
        codon = random.randrange(0, 64)
        if codon not in tables[table]:
            return codon


def random_sequence(tables: Dict[int, Set[int]], table: int, length: int) -> str:
    """
    Returns a random genetic sequence that doesn't contain a stop codon from the table or contain a stop codon at the end, according to tables
    """
    seq = "".join(codons[random_codon(tables, table)] for _ in range(length))
    if random.choice([True, False]):
        stopcodon = random.choice(tuple(tables[table]))
        seq = seq[:-3] + codons[stopcodon]
    return seq


def shift_frame(seq: str, frame: int) -> str:
    shift_width = abs(frame) - 1
    seq = seq[:shift_width] + seq
    if frame < 0:
        return seq[::-1]
    else:
        return seq


def test_table_detection() -> None:
    tables = load_tables()
    assert tables is not None
    tables_list = list(tables.keys())
    print()
    for _ in range(16):
        table = random.choice(tables_list)
        frame = random.choice([1, -1, 2, -2, 3, -3])
        print(
            f"Generating sequence for translation table {table} with reading frame {frame}..."
        )
        seq = shift_frame(random_sequence(tables, table, 100), frame)
        print(seq)
        print()
        print("Trying to detect the reading frame...")
        detected_frames = detect_reading_frame(seq)
        if detected_frames:
            if len(detected_frames) == 1:
                detected_frame = detected_frames[0]
                assert detected_frame == frame
                print("Detection successful")
            else:
                print("Non-unique reading frame: ", detected_frames)
        else:
            print("Cannot detect the reading frame")
        print()
        print()
