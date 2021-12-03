from typing import List, Tuple, Optional, Dict
from enum import IntEnum, Enum
import random
import itertools
import json

import pandas as pd
import pytest

from itaxotools.concatenator.library.model import GeneSeries
from itaxotools.concatenator.library.codons import (
    GeneticCode, ReadingFrame, NoReadingFrames,
    BadReadingFrame, AmbiguousReadingFrame, last_codon_slice,
    detect_reading_combinations)
import itaxotools.concatenator.library.codons as lib_codons
from itaxotools.concatenator.library.operators import OpDetectReadingFrame


class FrameRejection(IntEnum):
    Ambiguous = 0
    Unambiguous = 1
    Total = 2

    def frames(self):
        frames = {
            FrameRejection.Ambiguous: [-3, -1, 2, 3],
            FrameRejection.Unambiguous: [-3, -2, -1, 2, 3],
            FrameRejection.Total: [-3, -2, -1, 1, 2, 3]
        }[self]
        return [ReadingFrame.from_int(frame) for frame in frames]


def load_test_sequences() -> Dict[Tuple[FrameRejection, bool], List[str]]:
    with open("tests/sequences_test_codons.json") as file:
        sequences = json.load(file)
    return {
        (FrameRejection(entry['rejection']), entry['hasLastStop']): entry['sequences']
        for entry in sequences
    }


TEST_SEQUENCES = load_test_sequences()


@pytest.mark.parametrize("frame_rejection,last_codon_frame, series_gc",
                         itertools.product(
                             tuple(FrameRejection),
                             (False, True),
                             (GeneticCode.Unknown,
                              GeneticCode(1))))
def test_generated_sequences(frame_rejection: FrameRejection,
                             last_codon_frame: bool,
                             series_gc: GeneticCode):
    for seq in TEST_SEQUENCES[(frame_rejection, last_codon_frame)]:
        reading_combinations = detect_reading_combinations(seq, series_gc)
        if frame_rejection == FrameRejection.Total:
            assert ({(gc, frame)
                     for gc, frame in reading_combinations if gc == GeneticCode(1)}
                    ==
                    set())
        else:
            assert (GeneticCode(1), 1) in reading_combinations


@pytest.mark.parametrize("frame_rejection,last_codon_frame,series_frame,series_gc",
                         itertools.product(
                             tuple(FrameRejection),
                             (False, True),
                             (ReadingFrame.Unknown,
                              ReadingFrame.from_int(1),
                              ReadingFrame.from_int(2)),
                             (GeneticCode.Unknown, GeneticCode(1))
                         )
                         )
def test_generated_series(frame_rejection: FrameRejection,
                          last_codon_frame: bool,
                          series_frame: ReadingFrame,
                          series_gc: GeneticCode) -> None:
    rejected_frames = frame_rejection.frames()
    series = pd.Series(TEST_SEQUENCES[frame_rejection, last_codon_frame])
    data = GeneSeries(series)
    data.reading_frame = series_frame
    data.genetic_code = series_gc
    expected_exception: Optional[Exception] = None
    possible_frames = [frame for frame in list(
        ReadingFrame) if frame not in rejected_frames and frame != ReadingFrame.Unknown]
    if not possible_frames:
        expected_exception = NoReadingFrames
    elif series_frame == ReadingFrame.from_int(2):
        expected_exception = BadReadingFrame
    elif (not last_codon_frame and len(possible_frames) > 1 and
          series_frame == ReadingFrame.Unknown):
        expected_exception = NoReadingFrames

    if expected_exception is not None:
        with pytest.raises(expected_exception):
            OpDetectReadingFrame()(data)
    else:
        result = OpDetectReadingFrame()(data)
        assert result.reading_frame == ReadingFrame.from_int(1)


@pytest.fixture
def series_good() -> pd.Series:
    series = pd.Series({
        'seqid1': 'GCAGTACCCTAA',
        'seqid2': 'GCAGTATAA',
    })
    series.name = 'gene1'
    return series


@pytest.fixture
def series_bad() -> pd.Series:
    series = pd.Series({
        'seqid1': 'GCAGTATAATAA',  # has two stop codons
        'seqid2': 'GCAGTATAA',
    })
    series.name = 'gene1'
    return series


@pytest.fixture
def gene_good_standard_1st(series_good) -> GeneSeries:
    data = GeneSeries(series_good)
    data.genetic_code = GeneticCode(1)
    data.reading_frame = ReadingFrame(1)
    return data


@pytest.fixture
def gene_good_standard_unknown(series_good) -> GeneSeries:
    data = GeneSeries(series_good)
    data.genetic_code = GeneticCode(1)
    data.reading_frame = ReadingFrame(0)
    return data


@pytest.fixture
def gene_good_unknown_unknown(series_good) -> GeneSeries:
    data = GeneSeries(series_good)
    data.genetic_code = GeneticCode(0)
    data.reading_frame = ReadingFrame(0)
    return data


@pytest.fixture
def gene_bad_standard_1st(series_bad) -> GeneSeries:
    data = GeneSeries(series_bad)
    data.genetic_code = GeneticCode(1)
    data.reading_frame = ReadingFrame(1)
    return data


def test_good_standard_1st(gene_good_standard_1st):
    result = OpDetectReadingFrame()(gene_good_standard_1st)
    assert result.reading_frame == 1


def test_good_standard_unknown(gene_good_standard_unknown):
    result = OpDetectReadingFrame()(gene_good_standard_unknown)
    assert result.reading_frame == 1


def test_good_unknown_unknown(gene_good_unknown_unknown):
    result = OpDetectReadingFrame()(gene_good_unknown_unknown)
    assert result.reading_frame == 1


def test_bad_standard_1st(gene_bad_standard_1st):
    with pytest.raises(BadReadingFrame):
        OpDetectReadingFrame()(gene_bad_standard_1st)
