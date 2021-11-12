from typing import List, Tuple, Optional
from enum import IntEnum, Enum
import random
import itertools

import pandas as pd
import pytest

from itaxotools.concatenator.library.model import GeneSeries
from itaxotools.concatenator.library.codons import (
    GeneticCode, ReadingFrame, NoReadingFrames,
    BadReadingFrame, AmbiguousReadingFrame, last_codon_slice,
    detect_reading_combinations)
from itaxotools.concatenator.library.operators import OpDetectReadingFrame


def no_stop_sequence(len: int, gc: GeneticCode) -> str:
    seq = "".join(random.choices("ACGT", k=len))
    for stop in gc.stops:  # type: ignore
        seq = seq.replace(stop, "")
    return seq


def generate_sequence(seq_len: int, gc: GeneticCode,
                      frames_with_two_stops: List[ReadingFrame],
                      frame_with_last_codon: Optional[ReadingFrame]) -> str:
    """
    Generates a random sequence with approximate length `len`.

    All frames in `frames_with_two_stops` will contains at two stops.

    The frame `frame_with_last_codon` will contain a stop in last codon.
    """
    assert ReadingFrame.Unknown not in frames_with_two_stops
    assert frame_with_last_codon != ReadingFrame.Unknown
    segment_codon_len = (seq_len // len(frames_with_two_stops)) // 3
    positive_frames = [frame for frame in frames_with_two_stops if frame > 0]
    negative_frames = [-frame for frame in frames_with_two_stops if frame < 0]
    positive_part = ""
    for frame in positive_frames:
        positive_part += no_stop_sequence(segment_codon_len * 3 + frame - 1, gc)
        positive_part += random.choice(gc.stops)
    negative_part = ""
    for frame in negative_frames:
        negative_part += no_stop_sequence(segment_codon_len * 3 + frame - 1, gc)
        negative_part += random.choice(gc.stops)
    seq = positive_part + negative_part[::-1]
    if frame_with_last_codon:
        the_slice = last_codon_slice(seq, frame_with_last_codon)
        seq = seq[:the_slice.start] + random.choice(gc.stops) + seq[the_slice.stop:]
    return seq


def generated_series(
    frames_with_two_stops: List[ReadingFrame],
    frame_with_last_codon: Optional[ReadingFrame]
) -> pd.Series:
    series = pd.Series((generate_sequence(300,
                                          GeneticCode(1),
                                          frames_with_two_stops,
                                          frame_with_last_codon)
                        for _ in range(8)))
    series.name = 'gene1'
    return series


class FrameRejection(Enum):
    Ambiguous = 0
    Unambiguos = 1
    Total = 2

    def frames(self):
        frames = {
            FrameRejection.Ambiguous: [-3, -1, 2, 3],
            FrameRejection.Unambiguos: [-3, -2, -1, 2, 3],
            FrameRejection.Total: [-3, -2, -1, 1, 2, 3]
        }[self]
        return [ReadingFrame.from_int(frame) for frame in frames]


@pytest.mark.skip(reason="Fail because of incorrect sequence generation")
@pytest.mark.parametrize("frame_rejection,last_codon_frame, series_gc",
                         itertools.product(
                             tuple(FrameRejection),
                             (None, ReadingFrame.from_int(1)),
                             (GeneticCode.Unknown,
                              GeneticCode(1))))
def test_generated_sequence(frame_rejection, last_codon_frame, series_gc):
    seq = generate_sequence(300, GeneticCode(
        1), frame_rejection.frames(), last_codon_frame)
    if FrameRejection == FrameRejection.Total:
        assert detect_reading_combinations(seq, series_gc) == set()
    else:
        assert (GeneticCode(1), 1) in detect_reading_combinations(seq, series_gc)


@pytest.mark.skip(reason="Fail because of incorrect sequence generation")
@pytest.mark.parametrize("frame_rejection,last_codon_frame,series_frame,series_gc",
                         itertools.product(
                             tuple(FrameRejection),
                             (None, ReadingFrame.from_int(1)),
                             (ReadingFrame.Unknown,
                              ReadingFrame.from_int(1),
                              ReadingFrame.from_int(2)),
                             (GeneticCode.Unknown, GeneticCode(1))
                         )
                         )
def test_generated_series(frame_rejection: FrameRejection,
                          last_codon_frame: Optional[ReadingFrame],
                          series_frame: ReadingFrame,
                          series_gc: GeneticCode) -> None:
    rejected_frames = frame_rejection.frames()
    series = generated_series(rejected_frames, last_codon_frame)
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


def test_model_reading_frame():
    assert not bool(ReadingFrame(0))
    assert ReadingFrame(1) == 1
    assert ReadingFrame(1).label == '+1'
    assert ReadingFrame(-2).label == '-2'


def test_model_genetic_code():
    assert not bool(GeneticCode(0))
    assert not GeneticCode(0).stops
    assert GeneticCode(1) == 1
    assert GeneticCode(1).stops
    assert GeneticCode(1).text == 'Standard'
    assert GeneticCode(0).text == 'Unknown'
