
from typing import Callable, Dict
from pathlib import Path

from .model import GeneSeries, GeneDataFrame, GeneStream, GeneIO
from .types import Justification, TextCase
from .utils import ConfigurableCallable, Field, Group
from .file_utils import ZipFile, ZipPath
from .file_types import FileType, FileFormat, get_extension
from .operators import (
    OpCheckValid, OpIndexMerge, OpMakeUniform, OpDropEmpty, OpDropIfAllEmpty,
    OpSequenceCase, OpExtractCharsets, OpTranslateGap, OpPadReadingFrames,
    OpSanitizeGeneNames, OpSanitizeSpeciesNames, OpSpreadsheetCompatibility,
    OpTranslateMissing, OpReverseNegativeReadingFrames)

from . import ali, fasta, phylip
from . import nexus, tabfile
from . import iqtree, partition_finder


class FileWriter(ConfigurableCallable):
    type: FileType = None
    format: FileFormat = None

    case = Field(
        key='case',
        label='Nucleotide Case',
        doc=('Transform all nucleotides to UPPERCASE/lowercase.'),
        type=TextCase,
        list={case: case.description for case in TextCase},
        default=TextCase.Unchanged)

    padding = Field(
        key='padding',
        label='Padding',
        doc=('Used when making sequence length uniform \n'
             'or filling empty sequences.'),
        type=str,
        list=dict(**{k: k for k in '-?*Nn'}, **{'': 'Unpadded'}),
        default='-')

    translate_missing = Field(
        key='translate_missing',
        label='Missing',
        doc=('Character used to indicate a missing nucleotide.\n'
             'Enforced globally.'),
        type=str,
        list={'': 'Unchanged', '?': '?', 'N': 'N', 'n': 'n'},
        default='')

    translate_gap = Field(
        key='translate_gap',
        label='Gap',
        doc=('Character used to indicate a sequence gap.\n'
             'Enforced globally.'),
        type=str,
        list={'': 'Unchanged', '-': '-', '*': '*'},
        default='')

    sanitize = Field(
        key='sanitize',
        label='Sanitize names',
        doc=('Replace unusual characters from species and gene names,\n'
             'either by an underscore or a similar ASCII character.'),
        type=bool,
        default=True)

    def __init__(self, *args, **kwargs):
        """You may append custom filters to self.filters"""
        super().__init__(*args, **kwargs)
        self.filters = list()
        self.filters.append(self.filter)

    def apply_filters(self, stream: GeneStream) -> GeneStream:
        """Apply all filters"""
        for filter in self.filters:
            stream = filter(stream)
        return stream

    def filter(self, stream: GeneStream) -> GeneStream:
        """Stream operations, obeys inheritance"""
        stream = stream.pipe(OpCheckValid())
        if self.translate_missing:
            stream = stream.pipe(OpTranslateMissing(self.translate_missing))
        if self.translate_gap:
            stream = stream.pipe(OpTranslateGap(self.translate_gap))
        if self.case != TextCase.Unchanged:
            stream = stream.pipe(OpSequenceCase(self.case))
        if self.sanitize:
            stream = stream.pipe(OpSanitizeGeneNames())
            stream = stream.pipe(OpSanitizeSpeciesNames())
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        """Write to path after applying filter operations"""
        raise NotImplementedError


file_writers: Dict[FileType, Dict[FileFormat, FileWriter]] = {
    type: dict() for type in FileType}


def file_writer(
    type: FileType, format: FileFormat
) -> Callable[[FileWriter], FileWriter]:
    def decorator(writer: FileWriter) -> FileWriter:
        file_writers[type][format] = writer
        writer.type = type
        writer.format = format
        return writer
    return decorator


class _GeneWriter(FileWriter):
    geneIO: GeneIO = None

    drop_empty = Field(
        key='drop_empty',
        label='Drop Empty Sequences',
        doc=('Exclude all taxa that consist only of gaps and missing data.'),
        type=bool,
        default=True)

    @property
    def params(self) -> Group:
        return Group(key='root', children=[
            self._params_[param] for param in [
                'case',
                'padding',
                'translate_missing',
                'translate_gap',
                'drop_empty',
                'sanitize',
            ]])

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = super().filter(stream)
        if self.drop_empty:
            stream = stream.pipe(OpDropEmpty())
        return stream

    def write(self, stream: GeneStream, path: Path) -> None:
        raise NotImplementedError

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.apply_filters(stream)
        self.write(stream, path)


class _ConcatenatedWriter(_GeneWriter):
    def filter(self, stream: GeneStream) -> GeneStream:
        stream = (
            super().filter(stream)
            .pipe(OpIndexMerge())
            .pipe(OpMakeUniform(self.padding)))
        return stream

    def write(self, stream: GeneStream, path: Path) -> None:
        stream = self.apply_filters(stream)
        joined = GeneDataFrame.from_stream(stream, filler=self.padding)
        data = joined.dataframe.apply(
            lambda row: ''.join(row.values.astype(str)), axis=1)
        gene = GeneSeries(data, missing='', gap='')
        self.geneIO.gene_to_path(gene, path)


class _MultiFileWriter(_GeneWriter):
    @staticmethod
    def create(path: Path) -> Path:
        raise NotImplementedError

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = (
            super().filter(stream)
            .pipe(OpIndexMerge())
            .pipe(OpMakeUniform(self.padding)))
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.apply_filters(stream)
        container = self.create(path)
        for gene in stream:
            name = gene.series.name + get_extension(FileType.File, self.format)
            part = container / name
            self.geneIO.gene_to_path(gene, part)


class _MultiDirWriter(_MultiFileWriter):
    @staticmethod
    def create(path: Path) -> Path:
        path.mkdir(exist_ok=True)
        return path


class _MultiZipWriter(_MultiFileWriter):
    @staticmethod
    def create(path: Path) -> Path:
        archive = ZipFile(path, 'w')
        return ZipPath(archive)


def _register_type_writer(
    ftype: FileType,
    writer: _GeneWriter,
) -> None:

    @file_writer(ftype, FileFormat.Fasta)
    class FastaWriter(writer):
        geneIO = fasta

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.params.translate_missing.value = '?'
            self.params.translate_gap.value = '-'
            self.params.padding.value = '-'

    @file_writer(ftype, FileFormat.Phylip)
    class PhylipWriter(writer):
        geneIO = phylip

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.params.translate_missing.value = '?'
            self.params.translate_gap.value = '-'
            self.params.padding.value = '-'

    @file_writer(ftype, FileFormat.Ali)
    class AliWriter(writer):
        geneIO = ali

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.params.translate_missing.value = '?'
            self.params.translate_gap.value = '*'
            self.params.padding.value = '*'


for ftype, writer in {
    FileType.File: _ConcatenatedWriter,
    FileType.Directory: _MultiDirWriter,
    FileType.ZipArchive: _MultiZipWriter,
}.items():
    _register_type_writer(ftype, writer)


class _ContainerWriter(FileWriter):

    adjust_frames = Field(
        key='adjust_frames',
        label='Adjust reading frames',
        doc=('Genes with negative reading frames are reverse-complemented.\n'
             'Pad coding sequences to start with first codon position.'),
        type=bool,
        default=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.op_charsets = None
        self.params.translate_missing.value = '?'
        self.params.translate_gap.value = '-'
        self.params.padding.value = '-'

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = (
            super().filter(stream)
            .pipe(OpIndexMerge())
            .pipe(OpMakeUniform(self.padding)))
        if self.adjust_frames:
            padding = self.translate_missing or self.padding
            stream = stream.pipe(OpReverseNegativeReadingFrames())
            stream = stream.pipe(OpPadReadingFrames(padding))
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        self.op_charsets = OpExtractCharsets()
        stream = self.apply_filters(stream)
        container = self.create(path)
        stream = stream.pipe(self.op_charsets)
        joined = GeneDataFrame.from_stream(stream, filler=self.padding)
        data = joined.dataframe.apply(
            lambda row: ''.join(row.values.astype(str)), axis=1)
        gene = GeneSeries(data, missing='', gap='')
        phylip.gene_to_path(gene, container / self.alignment)
        self.write_config(container)

    @staticmethod
    def create(path: Path) -> Path:
        raise NotImplementedError

    def write_config(self, container: Path) -> None:
        raise NotImplementedError


class _IQTreeWriter(_ContainerWriter):
    alignment = Field(
        key='alignment',
        label='Alignment File',
        doc=('The name for the PHYLIP alignment file.\n'
             'The NEXUS configuration file will have the same base name.'),
        type=str,
        default='alignment.phy')

    include_full_markers_with_codons = Field(
        key='include_full_markers_with_codons',
        label='Include full markers with codons',
        doc=('Include the full marker definition for markers \n'
             'that have codon subsets defined via reading frames.'),
        type=bool,
        default=False)

    @property
    def params(self) -> Group:
        return Group(key='root', children=[
            self._params_[param] for param in [
                'alignment',
                'case',
                'padding',
                'translate_missing',
                'translate_gap',
                'include_full_markers_with_codons',
                'adjust_frames',
                'sanitize',
            ]])

    def write_config(self, container: Path) -> None:
        cfg_name = self.alignment.split('.')[0] + '.nex'
        iqtree.write_nexus(
            container / cfg_name,
            self.op_charsets.charsets,
            self.include_full_markers_with_codons)


class _PartitionFinderWriter(_ContainerWriter):
    alignment = Field(
        key='alignment',
        label='Alignment File',
        doc=('The name for the PHYLIP alignment file.'),
        type=str,
        default='alignment.phy')

    cfg_file = Field(
        key='cfg_file',
        label='Configuration File',
        doc=('The name for the PartitionFinder configuration file.'),
        type=str,
        default='partition_finder.cfg')

    @property
    def params(self) -> Group:
        return Group(key='root', children=[
            self._params_[param] for param in [
                'alignment',
                'cfg_file',
                'case',
                'padding',
                'translate_missing',
                'translate_gap',
                'adjust_frames',
                'sanitize',
            ]])

    def write_config(self, container: Path) -> None:
        partition_finder.write_cfg(
            container / self.cfg_file,
            self.op_charsets.charsets,
            self.alignment)


def _register_container_writer(
    format: FileFormat,
    writer: _ContainerWriter,
) -> None:

    @file_writer(FileType.Directory, format)
    class _ContainerDirWriter(writer):
        @staticmethod
        def create(path: Path) -> Path:
            path.mkdir(exist_ok=True)
            return path

    @file_writer(FileType.ZipArchive, format)
    class _ContainerZipWriter(writer):
        @staticmethod
        def create(path: Path) -> Path:
            archive = ZipFile(path, 'w')
            return ZipPath(archive)


for format, writer in {
    FileFormat.PartitionFinder: _PartitionFinderWriter,
    FileFormat.IQTree: _IQTreeWriter,
}.items():
    _register_container_writer(format, writer)


@file_writer(FileType.File, FileFormat.Nexus)
class NexusWriter(FileWriter):
    justification = Field(
        key='justification',
        label='Justification',
        doc=('Text justification for the species names.'),
        type=Justification,
        list={just: just.description for just in Justification},
        default=Justification.Left)

    separator = Field(
        key='separator',
        label='Separator',
        doc=('Used to separate species names and sequence data.'),
        type=str,
        list={' ': 'Space', '\t': 'Tab'},
        default=' ')

    adjust_frames = Field(
        key='adjust_frames',
        label='Adjust reading frames',
        doc=('Genes with negative reading frames are reverse-complemented.\n'
             'Pad coding sequences to start with first codon position.'),
        type=bool,
        default=True)

    include_full_markers_with_codons = Field(
        key='include_full_markers_with_codons',
        label='Include full markers with codons',
        doc=('Include the full marker definition for markers \n'
             'that have codon subsets defined via reading frames.'),
        type=bool,
        default=False)

    @property
    def params(self) -> Group:
        return Group(key='root', children=[
            self._params_[param] for param in [
                'case',
                'padding',
                'translate_missing',
                'translate_gap',
                'justification',
                'separator',
                'include_full_markers_with_codons',
                'adjust_frames',
                'sanitize',
            ]])

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.params.translate_missing.value = '?'
        self.params.translate_gap.value = '-'

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = super().filter(stream)
        stream = stream.pipe(OpIndexMerge())
        stream = (
            GeneDataFrame
            .from_stream(stream, filler=self.padding)
            .to_stream()
            .pipe(OpDropIfAllEmpty()))
        if self.adjust_frames:
            padding = self.translate_missing or self.padding
            stream = stream.pipe(OpReverseNegativeReadingFrames())
            stream = stream.pipe(OpPadReadingFrames(padding))
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.apply_filters(stream)
        nexus.stream_to_path(
            stream, path,
            self.justification,
            self.separator,
            self.include_full_markers_with_codons)


@file_writer(FileType.File, FileFormat.Tab)
class TabWriter(FileWriter):
    sequence_prefix = Field(
        key='sequence_prefix',
        label='Sequence Prefix',
        doc=('Added to the beginning of all sequence names.'),
        type=str,
        default='sequence_')

    spreadsheet = Field(
        key='spreadsheet',
        label='Spreadsheet Compatibility',
        doc=('Coding sequences starting with "-" will start with "N" instead'),
        type=bool,
        default=False)

    @property
    def params(self) -> Group:
        return Group(key='root', children=[
            self._params_[param] for param in [
                'sequence_prefix',
                'case',
                'translate_missing',
                'translate_gap',
                'sanitize',
                'spreadsheet',
            ]])

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.params.translate_missing.value = '?'
        self.params.translate_gap.value = '-'

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = super().filter(stream)
        if self.spreadsheet:
            stream = stream.pipe(OpSpreadsheetCompatibility())
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.apply_filters(stream)
        tabfile.stream_to_path(stream, path, self.sequence_prefix)


class WriterNotFound(Exception):
    def __init__(self, type: FileType, format: FileFormat):
        self.type = type
        self.format = format
        super().__init__(f'No writer for {str(type)} and {str(format)}')


def get_writer(type: FileType, format: FileFormat):
    if format not in file_writers[type]:
        raise WriterNotFound(type, format)
    return file_writers[type][format]()


def write_to_path(
    series: GeneStream,
    path: Path,
    type: FileType,
    format: FileFormat,
) -> None:
    writer = get_writer(type, format)
    return writer(series, path)
