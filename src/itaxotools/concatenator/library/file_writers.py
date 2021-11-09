
from typing import Callable, Dict
from pathlib import Path

from .model import GeneSeries, GeneDataFrame, GeneStream, GeneIO
from .types import Justification, TextCase
from .utils import ConfigurableCallable, Field, Group
from .file_utils import ZipFile, ZipPath
from .file_types import FileType, FileFormat, get_extension
from .operators import (
    OpCheckValid, OpIndexMerge, OpPadRight, OpDropEmpty,
    OpSequenceCase, OpExtractCharsets, OpApplyToSeries,
    OpSanitizeGeneNames, OpSanitizeSpeciesNames, OpSpreadsheetCompatibility,
    OpTranslateMissing, OpTranslateGap)

from . import ali, fasta, phylip
from . import nexus, tabfile
from . import iqtree, partition_finder


class FileWriter(ConfigurableCallable):
    type: FileType = None
    format: FileFormat = None

    case = Field(
        key='case',
        label='Nucleotide Case',
        doc=('...'),
        type=TextCase,
        list={case: case.description for case in TextCase},
        default=TextCase.Unchanged)

    padding = Field(
        key='padding',
        label='Padding',
        doc=('...'),
        type=str,
        list=dict(**{k: k for k in '-?*Nn'}, **{'': 'Unpadded'}),
        default='-')

    translate_missing = Field(
        key='translate_missing',
        label='Translate Missing',
        doc=('...'),
        type=str,
        list={'': 'Unchanged', '?': '?', 'N': 'N', 'n': 'n'},
        default='')

    translate_gap = Field(
        key='translate_gap',
        label='Translate Gap',
        doc=('...'),
        type=str,
        list={'': 'Unchanged', '-': '-', '*': '*'},
        default='')

    sanitize_genes = Field(
        key='sanitize_genes',
        label='Sanitize Gene Names',
        doc=('...'),
        type=bool,
        default=True)

    sanitize_species = Field(
        key='sanitize_species',
        label='Sanitize Species Names',
        doc=('...'),
        type=bool,
        default=True)

    def filter(self, stream: GeneStream) -> GeneStream:
        """Stream operations, obeys inheritance"""
        stream = stream.pipe(OpCheckValid())
        if self.translate_missing:
            stream = stream.pipe(OpTranslateMissing(self.translate_missing))
        if self.translate_gap:
            stream = stream.pipe(OpTranslateGap(self.translate_gap))
        stream = stream.pipe(OpDropEmpty())
        stream = stream.pipe(OpSequenceCase(self.case))
        if self.sanitize_genes:
            stream = stream.pipe(OpSanitizeGeneNames())
        if self.sanitize_species:
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

    @property
    def params(self) -> Group:
        return Group(key='root', children=[
            self._params_[param] for param in [
                'case',
                'padding',
                'translate_missing',
                'translate_gap',
                'sanitize_genes',
                'sanitize_species',
            ]])

    def write(self, stream: GeneStream, path: Path) -> None:
        raise NotImplementedError

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        self.write(stream, path)


class _ConcatenatedWriter(_GeneWriter):
    def filter(self, stream: GeneStream) -> GeneStream:
        stream = (
            super().filter(stream)
            .pipe(OpIndexMerge())
            .pipe(OpPadRight(self.padding)))
        return stream

    def write(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
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
            .pipe(OpDropEmpty())
            .pipe(OpIndexMerge())
            .pipe(OpPadRight(self.padding)))
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
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
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.op_charsets = OpExtractCharsets()
        self.params.translate_missing.value = '?'
        self.params.translate_gap.value = '-'
        self.params.padding.value = '-'

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = (
            super().filter(stream)
            .pipe(OpIndexMerge())
            .pipe(OpPadRight(self.padding)))
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
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
        doc=('...'),
        type=str,
        default='alignment.phy')

    @property
    def params(self) -> Group:
        return Group(key='root', children=[
            self._params_[param] for param in [
                'alignment',
                'case',
                'padding',
                'translate_missing',
                'translate_gap',
                'sanitize_genes',
                'sanitize_species',
            ]])

    def write_config(self, container: Path) -> None:
        cfg_name = self.alignment.split('.')[0] + '.nex'
        iqtree.write_nexus(
            container / cfg_name,
            self.op_charsets.charsets)


class _PartitionFinderWriter(_ContainerWriter):
    alignment = Field(
        key='alignment',
        label='Alignment File',
        doc=('...'),
        type=str,
        default='alignment.phy')

    cfg_file = Field(
        key='cfg_file',
        label='Configuration File',
        doc=('...'),
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
                'sanitize_genes',
                'sanitize_species',
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
        doc=('...'),
        type=Justification,
        list={just: just.description for just in Justification},
        default=Justification.Left)

    separator = Field(
        key='separator',
        label='Separator',
        doc=('...'),
        type=str,
        list={' ': 'Space', '\t': 'Tab'},
        default=' ')

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
                'sanitize_genes',
                'sanitize_species',
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
            .to_stream())
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        nexus.stream_to_path(stream, path, self.justification, self.separator)


@file_writer(FileType.File, FileFormat.Tab)
class TabWriter(FileWriter):
    sequence_prefix = Field(
        key='sequence_prefix',
        label='Sequence Prefix',
        doc=('...'),
        type=str,
        default='sequence_')

    spreadsheet = Field(
        key='spreadsheet',
        label='Spreadsheet Compatibility',
        doc=('...'),
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
                'sanitize_genes',
                'sanitize_species',
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
        stream = self.filter(stream)
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
