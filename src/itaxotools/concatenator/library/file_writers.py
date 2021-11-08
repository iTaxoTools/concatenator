
from typing import Callable, Dict
from pathlib import Path

from .model import GeneSeries, GeneDataFrame, GeneStream, GeneIO
from .types import Justification, TextCase
from .utils import ConfigurableCallable, Field, Group
from .file_utils import ZipFile, ZipPath
from .file_types import FileType, FileFormat, get_extension
from .operators import (
    OpCheckValid, OpIndexMerge, OpPadRight, OpDropEmpty,
    OpSequenceCase, OpExtractCharsets, OpApplyToSeries)

from . import ali, fasta, phylip
from . import nexus, tabfile
from . import iqtree, partition_finder


class FileWriter(ConfigurableCallable):
    type: FileType = None
    format: FileFormat = None

    case = Field(
        key='case',
        label='Sequence Case',
        doc=('...'),
        type=TextCase,
        list={case: case.description for case in TextCase},
        default=TextCase.Unchanged)

    def filter(self, stream: GeneStream) -> GeneStream:
        """Stream operations, obeys inheritance"""
        return (
            stream
            .pipe(OpCheckValid())
            .pipe(OpSequenceCase(self.case))
            )

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

    def write(self, stream: GeneStream, path: Path) -> None:
        raise NotImplementedError

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        self.write(stream, path)


class _ConcatenatedWriter(_GeneWriter):
    padding = Field('padding', value='')

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = (
            super().filter(stream)
            .pipe(OpIndexMerge())
            .pipe(OpPadRight(self.padding)))
        return stream

    def write(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        joined = GeneDataFrame.from_stream(stream)
        data = joined.dataframe.apply(
            lambda row: ''.join(row.values.astype(str)), axis=1)
        gene = GeneSeries(data, missing='', gap='')
        self.geneIO.gene_to_path(gene, path)


class _MultiFileWriter(_GeneWriter):
    padding = Field('padding', value='')

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

    @file_writer(ftype, FileFormat.Phylip)
    class PhylipWriter(writer):
        geneIO = phylip

    @file_writer(ftype, FileFormat.Ali)
    class AliWriter(writer):
        geneIO = ali


for ftype, writer in {
    FileType.File: _ConcatenatedWriter,
    FileType.Directory: _MultiDirWriter,
    FileType.ZipArchive: _MultiZipWriter,
}.items():
    _register_type_writer(ftype, writer)


class _ContainerWriter(FileWriter):
    padding = Field('padding', value='-')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.op_charsets = OpExtractCharsets()

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
        stream = stream.pipe(self.op_charsets)
        joined = GeneDataFrame.from_stream(stream)
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
    alignment = Field('alignment', value='alignment.phy')

    def write_config(self, container: Path) -> None:
        cfg_name = self.alignment.split('.')[0] + '.nex'
        iqtree.write_nexus(
            container / cfg_name,
            self.op_charsets.charsets)


class _PartitionFinderWriter(_ContainerWriter):
    alignment = Field('alignment', value='alignment.phy')
    cfg_file = Field('cfg_file', value='partition_finder.cfg')

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
    padding = Field(
        key='padding',
        label='Padding',
        doc=('...'),
        type=str,
        list=list('-?*Nn'),
        default='-')
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
                'justification',
                'separator',
                'padding',
            ]])

    def filter(self, stream: GeneStream) -> GeneStream:
        stream = super().filter(stream)
        stream = (
            GeneDataFrame.from_stream(stream).to_stream()
            .pipe(OpIndexMerge())
            .pipe(OpApplyToSeries(lambda x: x.fillna('')))
            .pipe(OpPadRight(self.padding)))
        return stream

    def call(self, stream: GeneStream, path: Path) -> None:
        stream = self.filter(stream)
        nexus.stream_to_path(stream, path, self.justification, self.separator)


@file_writer(FileType.File, FileFormat.Tab)
class TabWriter(FileWriter):
    sequence_prefix = Field('sequence_prefix', value='sequence_')

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
