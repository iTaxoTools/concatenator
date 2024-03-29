from typing import Callable, Dict, List, Iterator, Iterable, Optional, Set
from collections import defaultdict

import pandas as pd
import re

from .model import Operator, GeneSeries, GeneStream, GeneDataFrame
from .types import TextCase, Charset
from .utils import (
    Translation,
    Field,
    OrderedSet,
    removeprefix,
    sanitize,
    reverse_complement,
    has_uniform_length,
)
from .codons import final_column_reading_frame, ReadingFrame
from .general_info import (
    GeneralInfo,
    InfoColumns,
    FileGeneralInfo,
    GeneInfoColumns,
    GeneInfo,
)


class InvalidGeneSeries(Exception):
    def __init__(self, gene: GeneSeries, error: str):
        self.gene = gene
        super().__init__(error)


class OpPass(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        return gene


class OpBlock(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        return None


class OpCheckDuplicates(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._memory = set()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name in self._memory:
            raise InvalidGeneSeries(gene, f"Duplicate gene: {repr(gene.name)}")
        self._memory.add(gene.name)
        return gene


class OpCheckValid(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.op_check_duplicates = OpCheckDuplicates()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if not isinstance(gene, GeneSeries):
            raise InvalidGeneSeries(gene, "Not a GeneSeries!")
        if not isinstance(gene.series, pd.Series):
            raise InvalidGeneSeries(gene, "Gene series is not pandas.Series!")
        if gene.series.index.duplicated().any():
            raise InvalidGeneSeries(gene, "Duplicate indices")
        if not gene.name:
            raise InvalidGeneSeries(gene, 'Missing gene data: "name"')
        if not gene.missing:
            raise InvalidGeneSeries(gene, 'Missing gene data: "missing"')
        if not gene.gap:
            raise InvalidGeneSeries(gene, 'Missing gene data: "gap"')
        return self.op_check_duplicates(gene)


class OpTagSet(Operator):
    tag: str = Field("tag", value="")
    value: object = Field("value", value=None)

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        assert self.tag.isidentifier()
        gene = gene.copy()
        gene.tags[self.tag] = self.value
        return gene


class OpTagDelete(Operator):
    tag: str = Field("tag", value="")

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        assert self.tag.isidentifier()
        gene = gene.copy()
        del gene.tags[self.tag]
        return gene


class OpIndexMerge(Operator):
    index: str = Field("index", value="seqid")
    glue: str = Field("glue", value="_")

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        indices = gene.series.index.to_frame().fillna("", inplace=False)
        gene.series.index = pd.Index(
            indices.apply(self.glue.join, axis=1), name=self.index
        )
        return gene


class OpIndexFilter(Operator):
    index: str = "seqid"

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series.index = gene.series.index.to_frame()[self.index]
        return gene


class OpDropEmpty(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series = gene.series.dropna(inplace=False)
        empty = "".join(set(gene.missing + gene.gap))
        if empty:
            gene.series = gene.series[
                ~gene.series.str.fullmatch(f"[{re.escape(empty)}]+")
            ]
        gene.series = gene.series[gene.series.str.len() > 0]
        if not len(gene.series):
            return None
        return gene


class OpDropIfAllEmpty(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        temp = gene.copy()
        temp = OpDropEmpty()(temp)
        if temp is None:
            return None
        return gene


class OpMakeUniform(Operator):
    padding: str = Field("padding", value="")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert len(self.padding) in [0, 1]

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if not self.padding:
            return gene
        gene = gene.copy()
        gene.series = gene.series.fillna("")
        max_length = gene.series.str.len().max()
        if gene.reading_frame < 0:
            gene.series = gene.series.str.rjust(max_length, self.padding)
        else:
            gene.series = gene.series.str.ljust(max_length, self.padding)
        return gene


class OpTranslateSequences(Operator):
    translation: Translation = Field("translation", value={})

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series = gene.series.str.translate(self.translation)
        return gene


class OpTranslateMissing(Operator):
    missing: str = Field("missing", value="?")

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if not self.missing or not gene.missing:
            return gene
        assert len(self.missing) == 1
        mapping = {char: self.missing for char in gene.missing}
        translation = str.maketrans(mapping)
        gene = OpTranslateSequences(translation)(gene)
        gene.missing = self.missing
        return gene


class OpTranslateGap(Operator):
    gap: str = Field("gap", value="-")

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if not self.gap or not gene.gap:
            return gene
        assert len(self.gap) == 1
        mapping = {char: self.gap for char in gene.gap}
        translation = str.maketrans(mapping)
        gene = OpTranslateSequences(translation)(gene)
        gene.gap = self.gap
        return gene


class OpFilterGenes(Operator):
    genes: Set = Field("genes", value=set())

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name not in self.genes:
            return None
        return gene


class OpStencilGenes(Operator):
    operator: Operator = Field("operator", value=OpPass())
    genes: Iterable = Field("genes", value=list())

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name in self.genes:
            return self.operator(gene)
        return gene


class OpTranslateGenes(Operator):
    translation: Dict = Field("translation", value=dict())

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name in self.translation:
            if self.translation[gene.name] is None:
                return None
            gene = gene.copy()
            gene.name = self.translation[gene.name]
        return gene


class OpChainGenes(Operator):
    allow_duplicates: bool = Field("allow_duplicates", value=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._memory = set()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if self.allow_duplicates:
            return gene
        if gene.name in self._memory:
            return None
        self._memory.add(gene.name)
        return gene


class OpDetectReadingFrame(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        final_reading_frame = final_column_reading_frame(
            gene.series, gene.genetic_code, gene.reading_frame
        )
        gene.reading_frame = final_reading_frame
        return gene


class OpSanitizeGeneNames(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.name = sanitize(gene.name)
        return gene


class OpSanitizeSpeciesNames(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        indices = gene.series.index.names
        data = gene.series.reset_index()
        try:
            data[indices] = data[indices].map(sanitize)
        except Exception as e:
            print(data[indices])
            raise e
        gene.series = data.set_index(indices).iloc[:, 0]
        return gene


class OpSequenceCase(Operator):
    case: TextCase = Field("case", value=TextCase.Unchanged)

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if self.case.method is None:
            return gene
        gene = gene.copy()
        gene.series = gene.series.apply(self.case.method)
        return gene


class OpSpreadsheetCompatibility(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series = gene.series.str.replace("^-", "N", regex=True)
        return gene


class OpReverseComplement(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        gene = gene.copy()
        gene.series = gene.series.map(reverse_complement)
        return gene


class OpReverseNegativeReadingFrames(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.reading_frame >= 0:
            return gene
        gene = gene.copy()
        gene.reading_frame = ReadingFrame(-gene.reading_frame)
        return OpReverseComplement()(gene)


class OpPadReadingFrames(Operator):
    padding: str = Field("padding", value="")
    lengths = {2: 2, 3: 1}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert len(self.padding) in [0, 1]

    @staticmethod
    def pad_left(seq: str, pad: str) -> str:
        return f"{pad}{seq}"

    @staticmethod
    def pad_right(seq: str, pad: str) -> str:
        return f"{seq}{pad}"

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.reading_frame in [0, 1, -1]:
            return gene
        gene = gene.copy()
        func = self.pad_left if gene.reading_frame > 0 else self.pad_right
        reading_frame = abs(gene.reading_frame)
        if self.padding:
            padding = self.padding
        else:
            padding = "?" if "?" in gene.missing else gene.missing[0]
        padding = padding * self.lengths[reading_frame]
        gene.series = gene.series.apply(lambda seq: func(seq, padding))
        sign = 1 - 2 * int(gene.reading_frame < 0)
        gene.reading_frame = ReadingFrame(sign)
        return gene


class OpExtractCharsets(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.charsets = list()
        self.cursor = 1

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        assert has_uniform_length(gene.series)
        length = len(gene.series.iat[0])
        charset = Charset(
            gene.name, self.cursor, length, gene.reading_frame, gene.codon_names
        )
        print('####', charset)
        self.charsets.append(charset)
        self.cursor += length
        return gene


class OpUpdateMetadata(Operator):
    metas: Dict[str, Dict] = Field("metas", value={})

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if gene.name not in self.metas:
            return gene
        gene = gene.copy()
        for k, v in self.metas[gene.name].items():
            setattr(gene, k, v)
        return gene


class OpApplyToGene(Operator):
    func: Callable[[GeneSeries], GeneSeries] = Field("func", value=None)

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if self.func is None:
            return gene
        return self.func(gene)


class OpApplyToSeries(Operator):
    func: Callable[[pd.Series], pd.Series] = Field("func", value=None)

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        if self.func is None:
            return gene
        gene = gene.copy()
        gene.series = self.func(gene.series)
        return gene


class OpGeneralInfo(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.table = GeneralInfo.empty()

    def call(self, original: GeneSeries) -> Optional[GeneSeries]:
        gene = OpIndexMerge(index="taxon")(original)
        assert gene.series is not None
        dataframe = pd.DataFrame(gene.series.str.len()).rename(
            columns=(lambda _: InfoColumns.NucleotideCount)
        )
        missing_regex = re.compile(
            "|".join(re.escape(c) for c in gene.missing + gene.gap)
        )
        dataframe[InfoColumns.MissingCount] = gene.series.str.count(missing_regex)
        # NucleotideCount only counts non-missing nucleotides
        dataframe[InfoColumns.NucleotideCount] -= dataframe[InfoColumns.MissingCount]
        dataframe[InfoColumns.SeqCount] = 1
        dataframe[InfoColumns.SeqLenMax] = dataframe[InfoColumns.NucleotideCount]
        dataframe[InfoColumns.SeqLenMin] = dataframe[InfoColumns.NucleotideCount]
        dataframe[InfoColumns.GCCount] = gene.series.str.count(re.compile(r"G|g|C|c"))
        dataframe = dataframe.reset_index().rename(columns={"taxon": InfoColumns.Taxon})
        dataframe[InfoColumns.Gene] = gene.name
        dataframe.set_index(
            [InfoColumns.Gene, InfoColumns.Taxon], inplace=True, verify_integrity=True
        )
        # drop rows where the sequence is completely missing
        dataframe = dataframe.loc[dataframe[InfoColumns.NucleotideCount] > 0].copy()
        self.table += GeneralInfo(dataframe)
        return original


class OpGeneralInfoPerFile(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ops = defaultdict(lambda: OpGeneralInfo())
        self.sources = dict()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        assert gene.stream is not None
        if gene.stream.id not in self.sources:
            self.sources[gene.stream.id] = gene.stream.source
        op = self.ops[gene.stream.id]
        return op(gene)

    def get_info(self) -> pd.DataFrame:
        return GeneralInfo.by_input_file(self._file_infos())

    def _file_infos(self) -> Iterator[FileGeneralInfo]:
        for id in self.ops.keys():
            filename = self.sources[id].path.name
            file_format = self.sources[id].format
            table = self.ops[id].table
            yield FileGeneralInfo(filename, file_format, table)


class OpGeneralInfoPerGene(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.gene_info = pd.DataFrame(columns=[col for col in GeneInfoColumns])
        self.gene_info.index.name = InfoColumns.Gene

    def _append_gene_tags(self, gene: str, mafft: bool, length: bool, codon: bool):
        self.gene_info.loc[gene] = [mafft, length, codon]

    def call(self, original: GeneSeries) -> Optional[GeneSeries]:
        gene = OpIndexMerge(index="taxon")(original)
        self._append_gene_tags(
            gene.name,
            gene.tags.get("MafftRealigned", False),
            gene.tags.get("PaddedLength", False),
            gene.tags.get("PaddedCodonPosition", False),
        )
        return original

    def get_info(self, general_info: GeneralInfo) -> pd.DataFrame:
        gene_info = GeneInfo(self.gene_info)
        return general_info.by_gene(gene_info)


class OpGeneralInfoTagMafftRealigned(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        return OpTagSet('MafftRealigned', True)(gene)


class OpGeneralInfoTagPaddedLength(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        padded = not has_uniform_length(gene.series)
        return OpTagSet('PaddedLength', padded)(gene)


class OpGeneralInfoTagPaddedCodonPosition(Operator):
    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        padded = int(gene.reading_frame) not in [0, 1, -1]
        return OpTagSet('PaddedCodonPosition', padded)(gene)


# Pending removal, functionality to be merged into GeneDataFrame?
def _join_any(stream: GeneStream) -> GeneDataFrame:
    sentinel = "\u0000"
    all_keys = OrderedSet()

    def guarded(names: Iterator[str]) -> List[str]:
        return [name for name in names if name.startswith(sentinel)]

    def guard(names: Iterator[str]) -> List[str]:
        return [sentinel + name for name in names]

    def unguard(names: Iterator[str]) -> List[str]:
        return [removeprefix(name, sentinel) for name in names]

    def fold_keys(stream: GeneStream) -> Iterator[GeneSeries]:
        for gene in stream:
            gene = gene.copy()
            keys = guard(gene.series.index.names)
            gene.series.index = pd.MultiIndex.from_frame(
                gene.series.index.to_frame(), names=keys
            )
            all_keys.update(keys)
            gene.series = gene.series.reset_index(keys)
            yield gene

    genes = fold_keys(stream)
    gdf = GeneDataFrame.from_gene(next(genes))
    all = gdf.dataframe
    for gene in genes:
        merge_keys = all_keys & guarded(gene.series.columns)
        all = pd.merge(all, gene.series, how="outer", on=list(merge_keys))
        # We assume all gene series are compatible
        print(gene, gene.missing)
        gdf.missing = gene.missing
        gdf.gap = gene.gap
    print("end", gdf.missing)
    all.set_index(list(all_keys), inplace=True)
    all.index.names = unguard(all.index.names)
    gdf.dataframe = all
    return gdf
