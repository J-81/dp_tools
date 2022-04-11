""" Functions for loading into data model objects 

Loaders SHOULD:
  - run validation on the dataSystem object
  - return a dataSystem object
"""
import os
from pathlib import Path
import logging
from typing import List, Protocol, Union
import pkg_resources

from schema import Schema, Or, SchemaMissingKeyError
import pandas as pd
from dp_tools.components.components import (
    DatasetGeneCounts,
    DifferentialGeneExpression,
    GeneCounts,
    GenomeAlignments,
    NormalizedGeneCounts,
    RSeQCAnalysis,
    TrimReadsComponent,
)
import yaml

log = logging.getLogger(__name__)

from dp_tools.core.entity_model import (
    BaseComponent,
    BaseSample,
    BaseDataSystem,
    BaseDataset,
    DataDir,
    DataFile,
    GLDSDataSystem,
    TemplateDataSystem,
)

from dp_tools.bulkRNASeq.entity import BulkRNASeqDataset, BulkRNASeqSample
from dp_tools.bulkRNASeq.locaters import find_data_asset_path
from dp_tools.components import RawReadsComponent, BulkRNASeqMetadataComponent

def load_config(config: Union[str, Path]) -> dict:
    if isinstance(config, str):
        conf_data_assets = yaml.safe_load(pkg_resources.resource_string(__name__, f"../config/bulkRNASeq_v{config}.yaml"))["data assets"]
    elif isinstance(config, Path):
        conf_data_assets = yaml.safe_load(config.open())["data assets"]

    # validate with schema
    config_schema = Schema({
        "processed location": [str],
        "resource categories": {
            "subcategory": Or(None, str),
            "subdirectory": Or(None, str),
            "publish to repo": bool
        }
    })
    for key, conf_data_asset in conf_data_assets.items():
        try:
            config_schema.validate(conf_data_asset)
        except SchemaMissingKeyError as e:
            raise ValueError(f"Data asset config: '{key}' failed validation") from e
    
    return conf_data_assets

# TODO: Attach/associate data assets config with dataset/datasystem
def load_BulkRNASeq_STAGE_00(
    root_path: Path, config: Union[str, Path] = "Latest", dataSystem_name: str = None, stack: bool = False
):
    """ config is either an string version referring to a packaged config file or a path to a local config file. """
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # if not specified, assume root_path name is the dataSystem name
    if not dataSystem_name:
        dataSystem_name = root_path.name
    # create datasystem
    log.info(f"Attempting to load data model for raw directory: {str(root_path)}")
    dataSystem = GLDSDataSystem(base=BaseDataSystem(name=dataSystem_name))

    # load data assets config
    conf_data_assets = load_config(config)

    # create dataset
    dataset = BulkRNASeqDataset(base=BaseDataset(name=dataSystem_name))

    dataSystem.attach_dataset(dataset)
    # attach dataset components
    dataSystem.dataset.attach_component(
        BulkRNASeqMetadataComponent(
            base=BaseComponent(description="Metadata in a runsheet csv file"),
            runsheet=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["runsheet"], dataset=dataSystem_name)),
        ),
        attr="metadata",
    )

    # alias metadata for convenience
    metadata = dataSystem.dataset.all_components["metadata"]

    # create shared sample datafiles
    datf_readsMQC = DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw MultiQC directory ZIP"], dataset=dataSystem_name))

    # create samples
    for sample_name in metadata.samples:
        sample = BulkRNASeqSample(BaseSample(name=sample_name))
        if metadata.paired_end:
            raw_fwd_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Forward Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw forward reads fastq GZ"], sample=sample_name)),
                multiQCDirZIP=datf_readsMQC,
                fastqcReportHTML=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw forward reads fastQC HTML"], sample=sample_name)
                ),
                fastqcReportZIP=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw forward reads fastQC ZIP"], sample=sample_name)
                ),
            )
            raw_rev_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Reverse Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw reverse reads fastq GZ"], sample=sample_name)),
                multiQCDirZIP=datf_readsMQC,
                fastqcReportHTML=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw reverse reads fastQC HTML"], sample=sample_name)),
                fastqcReportZIP=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw reverse reads fastQC ZIP"], sample=sample_name))
            )
            log.debug(f"Attaching components to {sample.name}")
            sample.attach_component(raw_fwd_reads, attr="rawForwardReads")
            sample.attach_component(raw_rev_reads, attr="rawReverseReads")
        else:
            raw_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw reads fastq GZ"], sample=sample_name)),
                multiQCDirZIP=datf_readsMQC,
                fastqcReportHTML=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw reads fastQC HTML"], sample=sample_name)),
                fastqcReportZIP=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw reads fastQC ZIP"], sample=sample_name))
            )
            log.debug(f"Attaching components to {sample.name}")
            sample.attach_component(raw_reads, attr="rawReads")
        # attach components
        log.debug(f"Attaching {sample.name} to {dataset.name}")
        dataset.attach_sample(sample)

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


# TODO: docstring
def load_BulkRNASeq_STAGE_01(
    root_path: Path, dataSystem: TemplateDataSystem, config: Union[str, Path] = "Latest", stack: bool = False
):
    """Load data that should be present into stage 01 (TG-PreProc)

    :param dataSystem: The dataSystem as loaded through STAGE 00
    :type dataSystem: TemplateDataSystem
    """
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # log about datasystem being built upon
    log.info(f"Loading STAGE 01 BulkRNASeq data into: \n\t{dataSystem}")

    # load data assets config
    conf_data_assets = load_config(config)

    # attach dataset components
    # dataSystem.dataset.attach_component(
    #     BulkRNASeqMetadataComponent(
    #         base=BaseComponent(description="Metadata in a runsheet csv file"),
    #         runsheet=DataFile(runsheet.find(datasystem_name=dataSystem_name)),
    #     ),
    #     attr="metadata",
    # )

    # alias for convenience
    metadata = dataSystem.dataset.all_components["metadata"]
    dataset = dataSystem.dataset

    # create shared sample datafiles
    datf_readsMQC = DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed fastQC MultiQC directory ZIP"]))
    datf_trimmingMQC = DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimming MultiQC directory ZIP"]))

    # update samples
    for sample_name, sample in dataset.samples.items():
        if metadata.paired_end:
            trim_fwd_reads = TrimReadsComponent(
                base=BaseComponent(description="Trimmed Forward Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed forward reads fastq GZ"], sample=sample_name)),
                fastQCmultiQCDirZIP=datf_readsMQC,
                trimmingMultiQCDirZIP=datf_trimmingMQC,
                fastqcReportHTML=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed forward reads fastQC HTML"], sample=sample_name)),
                fastqcReportZIP=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed forward reads fastQC ZIP"], sample=sample_name)),
                trimmingReportTXT=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["forward reads trimming report"], sample=sample_name)),
            )
            trim_rev_reads = TrimReadsComponent(
                base=BaseComponent(description="Trimmed Reverse Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed reverse reads fastq GZ"], sample=sample_name)),
                fastQCmultiQCDirZIP=datf_readsMQC,
                trimmingMultiQCDirZIP=datf_trimmingMQC,
                fastqcReportHTML=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed reverse reads fastQC HTML"], sample=sample_name)),
                fastqcReportZIP=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed reverse reads fastQC ZIP"], sample=sample_name)),
                trimmingReportTXT=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["reverse reads trimming report"], sample=sample_name)),
            )
            sample.attach_component(trim_fwd_reads, attr="trimForwardReads")
            sample.attach_component(trim_rev_reads, attr="trimReverseReads")
        else:
            trim_reads = TrimReadsComponent(
                base=BaseComponent(description="Trimmed Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed reads fastq GZ"], sample=sample_name)),
                fastQCmultiQCDirZIP=datf_readsMQC,
                trimmingMultiQCDirZIP=datf_trimmingMQC,
                fastqcReportHTML=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed reads fastQC HTML"], sample=sample_name)),
                fastqcReportZIP=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed reads fastQC ZIP"], sample=sample_name)),
                trimmingReportTXT=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["reads trimming report"], sample=sample_name)),
            )
            sample.attach_component(trim_reads, attr="trimReads")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_02(
    root_path: Path, dataSystem: TemplateDataSystem, config: Union[str, Path] = "Latest", stack: bool = False
):
    """Load data that should be present into stage 02 (GenomeAlignment)

    :param dataSystem: The dataSystem as loaded through STAGE 02
    :type dataSystem: TemplateDataSystem
    """
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # log about datasystem being built upon
    log.info(f"Loading STAGE 02 BulkRNASeq data into: \n\t{dataSystem}")

    # load data assets config
    conf_data_assets = load_config(config)

    # alias for convenience
    dataset = dataSystem.dataset

    # create shared sample datafiles
    datf_alignMQC = DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["aligned MultiQC directory ZIP"]))

    # update samples
    for sample_name, sample in dataset.samples.items():
        # attach txBam
        genomeAlignments = GenomeAlignments(
            base=BaseComponent(description="Genome alignments"),
            alignedToTranscriptomeBam=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["aligned ToTranscriptome Bam"], sample=sample_name)),
            alignedSortedByCoordBam=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["aligned SortedByCoord Bam"], sample=sample_name)),
            alignedSortedByCoordResortedBam=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["aligned SortedByCoord ResortedBam"], sample=sample_name)),
            alignedSortedByCoordResortedBamIndex=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["aligned SortedByCoord ResortedBamIndex"], sample=sample_name)),
            logFinal=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["aligned log Final"], sample=sample_name)),
            logProgress=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["aligned log Progress"], sample=sample_name)),
            logFull=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["aligned log Full"], sample=sample_name)),
            sjTab=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["aligned sjTab"], sample=sample_name)),
            multiQCDirZIP=datf_alignMQC,
        )
        sample.attach_component(genomeAlignments, attr="genomeAlignments")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_0201(
    root_path: Path, dataSystem: TemplateDataSystem, config: Union[str, Path] = "Latest", stack: bool = False
):
    """Load data that should be present into stage 0201 (RSeQCAnalysis)

    :param dataSystem: The dataSystem as loaded through STAGE 0201
    :type dataSystem: TemplateDataSystem
    """
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # log about datasystem being built upon
    log.info(f"Loading STAGE 0201 BulkRNASeq data into: \n\t{dataSystem}")

    # load data assets config
    conf_data_assets = load_config(config)

    # alias for convenience
    dataset = dataSystem.dataset
    metadata = dataset.metadata

    # create shared sample datafiles
    genebody_coverage_MQC = DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["genebody coverage MultiQC directory ZIP"]))
    infer_experiment_MQC = DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["infer experiment MultiQC directory ZIP"]))
    # idMQC moved to metadata controlled block
    read_distribution_MQC = DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["read distribution MultiQC directory ZIP"]))

    # update samples
    for sample_name, sample in dataset.samples.items():
        # create common component
        rSeQCAnalysis = RSeQCAnalysis(
            base=BaseComponent(
                description="RSeQC Analysis based on reads aligned to genome in context of gene annotations"
            ),
            geneBodyCoverageMultiQCDirZIP=genebody_coverage_MQC,
            geneBodyCoverageOut=DataDir(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["genebody coverage out"], sample=sample_name)),
            inferExperimentMultiQCDirZIP=infer_experiment_MQC,
            inferExperimentOut=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["infer experiment out"], sample=sample_name)),
            readDistributionMultiQCDirZIP=read_distribution_MQC,
            readDistributionOut=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["read distribution out"], sample=sample_name)),
        )
        # add pair end specific datafile
        if metadata.paired_end:
            rSeQCAnalysis.innerDistanceMultiQCDirZIP = DataDir(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["inner distance MultiQC directory"], sample=sample_name))
            rSeQCAnalysis.innerDistanceOut = DataDir(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["inner distance out"], sample=sample_name))
        sample.attach_component(rSeQCAnalysis, attr="rSeQCAnalysis")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_03(
    root_path: Path, dataSystem: TemplateDataSystem, config: Union[str, Path] = "Latest", stack: bool = False
):
    """Load data that should be present into stage 03 (RSEM Counts)

    :param dataSystem: The dataSystem as loaded through STAGE 03
    :type dataSystem: TemplateDataSystem
    """
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # log about datasystem being built upon
    log.info(f"Loading STAGE 03 BulkRNASeq data into: \n\t{dataSystem}")

    # load data assets config
    conf_data_assets = load_config(config)

    # alias for convenience
    dataset = dataSystem.dataset
    metadata = dataset.metadata

    # create shared sample datafiles
    countsMQC = DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["RSEM counts MultiQC directory ZIP"]))

    # update dataset
    dataset.geneCounts = DatasetGeneCounts(
        base=BaseComponent(
            description="Gene counts at a dataset level from RSEM and DESeq2"
        ),
        numNonZero=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["number non-zero count genes table"])),
        unnormalizedCounts=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["unnormalized counts table"])),
    )

    # update samples
    for sample_name, sample in dataset.samples.items():
        # create common component
        geneCounts = GeneCounts(
            base=BaseComponent(description="Gene counts for the sample"),
            multiQCDirZIP=countsMQC,
            genesResults=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["sample gene counts table"], sample=sample_name)),
            isoformsResults=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["sample isoform counts table"], sample=sample_name)),
            statDir=DataDir(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["sample counts stats directory"], sample=sample_name)),
        )
        sample.attach_component(geneCounts, attr="geneCounts")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_04(
    root_path: Path, dataSystem: TemplateDataSystem, config: Union[str, Path] = "Latest", stack: bool = False
):
    """Load data that should be present into stage 03 (RSEM Counts)

    :param dataSystem: The dataSystem as loaded through STAGE 03
    :type dataSystem: TemplateDataSystem
    """
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # log about datasystem being built upon
    log.info(f"Loading STAGE 03 BulkRNASeq data into: \n\t{dataSystem}")

    # load data assets config
    conf_data_assets = load_config(config)

    # alias for convenience
    dataset = dataSystem.dataset

    # create shared sample datafiles
    ...

    # update dataset
    # first
    dataset.attach_component(
        NormalizedGeneCounts(
            base=BaseComponent(description="Normalized counts from Deseq2"),
            # erccNormalizedCountsCSV=DataFile(enc.find()), moved to ERCC block
            normalizedCountsCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["DESeq2 normalized counts table"])),
            sampleTableCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["sample table"])),
            unnormalizedCountsCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["DESeq2 unnormalized counts table"])),
        ),
        attr="normalizedGeneCounts",
    )
    dataset.attach_component(
        DifferentialGeneExpression(
            base=BaseComponent(description="Normalized counts from Deseq2"),
            contrastsCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["DESeq2 contrasts table"])),
            annotatedTableCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["DESeq2 annotated DGE table"])),
            visualizationTableCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["DESeq2 annotated DGE extended for viz table"])),
            visualizationPCATableCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["DESeq2 viz PCA table"])),
        ),
        attr="differentialGeneExpression",
    )
    if dataset.metadata.has_ercc:
        dataset.attach_component(
            DifferentialGeneExpression(
                base=BaseComponent(description="Normalized counts from Deseq2"),
                contrastsCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["ERCC normalized DESeq2 contrasts table"])),
                annotatedTableCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["ERCC normalized DESeq2 annotated DGE table"])),
                visualizationTableCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["ERCC normalized DESeq2 annotated DGE extended for viz table"])),
                visualizationPCATableCSV=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["ERCC normalized DESeq2 viz PCA table"])),
            ),
            attr="differentialGeneExpressionERCC",
        )
        dataset.normalizedGeneCounts.erccNormalizedCountsCSV = DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["ERCC normalized DESeq2 normalized counts table"])),

    # update samples
    ...

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem
