""" Functions for loading into data model objects 

Loaders SHOULD:
  - run validation on the dataSystem object
  - return a dataSystem object
"""
import functools
import os
from pathlib import Path
import logging
from typing import Union
import pkg_resources

from schema import Schema, Or, SchemaMissingKeyError, Optional
import pandas as pd
from dp_tools.components.components import (
    DatasetGeneCounts,
    DifferentialGeneExpression,
    ERCCAnalysis,
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
from dp_tools.bulkRNASeq.locaters import Locator
from dp_tools.components import RawReadsComponent, BulkRNASeqMetadataComponent
from dp_tools.core.configuration import load_full_config


@functools.cache
def load_config(config: Union[str, Path]) -> dict:
    if isinstance(config, str):
        resolved_config_path = os.path.join(
            "..", "config", f"bulkRNASeq_v{config}.yaml"
        )
        log.info(f"Loading config (relative to package): {resolved_config_path}")
        conf_data_assets = yaml.safe_load(
            pkg_resources.resource_string(__name__, resolved_config_path)
        )["data assets"]
    elif isinstance(config, Path):
        log.info(f"Loading config (direct path): {config}")
        conf_data_assets = yaml.safe_load(config.open())["data assets"]

    # validate with schema
    config_schema = Schema(
        {
            # a list of paths designating the
            # expected full data asset path
            "processed location": [str],
            # used to control how curation tables are formatted
            # and what data assets are included
            "resource categories": {
                # top level directory in S3 archive,
                # first and required Parameter Value name component in assay table
                "subcategory": Or(None, str),
                # second level directory in S3 archive,
                # second and optional Parameter Value name component in assay table
                "subdirectory": Or(None, str),
                # controls publishing data asset to GLDS repo
                # controls inclusion of data asset in assay table
                "publish to repo": bool,
                # controls if the subdirectory label should be included in the assay table
                "include subdirectory in table": bool,
                # controls the order of PROCESSED Parameter Value columns added to
                # existing RAW DATA assay table
                "table order": int,
            },
            # used to denote human description of the data asset
            # e.g. this is 'raw' data, or this is 'experimental' data
            # Can be leveraged by downstream applications
            # An example here includes dividing the asset
            # into different tables
            "tags": [str],
            # controls validation of expected data asset path,
            # disabling this means the path will not be checked for existence
            # some example use cases:
            #   - a table needs to be updated and the original
            #       local data is not readily available
            #
            #   - a dry run needs to be performed during development
            Optional("validate exists"): bool,
        }
    )
    for key, conf_data_asset in conf_data_assets.items():
        try:
            config_schema.validate(conf_data_asset)
        except SchemaMissingKeyError as e:
            raise ValueError(f"Data asset config: '{key}' failed validation") from e

    log.debug(f"Final config loaded: {conf_data_assets}")

    return conf_data_assets


# TODO: Attach/associate data assets config with dataset/datasystem
def load_BulkRNASeq_STAGE_00(
    root_path: Path,
    config: Union[str, Path] = "Latest",
    dataSystem_name: str = None,
    stack: bool = False,
    validation_enabled: bool = True,
):
    """config is either an string version referring to a packaged config file or a path to a local config file."""
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
    dataset = BulkRNASeqDataset(
        base=BaseDataset(name=dataSystem_name, config=load_full_config(config))
    )

    # create locator
    loc = Locator(
        root_dir=root_path,
        data_asset_config=conf_data_assets,
        validation_enabled=validation_enabled,
    )

    dataSystem.attach_dataset(dataset)
    # attach dataset components
    dataSystem.dataset.attach_component(
        BulkRNASeqMetadataComponent(
            base=BaseComponent(description="Metadata in a runsheet csv file"),
            runsheet=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="runsheet",
                    dataset=dataSystem_name,
                ),
            ),
            ISAarchive=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="ISA Archive",
                    glob=True,
                    dataset=dataSystem_name,
                ),
            ),
        ),
        attr="metadata",
    )

    # alias metadata for convenience
    metadata = dataSystem.dataset.all_components["metadata"]

    # create shared sample datafiles
    datf_readsMQC = DataFile(
        check_exists=validation_enabled,
        **loc.find_data_asset_path(
            config_key="raw MultiQC directory ZIP",
            dataset=dataSystem_name,
        ),
    )

    # create samples
    for sample_name in metadata.samples:
        sample = BulkRNASeqSample(BaseSample(name=sample_name))
        if metadata.paired_end:
            raw_fwd_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Forward Reads"),
                fastqGZ=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="raw forward reads fastq GZ",
                        sample=sample_name,
                    ),
                ),
                fastQCmultiQCDirZIP=datf_readsMQC,
                fastqcReportHTML=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="raw forward reads fastQC HTML",
                        sample=sample_name,
                    ),
                ),
                fastqcReportZIP=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="raw forward reads fastQC ZIP",
                        sample=sample_name,
                    ),
                ),
            )
            raw_rev_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Reverse Reads"),
                fastqGZ=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="raw reverse reads fastq GZ",
                        sample=sample_name,
                    ),
                ),
                fastQCmultiQCDirZIP=datf_readsMQC,
                fastqcReportHTML=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="raw reverse reads fastQC HTML",
                        sample=sample_name,
                    ),
                ),
                fastqcReportZIP=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="raw reverse reads fastQC ZIP",
                        sample=sample_name,
                    ),
                ),
            )
            log.debug(f"Attaching components to {sample.name}")
            sample.attach_component(raw_fwd_reads, attr="rawForwardReads")
            sample.attach_component(raw_rev_reads, attr="rawReverseReads")
        else:
            raw_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Reads"),
                fastqGZ=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="raw reads fastq GZ",
                        sample=sample_name,
                    ),
                ),
                fastQCmultiQCDirZIP=datf_readsMQC,
                fastqcReportHTML=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="raw reads fastQC HTML",
                        sample=sample_name,
                    ),
                ),
                fastqcReportZIP=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="raw reads fastQC ZIP",
                        sample=sample_name,
                    ),
                ),
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
    root_path: Path,
    dataSystem: TemplateDataSystem,
    config: Union[str, Path] = "Latest",
    stack: bool = False,
    validation_enabled: bool = True,
):
    """Load data that should be present into stage 01 (TG-PreProc)

    :param dataSystem: The dataSystem as loaded through STAGE 00
    :type dataSystem: TemplateDataSystem
    """
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # log about datasystem being built upon
    log.info(
        f"Loading STAGE 01 BulkRNASeq data into: {dataSystem.name}:{dataSystem.dataset.name}"
    )

    # load data assets config
    conf_data_assets = load_config(config)

    # create locator
    loc = Locator(
        root_dir=root_path,
        data_asset_config=conf_data_assets,
        validation_enabled=validation_enabled,
    )

    # attach dataset components
    # dataSystem.dataset.attach_component(
    #     BulkRNASeqMetadataComponent(
    #         base=BaseComponent(description="Metadata in a runsheet csv file"),
    #         runsheet=DataFile(check_exists=validation_enabled,runsheet.find(datasystem_name=dataSystem_name)),
    #     ),
    #     attr="metadata",
    # )

    # alias for convenience
    metadata = dataSystem.dataset.all_components["metadata"]
    dataset = dataSystem.dataset

    # create shared sample datafiles
    datf_readsMQC = DataFile(
        check_exists=validation_enabled,
        **loc.find_data_asset_path(
            config_key="trimmed fastQC MultiQC directory ZIP",
        ),
    )
    """ Archived data asset
    datf_trimmingMQC = DataFile(check_exists=validation_enabled,
        **loc.find_data_asset_path(config_key = "trimming MultiQC directory ZIP"],
        )
    )
    """

    # update samples
    for sample_name, sample in dataset.samples.items():
        if metadata.paired_end:
            trim_fwd_reads = TrimReadsComponent(
                base=BaseComponent(description="Trimmed Forward Reads"),
                fastqGZ=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="trimmed forward reads fastq GZ",
                        sample=sample_name,
                    ),
                ),
                fastQCmultiQCDirZIP=datf_readsMQC,
                # Archived data asset trimmingMultiQCDirZIP=datf_trimmingMQC,
                fastqcReportHTML=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="trimmed forward reads fastQC HTML",
                        sample=sample_name,
                    ),
                ),
                fastqcReportZIP=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="trimmed forward reads fastQC ZIP",
                        sample=sample_name,
                    ),
                ),
                trimmingReportTXT=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="forward reads trimming report",
                        sample=sample_name,
                    ),
                ),
            )
            trim_rev_reads = TrimReadsComponent(
                base=BaseComponent(description="Trimmed Reverse Reads"),
                fastqGZ=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="trimmed reverse reads fastq GZ",
                        sample=sample_name,
                    ),
                ),
                fastQCmultiQCDirZIP=datf_readsMQC,
                # Archived data asset trimmingMultiQCDirZIP=datf_trimmingMQC,
                fastqcReportHTML=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="trimmed reverse reads fastQC HTML",
                        sample=sample_name,
                    ),
                ),
                fastqcReportZIP=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="trimmed reverse reads fastQC ZIP",
                        sample=sample_name,
                    ),
                ),
                trimmingReportTXT=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="reverse reads trimming report",
                        sample=sample_name,
                    ),
                ),
            )
            sample.attach_component(trim_fwd_reads, attr="trimForwardReads")
            sample.attach_component(trim_rev_reads, attr="trimReverseReads")
        else:
            trim_reads = TrimReadsComponent(
                base=BaseComponent(description="Trimmed Reads"),
                fastqGZ=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="trimmed reads fastq GZ",
                        sample=sample_name,
                    ),
                ),
                fastQCmultiQCDirZIP=datf_readsMQC,
                # Archived data asset trimmingMultiQCDirZIP=datf_trimmingMQC,
                fastqcReportHTML=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="trimmed reads fastQC HTML",
                        sample=sample_name,
                    ),
                ),
                fastqcReportZIP=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="trimmed reads fastQC ZIP",
                        sample=sample_name,
                    ),
                ),
                trimmingReportTXT=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="reads trimming report",
                        sample=sample_name,
                    ),
                ),
            )
            sample.attach_component(trim_reads, attr="trimReads")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_02(
    root_path: Path,
    dataSystem: TemplateDataSystem,
    config: Union[str, Path] = "Latest",
    stack: bool = False,
    validation_enabled: bool = True,
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

    # create locator
    loc = Locator(
        root_dir=root_path,
        data_asset_config=conf_data_assets,
        validation_enabled=validation_enabled,
    )

    # alias for convenience
    dataset = dataSystem.dataset

    # create shared sample datafiles
    datf_alignMQC = DataFile(
        check_exists=validation_enabled,
        **loc.find_data_asset_path(
            config_key="aligned MultiQC directory ZIP",
        ),
    )

    # update samples
    for sample_name, sample in dataset.samples.items():
        # attach txBam
        genomeAlignments = GenomeAlignments(
            base=BaseComponent(description="Genome alignments"),
            alignedToTranscriptomeBam=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="aligned ToTranscriptome Bam",
                    sample=sample_name,
                ),
            ),
            alignedSortedByCoordBam=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="aligned SortedByCoord Bam",
                    sample=sample_name,
                ),
            ),
            alignedSortedByCoordResortedBam=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="aligned SortedByCoord ResortedBam",
                    sample=sample_name,
                ),
            ),
            alignedSortedByCoordResortedBamIndex=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="aligned SortedByCoord ResortedBamIndex",
                    sample=sample_name,
                ),
            ),
            logFinal=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="aligned log Final",
                    sample=sample_name,
                ),
            ),
            logProgress=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="aligned log Progress",
                    sample=sample_name,
                ),
            ),
            logFull=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="aligned log Full",
                    sample=sample_name,
                ),
            ),
            sjTab=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="aligned sjTab",
                    sample=sample_name,
                ),
            ),
            readsPerGeneTable=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="sample reads per gene table",
                    sample=sample_name,
                ),
            ),
            multiQCDirZIP=datf_alignMQC,
        )
        sample.attach_component(genomeAlignments, attr="genomeAlignments")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_0201(
    root_path: Path,
    dataSystem: TemplateDataSystem,
    config: Union[str, Path] = "Latest",
    stack: bool = False,
    validation_enabled: bool = True,
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

    # create locator
    loc = Locator(
        root_dir=root_path,
        data_asset_config=conf_data_assets,
        validation_enabled=validation_enabled,
    )

    # alias for convenience
    dataset = dataSystem.dataset
    metadata = dataset.metadata

    # create shared sample datafiles
    genebody_coverage_MQC = DataFile(
        check_exists=validation_enabled,
        **loc.find_data_asset_path(
            config_key="genebody coverage MultiQC directory ZIP",
        ),
    )
    infer_experiment_MQC = DataFile(
        check_exists=validation_enabled,
        **loc.find_data_asset_path(
            config_key="infer experiment MultiQC directory ZIP",
        ),
    )
    # idMQC moved to metadata controlled block
    read_distribution_MQC = DataFile(
        check_exists=validation_enabled,
        **loc.find_data_asset_path(
            config_key="read distribution MultiQC directory ZIP",
        ),
    )

    # update samples
    for sample_name, sample in dataset.samples.items():
        # create common component
        rSeQCAnalysis = RSeQCAnalysis(
            base=BaseComponent(
                description="RSeQC Analysis based on reads aligned to genome in context of gene annotations"
            ),
            geneBodyCoverageMultiQCDirZIP=genebody_coverage_MQC,
            geneBodyCoverageOut=DataDir(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="genebody coverage out",
                    sample=sample_name,
                ),
            ),
            inferExperimentMultiQCDirZIP=infer_experiment_MQC,
            inferExperimentOut=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="infer experiment out",
                    sample=sample_name,
                ),
            ),
            readDistributionMultiQCDirZIP=read_distribution_MQC,
            readDistributionOut=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="read distribution out",
                    sample=sample_name,
                ),
            ),
        )
        # add pair end specific datafile
        if metadata.paired_end:
            rSeQCAnalysis.innerDistanceMultiQCDirZIP = DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="inner distance MultiQC directory ZIP",
                    sample=sample_name,
                ),
            )
            rSeQCAnalysis.innerDistanceOut = DataDir(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="inner distance out",
                    sample=sample_name,
                ),
            )
        sample.attach_component(rSeQCAnalysis, attr="rSeQCAnalysis")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_03(
    root_path: Path,
    dataSystem: TemplateDataSystem,
    config: Union[str, Path] = "Latest",
    stack: bool = False,
    validation_enabled: bool = True,
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

    # create locator
    loc = Locator(
        root_dir=root_path,
        data_asset_config=conf_data_assets,
        validation_enabled=validation_enabled,
    )

    # alias for convenience
    dataset = dataSystem.dataset
    metadata = dataset.metadata

    # create shared sample datafiles
    countsMQC = DataFile(
        check_exists=validation_enabled,
        **loc.find_data_asset_path(
            config_key="RSEM counts MultiQC directory ZIP",
        ),
    )

    # update dataset
    dataset.rsemGeneCounts = DatasetGeneCounts(
        base=BaseComponent(
            description="Gene counts at a dataset level from RSEM and DESeq2"
        ),
        numNonZero=DataFile(
            check_exists=validation_enabled,
            **loc.find_data_asset_path(
                config_key="rsem number non-zero count genes table",
            ),
        ),
        unnormalizedCounts=DataFile(
            check_exists=validation_enabled,
            **loc.find_data_asset_path(
                config_key="rsem unnormalized counts table",
            ),
        ),
    )

    dataset.starGeneCounts = DatasetGeneCounts(
        base=BaseComponent(
            description="Gene counts at a dataset level from RSEM and DESeq2"
        ),
        numNonZero=DataFile(
            check_exists=validation_enabled,
            **loc.find_data_asset_path(
                config_key="star number non-zero count genes table",
            ),
        ),
        unnormalizedCounts=DataFile(
            check_exists=validation_enabled,
            **loc.find_data_asset_path(
                config_key="star unnormalized counts table",
            ),
        ),
    )

    # update samples
    for sample_name, sample in dataset.samples.items():
        # create common component
        geneCounts = GeneCounts(
            base=BaseComponent(description="Gene counts for the sample"),
            multiQCDirZIP=countsMQC,
            genesResults=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="sample gene counts table",
                    sample=sample_name,
                ),
            ),
            isoformsResults=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="sample isoform counts table",
                    sample=sample_name,
                ),
            ),
            statDir=DataDir(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="sample counts stats directory",
                    sample=sample_name,
                ),
            ),
        )
        sample.attach_component(geneCounts, attr="geneCounts")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_04(
    root_path: Path,
    dataSystem: TemplateDataSystem,
    config: Union[str, Path] = "Latest",
    stack: bool = False,
    validation_enabled: bool = True,
):
    """Load data that should be present into stage 03 (RSEM Counts)

    :param dataSystem: The dataSystem as loaded through STAGE 03
    :type dataSystem: TemplateDataSystem
    """
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # log about datasystem being built upon
    log.info(f"Loading STAGE 04 BulkRNASeq data into: \n\t{dataSystem}")

    # load data assets config
    conf_data_assets = load_config(config)

    # create locator
    loc = Locator(
        root_dir=root_path,
        data_asset_config=conf_data_assets,
        validation_enabled=validation_enabled,
    )

    # alias for convenience
    dataset = dataSystem.dataset

    # create shared sample datafiles
    ...

    # update dataset
    # first
    dataset.attach_component(
        NormalizedGeneCounts(
            base=BaseComponent(description="Normalized counts from Deseq2"),
            # erccNormalizedCountsCSV=DataFile(check_exists=validation_enabled,enc.find()), moved to ERCC block
            normalizedCountsCSV=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(config_key="DESeq2 normalized counts table"),
            ),
            sampleTableCSV=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="sample table",
                ),
            ),
            unnormalizedCountsCSV=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="DESeq2 unnormalized counts table",
                ),
            ),
        ),
        attr="normalizedGeneCounts",
    )
    dataset.attach_component(
        DifferentialGeneExpression(
            base=BaseComponent(description="Normalized counts from Deseq2"),
            contrastsCSV=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="DESeq2 contrasts table",
                ),
            ),
            annotatedTableCSV=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="DESeq2 annotated DGE table",
                ),
            ),
            visualizationTableCSV=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="DESeq2 annotated DGE extended for viz table",
                ),
            ),
            visualizationPCATableCSV=DataFile(
                check_exists=validation_enabled,
                **loc.find_data_asset_path(
                    config_key="DESeq2 viz PCA table",
                ),
            ),
        ),
        attr="differentialGeneExpression",
    )
    if dataset.metadata.has_ercc:
        dataset.attach_component(
            DifferentialGeneExpression(
                base=BaseComponent(description="Normalized counts from Deseq2"),
                contrastsCSV=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="ERCC normalized DESeq2 contrasts table",
                    ),
                ),
                annotatedTableCSV=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="ERCC normalized DESeq2 annotated DGE table",
                    ),
                ),
                visualizationTableCSV=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="ERCC normalized DESeq2 annotated DGE extended for viz table",
                    ),
                ),
                visualizationPCATableCSV=DataFile(
                    check_exists=validation_enabled,
                    **loc.find_data_asset_path(
                        config_key="ERCC normalized DESeq2 viz PCA table",
                    ),
                ),
            ),
            attr="differentialGeneExpressionERCC",
        )
        dataset.normalizedGeneCounts.erccNormalizedCountsCSV = DataFile(
            check_exists=validation_enabled,
            **loc.find_data_asset_path(
                config_key="ERCC normalized DESeq2 normalized counts table",
            ),
        )
        dataset.normalizedGeneCounts.erccSampleTableCSV = DataFile(
            check_exists=validation_enabled,
            **loc.find_data_asset_path(
                config_key="ERCC sample table",
            ),
        )
        dataset.attach_component(
            ERCCAnalysis(
                base=BaseComponent(description="ERCC Notebook analysis related files"),
                # WARNING: This is the only data asset expected to be generated outside the workflow at this time
                # as such it is assumed the file will exist; however, neither the locator nor the data asset itself will validate this
                notebookHTML=DataFile(
                    check_exists=False,
                    **loc.find_data_asset_path(
                        config_key="ERCC analysis HTML", search=False
                    ),
                ),
            ),
            attr="erccAnalysis",
        )

    # update samples
    ...

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem
