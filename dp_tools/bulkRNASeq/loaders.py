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
    datf_readsMQC = DataDir(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw MultiQC directory"], dataset=dataSystem_name))

    # create samples
    for sample_name in metadata.samples:
        sample = BulkRNASeqSample(BaseSample(name=sample_name))
        if metadata.paired_end:
            raw_fwd_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Forward Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw forward reads fastq GZ"], sample=sample_name)),
                multiQCDir=datf_readsMQC,
                fastqcReportHTML=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw forward reads fastQC HTML"], sample=sample_name)
                ),
                fastqcReportZIP=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw forward reads fastQC ZIP"], sample=sample_name)
                ),
            )
            raw_rev_reads = RawReadsComponent(
                base=BaseComponent(description="Raw Reverse Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["raw reverse reads fastq GZ"], sample=sample_name)),
                multiQCDir=datf_readsMQC,
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
                multiQCDir=datf_readsMQC,
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
    datf_readsMQC = DataDir(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed MultiQC directory"]))

    # update samples
    for sample_name, sample in dataset.samples.items():
        if metadata.paired_end:
            trim_fwd_reads = TrimReadsComponent(
                base=BaseComponent(description="Trimmed Forward Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed forward reads fastq GZ"], sample=sample_name)),
                multiQCDir=datf_readsMQC,
                fastqcReportHTML=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed forward reads fastQC HTML"], sample=sample_name)),
                fastqcReportZIP=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed forward reads fastQC ZIP"], sample=sample_name)),
                trimmingReportTXT=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["forward reads trimming report"], sample=sample_name)),
            )
            trim_rev_reads = TrimReadsComponent(
                base=BaseComponent(description="Trimmed Reverse Reads"),
                fastqGZ=DataFile(find_data_asset_path(root_dir=root_path, data_asset_config=conf_data_assets["trimmed reverse reads fastq GZ"], sample=sample_name)),
                multiQCDir=datf_readsMQC,
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
                multiQCDir=datf_readsMQC,
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
    root_path: Path, dataSystem: TemplateDataSystem, stack: bool = False
):
    """Load data that should be present into stage 02 (GenomeAlignment)

    :param dataSystem: The dataSystem as loaded through STAGE 02
    :type dataSystem: TemplateDataSystem
    """
    # ensure root path exists!
    if not root_path.is_dir():
        raise FileNotFoundError(f"Root path doesn't exist!: {root_path}")

    # log about datasystem being built upon
    log.info(f"Loading STAGE 01 BulkRNASeq data into: \n\t{dataSystem}")

    # initate locaters on the root path
    txBam = AlignedToTranscriptomeBam(search_root=root_path)
    coordBam = AlignedSortedByCoordBam(search_root=root_path)
    coordBamResort = AlignedSortedByCoordResortedBam(search_root=root_path)
    coordBamResortIndex = AlignedSortedByCoordResortedBamIndex(search_root=root_path)
    logFinal = LogFinal(search_root=root_path)
    logProgress = LogProgress(search_root=root_path)
    logFull = LogFull(search_root=root_path)
    sjTab = SjTab(search_root=root_path)
    readsMQC = MultiQCDir(search_root=root_path)

    # alias for convenience
    dataset = dataSystem.dataset

    # create shared sample datafiles
    datf_alignMQC = DataDir(
        readsMQC.find(
            rel_dir=Path.joinpath(Path("02-STAR_Alignment")), mqc_label="align",
        )
    )

    # update samples
    for sample_name, sample in dataset.samples.items():
        # attach txBam
        genomeAlignments = GenomeAlignments(
            base=BaseComponent(description="Genome alignments"),
            alignedToTranscriptomeBam=DataFile(
                path=txBam.find(sample_name=sample_name)
            ),
            alignedSortedByCoordBam=DataFile(
                path=coordBam.find(sample_name=sample_name)
            ),
            alignedSortedByCoordResortedBam=DataFile(
                path=coordBamResort.find(sample_name=sample_name)
            ),
            alignedSortedByCoordResortedBamIndex=DataFile(
                path=coordBamResortIndex.find(sample_name=sample_name)
            ),
            logFinal=DataFile(path=logFinal.find(sample_name=sample_name)),
            logProgress=DataFile(path=logProgress.find(sample_name=sample_name)),
            logFull=DataFile(path=logFull.find(sample_name=sample_name)),
            sjTab=DataFile(path=sjTab.find(sample_name=sample_name)),
            multiQCDir=datf_alignMQC,
        )
        sample.attach_component(genomeAlignments, attr="genomeAlignments")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_0201(
    root_path: Path, dataSystem: TemplateDataSystem, stack: bool = False
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

    # initate locaters on the root path
    gbOut = GeneBodyCoverageOut(search_root=root_path)
    ieOut = InferExperimentOut(search_root=root_path)
    idOut = InnerDistanceOut(search_root=root_path)
    rdOut = ReadDistributionOut(search_root=root_path)
    mQC = MultiQCDir(search_root=root_path)

    # alias for convenience
    dataset = dataSystem.dataset
    metadata = dataset.metadata

    # create shared sample datafiles
    gbMQC = DataDir(
        mQC.find(
            rel_dir=(Path("RSeQC_Analyses") / "02_geneBody_coverage"),
            mqc_label="geneBody_cov",
        )
    )
    ieMQC = DataDir(
        mQC.find(
            rel_dir=(Path("RSeQC_Analyses") / "03_infer_experiment"),
            mqc_label="infer_exp",
        )
    )
    # idMQC moved to metadata controlled block
    rdMQC = DataDir(
        mQC.find(
            rel_dir=(Path("RSeQC_Analyses") / "05_read_distribution"),
            mqc_label="read_dist",
        )
    )

    # update samples
    for sample_name, sample in dataset.samples.items():
        # create common component
        rSeQCAnalysis = RSeQCAnalysis(
            base=BaseComponent(
                description="RSeQC Analysis based on reads aligned to genome in context of gene annotations"
            ),
            geneBodyCoverageMultiQCDir=gbMQC,
            geneBodyCoverageOut=DataDir(path=gbOut.find(sample_name=sample_name)),
            inferExperimentMultiQCDir=ieMQC,
            inferExperimentOut=DataFile(path=ieOut.find(sample_name=sample_name)),
            readDistributionMultiQCDir=rdMQC,
            readDistributionOut=DataFile(path=rdOut.find(sample_name=sample_name)),
        )
        # add pair end specific datafile
        if metadata.paired_end:
            rSeQCAnalysis.innerDistanceMultiQCDir = DataDir(
                mQC.find(
                    rel_dir=(Path("RSeQC_Analyses") / "04_inner_distance"),
                    mqc_label="inner_dist",
                )
            )
            rSeQCAnalysis.innerDistanceOut = DataDir(
                path=idOut.find(sample_name=sample_name)
            )
        sample.attach_component(rSeQCAnalysis, attr="rSeQCAnalysis")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_03(
    root_path: Path, dataSystem: TemplateDataSystem, stack: bool = False
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

    # initate locaters on the root path
    gr = GenesResults(search_root=root_path)
    ir = IsoformsResults(search_root=root_path)
    rs = RSEMStat(search_root=root_path)
    mQC = MultiQCDir(search_root=root_path)
    # for dataset usage
    nnz = NumNonZero(search_root=root_path)
    ruc = RSEMUnnormalizedCounts(search_root=root_path)

    # alias for convenience
    dataset = dataSystem.dataset
    metadata = dataset.metadata

    # create shared sample datafiles
    countsMQC = DataDir(mQC.find(rel_dir=Path("03-RSEM_Counts"), mqc_label="RSEM_count",))

    # update dataset
    dataset.geneCounts = DatasetGeneCounts(
        base=BaseComponent(
            description="Gene counts at a dataset level from RSEM and DESeq2"
        ),
        numNonZero=DataFile(nnz.find()),
        unnormalizedCounts=DataFile(ruc.find()),
    )

    # update samples
    for sample_name, sample in dataset.samples.items():
        # create common component
        geneCounts = GeneCounts(
            base=BaseComponent(description="Gene counts for the sample"),
            multiQCDir=countsMQC,
            genesResults=DataFile(gr.find(sample_name=sample_name)),
            isoformsResults=DataFile(ir.find(sample_name=sample_name)),
            statDir=DataDir(rs.find(sample_name=sample_name)),
        )
        sample.attach_component(geneCounts, attr="geneCounts")

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem


def load_BulkRNASeq_STAGE_04(
    root_path: Path, dataSystem: TemplateDataSystem, stack: bool = False
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

    # initate locaters on the root path
    enc = ERCCNormalizedCountsCSV(search_root=root_path)
    nc = NormalizedCountsCSV(search_root=root_path)
    st = SampleTableCSV(search_root=root_path)
    uc = UnnormalizedCountsCSV(search_root=root_path)
    contrasts = ContrastsCSV(search_root=root_path)
    at = AnnotatedTableCSV(search_root=root_path)
    vt = VisualizationTableCSV(search_root=root_path)
    vpt = VisualizationPCATableCSV(search_root=root_path)

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
            normalizedCountsCSV=DataFile(nc.find()),
            sampleTableCSV=DataFile(st.find()),
            unnormalizedCountsCSV=DataFile(uc.find()),
        ),
        attr="normalizedGeneCounts",
    )
    dataset.attach_component(
        DifferentialGeneExpression(
            base=BaseComponent(description="Normalized counts from Deseq2"),
            contrastsCSV=DataFile(contrasts.find()),
            annotatedTableCSV=DataFile(at.find()),
            visualizationTableCSV=DataFile(vt.find()),
            visualizationPCATableCSV=DataFile(vpt.find()),
        ),
        attr="differentialGeneExpression",
    )
    if dataset.metadata.has_ercc:
        dataset.attach_component(
            DifferentialGeneExpression(
                base=BaseComponent(description="Normalized counts from Deseq2"),
                contrastsCSV=DataFile(
                    contrasts.find(nested_dir="ERCC_NormDGE", prefix="ERCCnorm_")
                ),
                annotatedTableCSV=DataFile(
                    at.find(nested_dir="ERCC_NormDGE", prefix="ERCCnorm_")
                ),
                visualizationTableCSV=DataFile(
                    vt.find(nested_dir="ERCC_NormDGE", suffix="_ERCCnorm")
                ),
                visualizationPCATableCSV=DataFile(
                    vpt.find(nested_dir="ERCC_NormDGE", suffix="_ERCCnorm")
                ),
            ),
            attr="differentialGeneExpressionERCC",
        )
        dataset.normalizedGeneCounts.erccNormalizedCountsCSV = DataFile(enc.find())

    # update samples
    ...

    # return dataSystem only if not being used in a loading stack
    if stack:
        return root_path, dataSystem
    else:
        return dataSystem
