import os
from pathlib import Path
from dp_tools.bulkRNASeq.loaders import load_BulkRNASeq_STAGE_00
import pytest


def test_bulkRNASeq_STAGE00_Components_paired(glds194_dataSystem_STAGE00):
    # iterate through all raw read components
    for sample in glds194_dataSystem_STAGE00.dataset.samples.values():

        fwdReads = sample.rawForwardReads
        revReads = sample.rawReverseReads

        # test general stats from multiqc
        assert fwdReads.mqcData["FastQC"]["General_Stats"]["total_sequences"] == 2000
        assert revReads.mqcData["FastQC"]["General_Stats"]["total_sequences"] == 2000
        assert fwdReads.mqcData["FastQC"]["General_Stats"]["total_sequences"] == 2000
        assert revReads.mqcData["FastQC"]["General_Stats"]["total_sequences"] == 2000


def test_bulkRNASeq_STAGE00_Components_single(glds48_dataSystem_STAGE00):
    # iterate through all raw read components
    for sample in glds48_dataSystem_STAGE00.dataset.samples.values():

        reads = sample.rawReads

        # test general stats from multiqc
        assert reads.mqcData["FastQC"]["General_Stats"]["total_sequences"] == 2000

        # test plot data from multiqc
        # should return a dataframe representing plot data for the sample
        assert set(reads.mqcData["FastQC"]["Plots"].keys()).issuperset(
            {
                "Mean Quality Scores",
                "Overrepresented sequences",
                "Per Base N Content",
                "Per Sequence GC Content",
                "Per Sequence Quality Scores",
                "Sequence Counts",
                "Sequence Duplication Levels",
            }
        )


def test_bulkRNASeq_STAGE01_Components_paired(glds194_dataSystem_STAGE01):
    # iterate through all raw read components
    for sample in glds194_dataSystem_STAGE01.dataset.samples.values():

        fwdReads = sample.trimForwardReads
        revReads = sample.trimReverseReads

        # test general stats from multiqc
        assert fwdReads.mqcData["FastQC"]["General_Stats"]["total_sequences"] <= 2000
        assert revReads.mqcData["FastQC"]["General_Stats"]["total_sequences"] <= 2000

        # test plot data from multiqc
        # should return a dataframe representing plot data for the sample
        assert set(fwdReads.mqcData.keys()) == {"Cutadapt", "FastQC"}
        assert set(revReads.mqcData.keys()) == {"Cutadapt", "FastQC"}
        assert set(fwdReads.mqcData["FastQC"].keys()) == {"General_Stats", "Plots"}
        assert set(revReads.mqcData["FastQC"].keys()) == {"General_Stats", "Plots"}
        assert set(fwdReads.mqcData["Cutadapt"].keys()) == {"General_Stats", "Plots"}
        assert set(revReads.mqcData["Cutadapt"].keys()) == {"General_Stats", "Plots"}


def test_bulkRNASeq_STAGE01_Components_single(glds48_dataSystem_STAGE01):
    # iterate through all raw read components
    for sample in glds48_dataSystem_STAGE01.dataset.samples.values():

        reads = sample.trimReads

        # test general stats from multiqc
        assert reads.mqcData["FastQC"]["General_Stats"]["total_sequences"] <= 2000

        # test plot data from multiqc
        # should return a dataframe representing plot data for the sample
        assert set(reads.mqcData.keys()) == {"Cutadapt", "FastQC"}
        assert set(reads.mqcData["FastQC"].keys()) == {"General_Stats", "Plots"}
        assert set(reads.mqcData["Cutadapt"].keys()) == {"General_Stats", "Plots"}
