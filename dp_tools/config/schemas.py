""" Schemas for validation 
Uses Schema to allow usage of validation functions
"""
from schema import Schema, Optional

runsheet = {
    "bulkRNASeq": Schema(
        {
            str: {
                "Original Sample Name": str,
                "has_ERCC": bool,
                "organism": str,
                "paired_end": bool,
                "read1_path": str,
                Optional("read2_path"): str,
                str: object,  # this is used to pass other columns, chiefly Factor Value ones
            }
        }
    ),
    "microarray": Schema(
        {
            str: {
                "Original Sample Name": str,
                "organism": str,
                "Study Assay Measurement": str,
                "Study Assay Technology Type": str,
                "Study Assay Technology Platform": str,
                "Source Name": str,
                "Label": str,
                "Hybridization Assay Name": str,
                "Array Data File Name": str,
                "Array Data File Path": str,
                str: object,  # this is used to pass other columns, chiefly Factor Value ones
            }
        }
    ),
    "methylSeq": Schema(
        {
            str: {
                "Original Sample Name": str,
                "organism": str,
                "paired_end": bool,
                "read1_path": str,
                Optional("read2_path"): str,
                str: object,  # this is used to pass other columns, chiefly Factor Value ones
            }
        }
    ),
}
