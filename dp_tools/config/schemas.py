""" Schemas for validation 
Uses Schema to allow usage of validation functions
"""
from schema import Schema
from schema import Optional as schema_Optional
from typing import Optional
import pandera as pa

check_single_value = pa.Check(
    lambda x: len(x.unique()) == 1,
    title="Check that all values in the column are identical",
    description="Useful for columns that contain dataset level metadata like organism and paired_end.",
    error="Dataset level columns do NOT contain one unique value"
    )

check_read2_path_populated_if_paired_end = pa.Check(
    lambda df: ("read2_path" in df.columns and df['paired_end'][0] == True) or
               ("read2_path" not in df.columns and df['paired_end'][0] == False),
    title="Check 'read2_path' is only populated if paired_end is True",
    description="Failures here are likely either due to manual user error or inappropriate source file (e.g. ISA archive)",
    error="Expected 'read2_path' to be populated only if paired_end is True"
    )

runsheet = {
    "bulkRNASeq": pa.DataFrameSchema(
        columns={
            "Original Sample Name": pa.Column(str),
            "has_ERCC": pa.Column(bool, check_single_value),
            "organism": pa.Column(str, check_single_value),
            "paired_end": pa.Column(bool, check_single_value),
            "read1_path": pa.Column(str),
            "read2_path": pa.Column(str, required=False), # Expect if paired_end is True
        },
        # define checks at the DataFrameSchema-level
        checks=check_read2_path_populated_if_paired_end
    ),
    "methylSeq": pa.DataFrameSchema(
        columns={
            "Original Sample Name": pa.Column(str),
            "has_ERCC": pa.Column(bool, check_single_value),
            "organism": pa.Column(str, check_single_value),
            "paired_end": pa.Column(bool, check_single_value),
            "read1_path": pa.Column(str),
            "read2_path": pa.Column(str, required=False), # Expect if paired_end is True
        },
        # define checks at the DataFrameSchema-level
        checks=check_read2_path_populated_if_paired_end
    )
}
