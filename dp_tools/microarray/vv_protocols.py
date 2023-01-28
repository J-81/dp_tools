from pathlib import Path
import re
from typing import Union
import yaml
import logging

from dp_tools.core.entity_model import Dataset

log = logging.getLogger(__name__)

from dp_tools.core.check_model import ValidationProtocol
# ORDER is intentional to allow shared name functions to overload in favor of microarray
from dp_tools.bulkRNASeq.checks import *  # normally this isn't ideal, however, as a library of functions this seems reasonable
from dp_tools.microarray.checks import *  # normally this isn't ideal, however, as a library of functions this seems reasonable



CONFIG = {
    "Metadata-check_metadata_attributes_exist": {
        "expected_attrs": ["biomart_attribute", "organism"]
    },
}


def validate_Agilent1Channel(
    dataset: Dataset,
    config_path: Path = None,
    run_args: dict = None,
    report_args: dict = None,
    protocol_args: dict = None,
    defer_run: bool = False,
) -> Union[ValidationProtocol, ValidationProtocol.Report]:

    if config_path is not None:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
    else:
        config = CONFIG

    if run_args is None:
        run_args = dict()

    if report_args is None:
        report_args = dict()

    if protocol_args is None:
        protocol_args = dict()
    # init validation protocol
    vp = ValidationProtocol(**protocol_args)
    # fmt: on
    with vp.component_start(
        name=dataset.name,
        description="Validate microarray processed data",
    ):
        with vp.component_start(
            name="Metadata", description="Metadata file validation"
        ):
            with vp.payload(payloads=[{"dataset": dataset}]):
                vp.add(
                    check_metadata_attributes_exist,
                    config=config["Metadata-check_metadata_attributes_exist"],
                )
        with vp.component_start(
            name="DGE Metadata",
            description="",
        ):

            with vp.component_start(
                name="Sample Table",
                description="",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "runsheet": lambda: dataset.data_assets["runsheet"].path,
                            "sampleTable": lambda: dataset.data_assets[
                                "sample table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(
                        check_sample_table_against_runsheet,
                        config={"all_samples_required": True},
                    )
                    vp.add(check_sample_table_for_correct_group_assignments)

            with vp.component_start(
                name="Contrasts Tables",
                description="",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "runsheet": lambda: dataset.data_assets["runsheet"].path,
                            "contrasts_table": lambda: dataset.data_assets[
                                "DE contrasts table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(check_contrasts_table_headers)
                    vp.add(check_contrasts_table_rows)
        with vp.component_start(
            name="DGE Output",
            description="",
        ):
            with vp.payload(
                payloads=[
                    {
                        "organism": lambda: dataset.metadata["organism"],
                        "samples": lambda: set(dataset.samples),
                        "dge_table": lambda: dataset.data_assets[
                            "DE annotated table"
                        ].path,
                        "runsheet": lambda: dataset.data_assets["runsheet"].path,
                    }
                ]
            ):
                vp.add(check_dge_table_annotation_columns_exist)
                vp.add(check_dge_table_sample_columns_exist)
                vp.add(check_dge_table_sample_columns_constraints)
                vp.add(check_dge_table_group_columns_exist)
                vp.add(check_dge_table_group_columns_constraints)
                vp.add(check_dge_table_comparison_statistical_columns_exist)
                vp.add(check_dge_table_group_statistical_columns_constraints)
                vp.add(check_dge_table_fixed_statistical_columns_exist)
                vp.add(check_dge_table_fixed_statistical_columns_constraints)
                vp.add(check_dge_table_log2fc_within_reason)

            with vp.component_start(
                name="Viz Tables",
                description="Extended from the dge tables",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "organism": lambda: dataset.metadata["organism"],
                            "samples": lambda: set(dataset.samples),
                            "dge_table": lambda: dataset.data_assets[
                                "DE annotated extended for viz table"
                            ].path,
                            "runsheet": lambda: dataset.data_assets["runsheet"].path,
                        }
                    ]
                ):
                    vp.add(check_dge_table_annotation_columns_exist)
                    vp.add(check_dge_table_sample_columns_exist)
                    vp.add(check_dge_table_sample_columns_constraints)
                    vp.add(check_dge_table_group_columns_exist)
                    vp.add(check_dge_table_group_columns_constraints)
                    vp.add(check_dge_table_comparison_statistical_columns_exist)
                    vp.add(check_dge_table_group_statistical_columns_constraints)
                    vp.add(check_dge_table_fixed_statistical_columns_exist)
                    vp.add(check_dge_table_fixed_statistical_columns_constraints)
                    vp.add(check_dge_table_log2fc_within_reason)
                    vp.add(check_viz_table_columns_exist)
                    vp.add(check_viz_table_columns_constraints)

                with vp.payload(
                    payloads=[
                        {
                            "samples": lambda: set(dataset.samples),
                            "pca_table": lambda: dataset.data_assets[
                                "viz PCA table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(check_viz_pca_table_index_and_columns_exist)
    # return protocol object without running or generating a report
    if defer_run:
        return vp

    vp.run(**run_args)

    # return report
    return vp.report(**report_args, combine_with_flags=dataset.loaded_assets_dicts)
