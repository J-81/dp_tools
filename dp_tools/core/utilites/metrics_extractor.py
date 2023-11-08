# Enum for assay types
import ast
from enum import Enum
import json
import pathlib
from pathlib import Path
from dataclasses import dataclass
from typing import Literal
import numpy as np

import yaml
import pandas as pd
from loguru import logger

from dp_tools.scripts import convert
from dp_tools.core.utilites import multiqc_tools


class AssayType(Enum):
    bulkRNASeq = 1
    bulkRNASeq_VV = 2
    scRNASeq = 3
    spatialRNASeq = 4


@dataclass
class MultiQCTargetSection:
    targets: list[Path]
    section_name: str
    modules: list[str]
    jsonTarget: list[str] | Literal[False]


class MetricsExtractor:
    # Class Attributes
    metrics: pd.DataFrame = pd.DataFrame()
    mqc_metrics: pd.DataFrame = pd.DataFrame()
    samplewise_metrics: pd.DataFrame = pd.DataFrame()

    def __init__(self, targets: list[MultiQCTargetSection]):
        self.targets = targets

    # Ensure all column names are tuples
    @staticmethod
    def ensure_tuple(col_name, desired_length=4):
        if isinstance(col_name, tuple):
            # Pad with None if the tuple is smaller than desired_length
            return col_name + (None,) * (desired_length - len(col_name))
        else:
            # Create a tuple of desired_length with the col_name as the first element
            return (col_name,) + (None,) * (desired_length - 1)

    def extract_general_information(self, assay_type: AssayType, yaml_file: str):
        """This function parses data from a yaml file that is applicable on a dataset level.

        Examples include:
            OSD-#
            GLDS-#
            Sample Name
            Organism
            Tissue Type
            Library prep method (e.g. ribo-deplete (aka: totRNA) or polyA-enriched (aka: mRNA))
            % rRNA contamination
            PE or SE
            Stranded or Unstranded
            Library prep kit
            Data source (GeneLab generated, User-submitted, or Federated)
        """
        # Parse the yaml file
        with open(yaml_file) as file:
            data = yaml.load(file, Loader=yaml.FullLoader)

        EXPECTED_KEYS = [
            "OSD-#",
            "GLDS-#",
            "Sample Name",
            "Organism",
            "Tissue Type",
            "Library prep method",
            "% rRNA contamination",
            "PE or SE",
            "Stranded or Unstranded",
            "Library prep kit",
            "Data source",
        ]

        # Validate all keys and report all missing
        missing_keys: list[str] = list()
        for key in EXPECTED_KEYS:
            if key not in data:
                missing_keys.append(key)

        if missing_keys:
            raise ValueError(f"Missing keys: {missing_keys}")

    def extract_sections(self):
        def _extract_section_from_json(
            self, section_name: str, json_file: Path, module: str
        ):
            # Load json data
            with open(json_file) as file:
                data = json.load(file)

            # Note: Certain modules like RSeQC don't produce a general stats table
            # So we need to check if it exists before trying to extract it
            if data["report_general_stats_data"]:
                flat_data = multiqc_tools.get_reformated_source_dict(
                    data["report_general_stats_data"][
                        0
                    ]  # assumes only one module per multiQC json file
                )

                df_general_stats = pd.DataFrame(flat_data).T

                # Add a section as the first part of the column MultiIndex
                df_general_stats.columns = pd.MultiIndex.from_tuples(
                    [
                        (section_name, f"multiqc_{module}", "general_stats", col)
                        for col in df_general_stats.columns
                    ]
                )

                # Add section_name as last part of row MultiIndex
                df_general_stats.index = df_general_stats.index.set_names(
                    ["sample name", "sample subcomponent"]
                )

                df_general_stats.index = pd.MultiIndex.from_tuples(
                    list(
                        zip(
                            df_general_stats.index.get_level_values(
                                "sample name"
                            ).str.replace(
                                "-", "_"
                            ),  # Plots data performs this conversion so we match it here
                            df_general_stats.index.get_level_values(
                                "sample subcomponent"
                            ),
                            [section_name] * len(df_general_stats.index),
                        )
                    )
                ).set_names(["sample name", "sample subcomponent", "name"])

                df_updated_metrics = df_general_stats
                has_general_stats = True
            else:
                has_general_stats = False

            # Same for plot data
            df_plot_data = multiqc_tools.format_plots_as_dataframe(
                data["report_plot_data"]
            )
            # Add section_name as last part of row MultiIndex
            df_plot_data.index = df_plot_data.index.set_names(
                ["sample name", "sample subcomponent"]
            )

            df_plot_data.index = pd.MultiIndex.from_tuples(
                list(
                    zip(
                        df_plot_data.index.get_level_values("sample name"),
                        df_plot_data.index.get_level_values("sample subcomponent"),
                        [section_name] * len(df_plot_data.index),
                    )
                )
            ).set_names(["sample name", "sample subcomponent", "name"])

            # Add a section as the first part of the column MultiIndex
            df_plot_data.columns = pd.MultiIndex.from_tuples(
                [(section_name, *col) for col in df_plot_data.columns]
            )

            if has_general_stats:
                df_updated_metrics = df_updated_metrics.merge(
                    df_plot_data, left_index=True, right_index=True
                )
            else:
                df_updated_metrics = df_plot_data

            # Convert all columns to tuples
            columns_as_tuples = df_updated_metrics.columns.map(self.ensure_tuple)

            # Create MultiIndex
            df_updated_metrics.columns = pd.MultiIndex.from_tuples(columns_as_tuples)

            self.metrics = self.metrics.append(df_updated_metrics)

        def _extract_section(
            self, section_name: str, files: list[Path], modules: list[str]
        ):
            mqc_ret = multiqc_tools.get_parsed_data(
                input_f=[str(f) for f in files], modules=modules, as_dataframe=False
            )
            flat_data = multiqc_tools.flatten_raw_data(mqc_ret["report"])

            df_general_stats = pd.DataFrame(flat_data).T

            # Add section_name as last part of row MultiIndex
            df_general_stats.index = df_general_stats.index.set_names(
                ["sample name", "sample subcomponent"]
            )

            df_general_stats.index = pd.MultiIndex.from_tuples(
                list(
                    zip(
                        df_general_stats.index.get_level_values("sample name"),
                        df_general_stats.index.get_level_values("sample subcomponent"),
                        [section_name] * len(df_general_stats.index),
                    )
                )
            ).set_names(["sample name", "sample subcomponent", "name"])

            # Metrics names may include '-', whereas all multiQC names covert these to '_'
            # So here we create a temporary column with the '-' replaced with '_' for merging purposes
            # if isinstance(self.metrics.index, pd.MultiIndex):
            #     idx_sample_name = self.metrics.index.get_level_values('sample name')
            #     idx_sample_name = idx_sample_name.str.replace('-','_')
            #     self.metrics.index = pd.MultiIndex.from_arrays(
            #             [
            #             idx_sample_name,
            #             self.metrics.index.get_level_values('sample subcomponent'),
            #             self.metrics.index.get_level_values('name'),
            #             ],
            #             names = ['sample name','sample subcomponent','name']
            #         )
            # else:
            #     self.metrics.index = self.metrics.index.str.replace('-','_')
            # # self.metrics.index = self.metrics.index.str.replace('-','_')

            df_updated_metrics = df_general_stats

            # Same for plot data
            df_plot_data = multiqc_tools.format_plots_as_dataframe(mqc_ret)
            # Add section_name as last part of row MultiIndex
            df_plot_data.index = df_plot_data.index.set_names(
                ["sample name", "sample subcomponent"]
            )

            df_plot_data.index = pd.MultiIndex.from_tuples(
                list(
                    zip(
                        df_plot_data.index.get_level_values("sample name"),
                        df_plot_data.index.get_level_values("sample subcomponent"),
                        [section_name] * len(df_plot_data.index),
                    )
                )
            ).set_names(["sample name", "sample subcomponent", "name"])

            df_updated_metrics = df_updated_metrics.merge(
                df_plot_data, left_index=True, right_index=True
            )

            # Add a section as the first part of the column MultiIndex
            df_updated_metrics.columns = pd.MultiIndex.from_tuples(
                [(section_name, *col) for col in df_updated_metrics.columns]
            )

            self.metrics = self.metrics.append(df_updated_metrics)

        for target in self.targets:
            if target.jsonTarget == False:
                _extract_section(
                    self, target.section_name, target.targets, target.modules
                )
            else:
                _extract_section_from_json(
                    self, target.section_name, target.jsonTarget, target.modules[0]
                )  # Restriction: Can only handle one module multiQC json
        # Convert index of three part tuple to MultiIndex
        # Unnamed so must access part position in tuple
        # self.metrics.index = pd.MultiIndex.from_tuples(self.metrics.index, names = ['sample name','sample subcomponent','name'])
        # self.metrics.index = self.metrics.index.set_names(['sample name','sample subcomponent','name'])

        # Merge in samplewise metrics
        metrics_reset = self.metrics.reset_index(level=["sample subcomponent", "name"])

        samplewise_metrics_cleaned = self.samplewise_metrics.copy()
        samplewise_metrics_cleaned.index = samplewise_metrics_cleaned.index.str.replace(
            "-", "_"
        )

        merged = metrics_reset.merge(
            samplewise_metrics_cleaned,
            how="left",
            left_on="sample name",
            right_index=True,
        )
        # Rename based on length of coerced tuples
        merged = merged.rename(
            columns={
                ("sample subcomponent", "", "", ""): "sample subcomponent",
                ("name", "", "", ""): "name",
            }
        )
        merged = merged.set_index(["sample subcomponent", "name"], append=True)

        self.metrics = merged

    def extract_data_from_isa(
        self, accession: str, isa_archive: pathlib.Path, config: tuple[str, str]
    ):
        class mock_schema:
            @staticmethod
            def validate(df):
                pass

        samplewise_metrics = convert.isa_to_runsheet(
            accession,
            isa_archive,
            config,
            schema=mock_schema(),  # type: ignore
            assert_factor_values=False,
        )

        self.samplewise_metrics = samplewise_metrics

    def append_manual_yaml_data(self, target_yaml: Path):
        # Start with df_isa and add columns for each key value in yaml
        with open(target_yaml) as file:
            new_data = yaml.safe_load(file)

        # Add the new data to the existing data as new columns
        for key, value in new_data.items():
            self.samplewise_metrics[key] = value

    def process_metrics(self, assay_type: AssayType):
        match assay_type:
            case AssayType.bulkRNASeq:
                df_samplewise = pd.DataFrame()

                # Copy here is inefficient but useful to keep original dataframe unmodified
                df_interim = self.metrics.copy()

                # Ensure all column names are tuples
                def ensure_tuple(col_name, desired_length=4):
                    if isinstance(col_name, tuple):
                        # Pad with None if the tuple is smaller than desired_length
                        return col_name + (None,) * (desired_length - len(col_name))
                    else:
                        # Create a tuple of desired_length with the col_name as the first element
                        return (col_name,) + (None,) * (desired_length - 1)

                # Convert all columns to tuples
                columns_as_tuples = df_interim.columns.map(ensure_tuple)

                # Create MultiIndex
                df_interim.columns = pd.MultiIndex.from_tuples(columns_as_tuples)

                def _process_fastqc_data(df_full: pd.DataFrame, section_name: str):
                    df_samplewise = pd.DataFrame()
                    # Raw reads
                    df_fastqc_subset = df_full.xs(
                        key=(section_name, "_R1"),
                        axis="rows",
                        level=["name", "sample subcomponent"],
                    )

                    # M Seqs (read depth)
                    df_samplewise.index = df_fastqc_subset.index
                    try:
                        df_samplewise["Total Seqs"] = df_fastqc_subset[
                            (
                                section_name,
                                "multiqc_fastqc",
                                "general_stats",
                                "Total Sequences",
                            )
                        ].astype(int)
                    except (
                        KeyError
                    ):  # Sometimes the key is named differently (e.g. OSD-511)
                        df_samplewise["Total Seqs"] = df_fastqc_subset[
                            (
                                section_name,
                                "multiqc_fastqc",
                                "general_stats",
                                "total_sequences",
                            )
                        ].astype(int)

                    # Read length (This can be a range, especially for trimmed reads)
                    df_samplewise["Mean Read Length"] = df_fastqc_subset[
                        (
                            section_name,
                            "multiqc_fastqc",
                            "general_stats",
                            "avg_sequence_length",
                        )
                    ]

                    # # Check column type and perform action
                    # def _process_read_length_column(column):
                    #     if np.issubdtype(column.dtype, np.number):
                    #         return column.astype(int), column.astype(int)
                    #     elif column.dtype == object:
                    #         min_length = column.str.split("-").apply(
                    #             lambda read_length_list: read_length_list[0]
                    #         )
                    #         max_length = column.str.split("-").apply(
                    #             lambda read_length_list: read_length_list[1]
                    #         )
                    #         return min_length.astype(int), max_length.astype(int)
                    #     raise ValueError(
                    #         "Column type not recognized. Expected int or string like object (i.e. '20-151')"
                    #     )

                    # (
                    #     df_samplewise["Min Read Length"],
                    #     df_samplewise["Max Read Length"],
                    # ) = _process_read_length_column(read_length_series)

                    # Mean & Median Q Score (Across all bases)
                    df_samplewise = df_samplewise.merge(
                        (
                            df_full.xs(
                                key=section_name,
                                axis="rows",
                                level="name",
                            )
                            .xs(
                                key=(section_name, "FastQC: Mean Quality Scores"),
                                axis="columns",
                                level=[0, 1],
                            )
                            .agg(["mean", "median"], axis="columns")
                        ).rename(
                            columns={
                                "mean": "Average Q Score (Across all plotted base positions)",
                                "median": "Median Q Score (Across all plotted base positions)",
                            }
                        ),
                        left_index=True,
                        right_index=True,
                    )

                    # % Dups
                    try:
                        df_samplewise = 100 - df_samplewise.merge(
                            (
                                df_full.xs(
                                    key=(
                                        section_name,
                                        "multiqc_fastqc",
                                        "general_stats",
                                    ),
                                    axis="columns",
                                ).xs(key=(section_name), axis="rows", level="name")[
                                    "total_deduplicated_percentage"
                                ]
                            ).rename("% Dups"),
                            left_index=True,
                            right_index=True,
                        )
                    except:  # Sometimes key names are different
                        df_samplewise = df_samplewise.merge(
                            (
                                df_full.xs(
                                    key=(
                                        section_name,
                                        "multiqc_fastqc",
                                        "general_stats",
                                    ),
                                    axis="columns",
                                ).xs(key=(section_name), axis="rows", level="name")[
                                    "percent_duplicates"
                                ]
                            ).rename("% Dups"),
                            left_index=True,
                            right_index=True,
                        )

                    # Mean %GC
                    try:
                        df_samplewise = df_samplewise.merge(
                            df_full.xs(
                                key=(
                                    section_name,
                                    "multiqc_fastqc",
                                    "general_stats",
                                    "%GC",
                                ),
                                axis="columns",
                            )
                            .xs(key=(section_name), axis="rows", level="name")
                            .rename("Mean %GC"),
                            left_index=True,
                            right_index=True,
                        )
                    except KeyError:  # %GC <-> percent_gc
                        df_samplewise = df_samplewise.merge(
                            df_full.xs(
                                key=(
                                    section_name,
                                    "multiqc_fastqc",
                                    "general_stats",
                                    "percent_gc",
                                ),
                                axis="columns",
                            )
                            .xs(key=(section_name), axis="rows", level="name")
                            .rename("Mean %GC"),
                            left_index=True,
                            right_index=True,
                        )

                    def _first_col_reaching_min(row, min_value):
                        # Use boolean indexing to get columns that are >= min_value
                        valid = row[row >= min_value]

                        # Return the first column name that satisfies the condition, or None if not found
                        return valid.index[0] if not valid.empty else None

                    def _last_col_reaching_min(row, min_value):
                        # Use boolean indexing to get columns that are >= min_value
                        valid = row[row >= min_value]

                        # Return the last column name that satisfies the condition, or None if not found
                        return valid.index[-1] if not valid.empty else None

                    def _get_first_column_where_cumulative_sum_exceeds_proportion_of_row_sum(
                        row, proportion
                    ):
                        # Find all columns that exceed the proportion of the row sum
                        valid = row[row.cumsum() >= proportion * row.sum()]

                        # Return the first column name that satisfies the condition, or None if not found
                        return valid.index[0] if not valid.empty else None

                    df_gc_plots = df_full.xs(
                        key=section_name,
                        axis="rows",
                        level="name",
                    ).xs(
                        key=(
                            section_name,
                            "FastQC: Per Sequence GC Content",
                            "FastQC: Per Sequence GC Content",
                        ),
                        axis="columns",
                        level=[0, 1, 2],
                    )

                    # Min %GC Reaching %1 Counts
                    df_samplewise = df_samplewise.merge(
                        (
                            df_gc_plots.apply(
                                _first_col_reaching_min, axis="columns", min_value=1
                            ).apply(lambda s: s.split()[0])
                        ).rename("Min %GC reaching 1% Counts"),
                        left_index=True,
                        right_index=True,
                    )

                    # Max %GC Reaching %1 Counts
                    df_samplewise = df_samplewise.merge(
                        (
                            df_gc_plots.apply(
                                _last_col_reaching_min, axis="columns", min_value=1
                            ).apply(lambda s: s.split()[0])
                        ).rename("Max %GC reaching 1% Counts"),
                        left_index=True,
                        right_index=True,
                    )

                    # 25% Quartile AUC Point %GC
                    df_samplewise = df_samplewise.merge(
                        (
                            df_gc_plots.apply(
                                _get_first_column_where_cumulative_sum_exceeds_proportion_of_row_sum,
                                axis="columns",
                                proportion=0.25,
                            ).apply(lambda s: s.split()[0])
                        ).rename("25% Quartile AUC Point %GC"),
                        left_index=True,
                        right_index=True,
                    )
                    # 50% Quartile AUC Point %GC
                    df_samplewise = df_samplewise.merge(
                        (
                            df_gc_plots.apply(
                                _get_first_column_where_cumulative_sum_exceeds_proportion_of_row_sum,
                                axis="columns",
                                proportion=0.50,
                            ).apply(lambda s: s.split()[0])
                        ).rename("50% Quartile AUC Point %GC"),
                        left_index=True,
                        right_index=True,
                    )
                    # 75% Quartile AUC Point %GC
                    df_samplewise = df_samplewise.merge(
                        (
                            df_gc_plots.apply(
                                _get_first_column_where_cumulative_sum_exceeds_proportion_of_row_sum,
                                axis="columns",
                                proportion=0.75,
                            ).apply(lambda s: s.split()[0])
                        ).rename("75% Quartile AUC Point %GC"),
                        left_index=True,
                        right_index=True,
                    )

                    # % N Content
                    df_n_content_plots = df_full.xs(
                        key=section_name,
                        axis="rows",
                        level="name",
                    ).xs(
                        key=(
                            section_name,
                            "FastQC: Per Base N Content",
                            "FastQC: Per Base N Content",
                        ),
                        axis="columns",
                        level=[0, 1, 2],
                    )

                    # % N Content Summed Across All Plotted Bases Positions
                    df_samplewise = df_samplewise.merge(
                        df_n_content_plots.sum(axis="columns").rename(
                            "% N Content Summed Across All Plotted Bases Positions"
                        ),
                        left_index=True,
                        right_index=True,
                    )

                    return df_samplewise

                def _process_align_data(df_full: pd.DataFrame, section_name: str):
                    df_samplewise = pd.DataFrame()

                    df_align_subset = (
                        df_full.xs(
                            key=section_name,
                            axis="rows",
                            level="name",
                        )
                        .xs(key=section_name, axis="columns", level=0)
                        .droplevel("sample subcomponent", axis="rows")
                    )

                    df_samplewise.index = df_align_subset.index
                    df_samplewise["% Uniquely mapped"] = df_align_subset[
                        (
                            "multiqc_star",
                            "general_stats",
                            "uniquely_mapped_percent",
                        )
                    ].astype(float)

                    df_samplewise["% Mapped to multiple loci"] = df_align_subset[
                        ("multiqc_star", "general_stats", "multimapped_percent")
                    ].astype(float)

                    df_samplewise["% Mapped to too many loci"] = df_align_subset[
                        (
                            "multiqc_star",
                            "general_stats",
                            "multimapped_toomany_percent",
                        )
                    ].astype(float)

                    df_samplewise["% Unmapped too short"] = df_align_subset[
                        (
                            "multiqc_star",
                            "general_stats",
                            "unmapped_tooshort_percent",
                        )
                    ].astype(float)

                    df_samplewise["% Unmapped other"] = df_align_subset[
                        (
                            "multiqc_star",
                            "general_stats",
                            "unmapped_other",
                        )
                    ].astype(float)

                    return df_samplewise

                def _process_rseqc_genebody_coverage_data(
                    df_full: pd.DataFrame, section_name: str
                ):
                    df_samplewise = pd.DataFrame()

                    df_rseqc_subset = (
                        df_full.xs(
                            key=section_name,
                            axis="rows",
                            level="name",
                        )
                        .xs(key=section_name, axis="columns", level=0)
                        .droplevel("sample subcomponent", axis="rows")
                    )

                    df_samplewise.index = df_rseqc_subset.index

                    # Average % Coverage from 5-20 percentile (5' end coverage)
                    def _get_mean_for_percentile_range(min_range: int, max_range: int):
                        level_0_value = "RSeQC: Gene Body Coverage"
                        level_1_value = "RSeQC: Gene Body Coverage"
                        level_2_values_to_select = [
                            f"{i} Gene Body Percentile (5' -> 3') (% Coverage)"
                            for i in range(min_range, max_range + 1)
                        ]  # list of desired values for third level

                        mask = (
                            (
                                df_rseqc_subset.columns.get_level_values(0)
                                == level_0_value
                            )
                            & (
                                df_rseqc_subset.columns.get_level_values(1)
                                == level_1_value
                            )
                            & (
                                df_rseqc_subset.columns.get_level_values(2).isin(
                                    level_2_values_to_select
                                )
                            )
                        )

                        return (
                            df_rseqc_subset.loc[:, mask]
                            .astype(float)
                            .mean(axis="columns")
                        )

                    @dataclass
                    class TARGET_RANGE:
                        lower_bound: int
                        upper_bound: int
                        label: str

                    TARGET_RANGES: list[TARGET_RANGE] = [
                        TARGET_RANGE(5, 20, "5' end coverage"),
                        TARGET_RANGE(40, 60, "middle coverage"),
                        TARGET_RANGE(80, 95, "3' end coverage"),
                    ]

                    for target in TARGET_RANGES:
                        df_samplewise[
                            f"Average % Coverage from {target.lower_bound}-{target.upper_bound} percentile ({target.label})"
                        ] = _get_mean_for_percentile_range(
                            target.lower_bound, target.upper_bound
                        )

                    df_samplewise["Ratio of 3' end coverage to 5' end coverage"] = (
                        df_samplewise[
                            f"Average % Coverage from 80-95 percentile (3' end coverage)"
                        ]
                        / df_samplewise[
                            f"Average % Coverage from 5-20 percentile (5' end coverage)"
                        ]
                    )

                    return df_samplewise

                def _process_rseqc_infer_experiment_data(
                    df_full: pd.DataFrame, section_name: str
                ):
                    df_samplewise = pd.DataFrame()

                    df_rseqc_subset = (
                        df_full.xs(
                            key=section_name,
                            axis="rows",
                            level="name",
                        )
                        .xs(key=section_name, axis="columns", level=0)
                        .droplevel("sample subcomponent", axis="rows")
                    )

                    df_samplewise["% Sense"] = df_rseqc_subset[
                        (
                            "RSeQC: Infer experiment",
                            "RSeQC: Infer experiment",
                            "Sense (% Tags)",
                        )
                    ].astype(float)

                    df_samplewise["% Antisense"] = df_rseqc_subset[
                        (
                            "RSeQC: Infer experiment",
                            "RSeQC: Infer experiment",
                            "Antisense (% Tags)",
                        )
                    ].astype(float)

                    df_samplewise["% Undetermined"] = df_rseqc_subset[
                        (
                            "RSeQC: Infer experiment",
                            "RSeQC: Infer experiment",
                            "Undetermined (% Tags)",
                        )
                    ].astype(float)

                    return df_samplewise

                def _process_rseqc_inner_distance_data(
                    df_full: pd.DataFrame, section_name: str
                ):
                    df_samplewise = pd.DataFrame()

                    df_rseqc_subset = (
                        df_full.xs(
                            key=section_name,
                            axis="rows",
                            level="name",
                        )
                        .xs(key=section_name, axis="columns", level=0)
                        .droplevel("sample subcomponent", axis="rows")
                    )

                    # Inner Distance Peak Distance
                    # Extract from column tuple
                    # Example: ('RSeQC: Inner Distance', 'RSeQC: Inner Distance', '-117.5 Inner Distance (bp) (Counts)')
                    # Yields: -117.5
                    try:
                        df_samplewise["Peak Inner Distance"] = (
                            df_rseqc_subset.idxmax(axis="columns")
                            .apply(lambda col: col[2])
                            .astype(float)
                        )
                    except (
                        ValueError
                    ):  # e.g. ValueError: could not convert string to float: '-142.5 Inner Distance (bp) (Counts)'
                        df_samplewise["Peak Inner Distance"] = (
                            df_rseqc_subset.idxmax(axis="columns")
                            .apply(lambda col: col[2].split()[0])
                            .astype(float)
                        )

                    # % Reads At Inner Distance Peak Distance
                    df_samplewise["% Reads At Peak Inner Distance"] = (
                        df_rseqc_subset.max(axis="columns")
                        / df_rseqc_subset.sum(axis="columns")
                        * 100
                    )

                    # TAGUP: Inner distance at 1%

                    return df_samplewise

                def _process_rseqc_read_distribution_data(
                    df_full: pd.DataFrame, section_name: str
                ):
                    df_samplewise = pd.DataFrame()

                    df_rseqc_subset = (
                        df_full.xs(
                            key=section_name,
                            axis="rows",
                            level="name",
                        )
                        .xs(key=section_name, axis="columns", level=0)
                        .droplevel("sample subcomponent", axis="rows")
                    )

                    @dataclass
                    class TARGET_LABELS:
                        dataframe_name: str
                        metrics_name: str

                    TARGETS: list[TARGET_LABELS] = [
                        TARGET_LABELS("CDS_Exons (# Tags)", "% CDS_Exons"),
                        TARGET_LABELS("5'UTR_Exons (# Tags)", "% 5'UTR_Exons"),
                        TARGET_LABELS("3'UTR_Exons (# Tags)", "% 3'UTR_Exons"),
                        TARGET_LABELS("Introns (# Tags)", "% Introns"),
                        TARGET_LABELS("TSS_up_1kb (# Tags)", "% TSS_up_1kb"),
                        TARGET_LABELS("TSS_up_1kb-5kb (# Tags)", "% TSS_up_1kb-5kb"),
                        TARGET_LABELS("TSS_up_5kb-10kb (# Tags)", "% TSS_up_5kb-10kb"),
                        TARGET_LABELS("TES_down_1kb (# Tags)", "% TES_down_1kb"),
                        TARGET_LABELS("TES_down_1kb-50kb (# Tags)", "% TES_down_1kb-5kb"),
                        TARGET_LABELS("TES_down_5kb-10kb", "% TES_down_5kb-10kb"),
                        TARGET_LABELS(
                            "Other_intergenic (# Tags)", "% Other_intergenic"
                        ),
                    ]

                    for target in TARGETS:
                        try:
                            # Plot data
                            df_samplewise[target.metrics_name] = df_rseqc_subset[
                                (
                                    "RSeQC: Read Distribution",
                                    "RSeQC: Read Distribution",
                                    target.dataframe_name,
                                )
                            ]
                        except:
                            # No plot data means zero for the given tag
                            df_samplewise[target.metrics_name] = 0

                    # Convert all to percents by summing across row and dividing each by sum
                    df_samplewise = df_samplewise.apply(
                        lambda col: col / df_samplewise.sum(axis="columns") * 100
                    )

                    return df_samplewise

                def _process_rsem_data(df_full: pd.DataFrame, section_name: str):
                    df_samplewise = pd.DataFrame()
                    print(list(df_full.columns))
                    df_align_subset = (
                        df_full.xs(
                            key=section_name,
                            axis="rows",
                            level="name",
                        )
                        .xs(key=section_name, axis="columns", level=0)
                        .droplevel("sample subcomponent", axis="rows")
                    )

                    df_samplewise.index = df_align_subset.index
                    df_samplewise["% Aligned uniquely to a gene"] = df_align_subset[
                        ("multiqc_rsem", "general_stats", "Unique")
                    ].astype(float)

                    df_samplewise["% Aligned to multiple genes"] = df_align_subset[
                        ("multiqc_rsem", "general_stats", "Multi")
                    ].astype(float)

                    df_samplewise[
                        "% Filtered due to too many alignments"
                    ] = df_align_subset[
                        (
                            "multiqc_rsem",
                            "general_stats",
                            "Filtered",
                        )
                    ].astype(
                        float
                    )

                    df_samplewise["% Unalignable reads"] = df_align_subset[
                        (
                            "multiqc_rsem",
                            "general_stats",
                            "Unalignable",
                        )
                    ].astype(float)

                    # At this point, all in counts
                    # Convert all to percents by summing across row and dividing each by sum
                    df_samplewise = df_samplewise.apply(
                        lambda col: col / df_align_subset[
                        (
                            "multiqc_rsem",
                            "general_stats",
                            "Total",
                        )
                    ].astype(float) * 100
                    )

                    return df_samplewise
                
                df_samplewise_raw = _process_fastqc_data(df_interim, "raw reads")
                df_samplewise_trimmed = _process_fastqc_data(
                    df_interim, "trimmed reads"
                )

                df_samplewise_align = _process_align_data(df_interim, "aligned reads")

                df_samplewise_rseqc_genebody_coverage = (
                    _process_rseqc_genebody_coverage_data(
                        df_interim, "rseqc: genebody coverage"
                    )
                )

                df_samplewise_rseqc_infer_experiment = (
                    _process_rseqc_infer_experiment_data(
                        df_interim, "rseqc: infer experiment"
                    )
                )

                df_samplewise_rseqc_inner_distance = _process_rseqc_inner_distance_data(
                    df_interim, "rseqc: inner distance"
                )

                df_samplewise_rseqc_read_distribution = (
                    _process_rseqc_read_distribution_data(
                        df_interim, "rseqc: read distribution"
                    )
                )

                df_samplewise_rsem = _process_rsem_data(df_interim, "rsem count")

                # Merge all
                df_merged = (
                    df_samplewise_raw.merge(
                        df_samplewise_trimmed,
                        left_index=True,
                        right_index=True,
                        suffixes=(" Raw", " Trimmed"),
                    )
                    .merge(
                        df_samplewise_align,
                        left_index=True,
                        right_index=True,
                    )
                    .merge(
                        df_samplewise_rseqc_genebody_coverage,
                        left_index=True,
                        right_index=True,
                    )
                    .merge(
                        df_samplewise_rseqc_infer_experiment,
                        left_index=True,
                        right_index=True,
                    )
                    .merge(
                        df_samplewise_rseqc_inner_distance,
                        left_index=True,
                        right_index=True,
                    )
                    .merge(
                        df_samplewise_rseqc_read_distribution,
                        left_index=True,
                        right_index=True,
                    )
                    .merge(
                        df_samplewise_rsem,
                        left_index=True,
                        right_index=True,
                    )
                    .sort_index()
                )

                return df_merged

            case AssayType.bulkRNASeq_VV:
                df_samplewise = pd.DataFrame()

                # Copy here is inefficient but useful to keep original dataframe unmodified
                df_interim = self.metrics.copy()

                # Convert all columns to tuples
                columns_as_tuples = df_interim.columns.map(self.ensure_tuple)

                # Create MultiIndex
                df_interim.columns = pd.MultiIndex.from_tuples(columns_as_tuples)

                # Raw reads
                raw_reads = df_interim.xs(
                    key=("raw reads", "_R1"),
                    axis="rows",
                    level=["name", "sample subcomponent"],
                )

                # Read Depth Range
                df_samplewise.index = raw_reads.index
                df_samplewise["Total Seqs"] = raw_reads[
                    ("raw reads", "multiqc_fastqc", "general_stats", "Total Sequences")
                ].astype(int)

                # Read length
                df_samplewise["Read Length"] = raw_reads[
                    ("raw reads", "multiqc_fastqc", "general_stats", "Sequence length")
                ].astype(int)

                # Mean & Median Q Score (Across all bases)
                df_samplewise = df_samplewise.merge(
                    (
                        df_interim.xs(
                            key="raw reads",
                            axis="rows",
                            level="name",
                        )
                        .xs(
                            key=("raw reads", "FastQC: Mean Quality Scores"),
                            axis="columns",
                            level=[0, 1],
                        )
                        .agg(["mean", "median"], axis="columns")
                    ).rename(
                        columns={
                            "mean": "Average Q Score (Across all plotted base positions)",
                            "median": "Median Q Score (Across all plotted base positions)",
                        }
                    ),
                    left_index=True,
                    right_index=True,
                )

                # % Dups

                print("DONE")
            case _:
                raise NotImplementedError(
                    f"Assay type {assay_type} not implemented for summarization."
                )

    def load_metrics_csv(self, metrics_csv: Path):
        # check\ metrics hasn't been created yet or loaded
        assert self.metrics.equals(
            pd.DataFrame()
        ), "Metrics already loaded. Please create a new MetricsExtractor object."

        self.metrics = pd.read_csv(metrics_csv, index_col=[0, 1, 2])

        # Set index names
        self.metrics.index = self.metrics.index.set_names(
            ["sample name", "sample subcomponent", "name"]
        )

        # Convert column names to tuples if they represent valid tuples
        def _convert_to_tuple_if_valid(col_name):
            try:
                # Check if the column name can be evaluated to a tuple
                result = ast.literal_eval(col_name)
                if isinstance(result, tuple):
                    return result
            except (SyntaxError, ValueError):
                pass
            return col_name

        self.metrics.columns = [
            _convert_to_tuple_if_valid(col) for col in self.metrics.columns
        ]


def generate_extractor_from_yaml_config(config: Path) -> MetricsExtractor:
    with open(config) as file:
        config_data = yaml.safe_load(file)

    targets: list[MultiQCTargetSection] = list()

    for section in config_data["Extraction Settings"]["sections"]:
        if not section["enabled"]:
            logger.info(f"Skipping {section['name']} because it is disabled.")
            continue

        # Set up MultiQC targets
        search_dir = Path(
            config_data["Extraction Settings"]["root search directory"]
        ) / Path(*section["multiQC"]["logs directory"])

        if section["multiQC"].get("from json", False):
            jsonTarget = Path(
                config_data["Extraction Settings"]["root search directory"]
            ) / Path(*section["multiQC"]["from json"])
        else:
            jsonTarget = False

        found_files: list[Path] = list()
        for logs_pattern in section["multiQC"]["logs pattern(s)"]:
            if section["multiQC"]["search recursively"]:
                found_files.extend(list(search_dir.rglob(logs_pattern)))
            else:
                found_files.extend(list(search_dir.glob(logs_pattern)))

        # Catch empty lists
        if len(found_files) == 0 and not jsonTarget:
            raise ValueError(
                f"No files found for {section['name']}. Configuration may be broken or consider disabling section if data is not present."
            )

        targets.append(
            MultiQCTargetSection(
                targets=found_files,
                section_name=section["name"],
                modules=section["multiQC"]["modules"],
                jsonTarget=jsonTarget,
            )
        )

    return MetricsExtractor(targets=targets)
