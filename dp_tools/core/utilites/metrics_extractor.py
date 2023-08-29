# Enum for assay types
from enum import Enum
import pathlib
from pathlib import Path
from dataclasses import dataclass

import yaml
import pandas as pd
from loguru import logger

from dp_tools.scripts import convert
from dp_tools.core.utilites import multiqc_tools

class AssayType(Enum):
    bulkRNASeq = 1
    scRNASeq = 2
    spatialRNASeq = 3

@dataclass
class MultiQCTargetSection():
    targets: list[Path]
    section_name: str
    modules: list[str]

class MetricsExtractor():

    # Class Attributes
    metrics: pd.DataFrame = pd.DataFrame()
    mqc_metrics: pd.DataFrame = pd.DataFrame()
    samplewise_metrics: pd.DataFrame = pd.DataFrame()

    def __init__(self,
                 targets: list[MultiQCTargetSection]
                 ):
        self.targets = targets

    def extract_general_information(self, assay_type: AssayType, yaml_file: str):
        """ This function parses data from a yaml file that is applicable on a dataset level.

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
            'OSD-#',
            'GLDS-#',
            'Sample Name',
            'Organism',
            'Tissue Type',
            'Library prep method',
            '% rRNA contamination',
            'PE or SE',
            'Stranded or Unstranded',
            'Library prep kit',
            'Data source'
        ]

        # Validate all keys and report all missing
        missing_keys: list[str] = list()
        for key in EXPECTED_KEYS:
            if key not in data:
                missing_keys.append(key)
        
        if missing_keys:
            raise ValueError(f"Missing keys: {missing_keys}")


    def extract_sections(self):

        def _extract_section(self, section_name: str, files: list[Path], modules: list[str]):
                mqc_ret = multiqc_tools.get_parsed_data(input_f = 
                    [str(f) for f in files], 
                    modules = modules,
                    as_dataframe = False)
                flat_data = multiqc_tools.flatten_raw_data(mqc_ret['report'])

                df_general_stats = pd.DataFrame(flat_data).T

                # Add section_name as last part of row MultiIndex
                df_general_stats.index = df_general_stats.index.set_names(['sample name','sample subcomponent'])
                
                df_general_stats.index = pd.MultiIndex.from_tuples(
                    list(zip(df_general_stats.index.get_level_values('sample name'),
                            df_general_stats.index.get_level_values('sample subcomponent'),
                            [section_name] * len(df_general_stats.index)))
                ).set_names(
                    ['sample name','sample subcomponent','name']
                )


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
                df_plot_data.index = df_plot_data.index.set_names(['sample name','sample subcomponent'])
                
                df_plot_data.index = pd.MultiIndex.from_tuples(
                    list(zip(df_plot_data.index.get_level_values('sample name'),
                            df_plot_data.index.get_level_values('sample subcomponent'),
                            [section_name] * len(df_plot_data.index)))
                ).set_names(
                   ['sample name','sample subcomponent','name']
                )
                
                df_updated_metrics = df_updated_metrics.merge(df_plot_data, left_index=True, right_index=True)

                # Add a section as the first part of the column MultiIndex
                df_updated_metrics.columns = pd.MultiIndex.from_tuples([(section_name, *col) for col in df_updated_metrics.columns])

                self.metrics = self.metrics.append(df_updated_metrics)

        for target in self.targets:
            _extract_section(self, target.section_name, target.targets, target.modules)
        # Convert index of three part tuple to MultiIndex
        # Unnamed so must access part position in tuple
        # self.metrics.index = pd.MultiIndex.from_tuples(self.metrics.index, names = ['sample name','sample subcomponent','name'])
        # self.metrics.index = self.metrics.index.set_names(['sample name','sample subcomponent','name'])

        # Merge in samplewise metrics
        metrics_reset = self.metrics.reset_index(level=['sample subcomponent','name'])
        
        samplewise_metrics_cleaned = self.samplewise_metrics.copy()
        samplewise_metrics_cleaned.index = samplewise_metrics_cleaned.index.str.replace('-','_')
        
        merged = metrics_reset.merge(samplewise_metrics_cleaned, how='left', left_on='sample name', right_index=True)
        # Rename based on length of coerced tuples
        merged = merged.rename(
            columns = {
                ('sample subcomponent','','',''): 'sample subcomponent',
                ('name','','',''): 'name',

            }
        )
        merged = merged.set_index(
            [
                'sample subcomponent',
                'name'
            ], append=True)
        
        self.metrics = merged

    def extract_data_from_isa(self, accession: str, isa_archive: pathlib.Path, config: tuple[str, str]):

        class mock_schema():
            @staticmethod
            def validate(df):
                pass

        samplewise_metrics = convert.isa_to_runsheet(accession, 
                                     isa_archive, 
                                     config, 
                                     schema=mock_schema(), # type: ignore
                                     assert_factor_values=False
                                     )
        
        self.samplewise_metrics = samplewise_metrics

    def append_manual_yaml_data(self, target_yaml: Path):
        # Start with df_isa and add columns for each key value in yaml
        with open(target_yaml) as file:
            new_data = yaml.safe_load(file)

        # Add the new data to the existing data as new columns
        for key, value in new_data.items():
            self.samplewise_metrics[key] = value
        
def generate_extractor_from_yaml_config(config: Path) -> MetricsExtractor:

    with open(config) as file:
        config_data = yaml.safe_load(file)

    targets: list[MultiQCTargetSection] = list()

    for section in config_data['Extraction Settings']['sections']:
        if not section['enabled']:
            logger.info(f"Skipping {section['name']} because it is disabled.")
            continue

        # Set up MultiQC targets
        search_dir = Path(config_data['Extraction Settings']['root search directory']) / Path(*section['multiQC']['logs directory'])
        found_files: list[Path] = list()
        for logs_pattern in section['multiQC']['logs pattern(s)']:
            if section['multiQC']['search recursively']:
                found_files.extend(list(search_dir.rglob(logs_pattern)))
            else:
                found_files.extend(list(search_dir.glob(logs_pattern)))
        
        # Catch empty lists
        if len(found_files) == 0:
            raise ValueError(f"No files found for {section['name']}. Configuration may be broken or consider disabling section if data is not present.")

        targets.append(
            MultiQCTargetSection(
                targets = found_files,
                section_name = section['name'],
                modules = section['multiQC']['modules']
            )
        )

    return MetricsExtractor(
        targets = targets
    )