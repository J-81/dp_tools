# function for running multiqc and getting data objects back
from collections import defaultdict
from pathlib import Path
import logging
from types import ModuleType
from typing import List, TypedDict

import multiqc
import pandas as pd

# iterable to remove suffixes and add them as subsource descriptors
SUBSOURCES = [
    "_R1_raw",
    "_R2_raw",
    "_R1",
    "_R2",
    "__STARpass1",
]

# iterable to remove suffixes that does NOT add them as subsource descriptors (often due to the name being redundantly associated with columns)
SCRUB_SAMPLES = ["_read_dist", "_infer_expt"]


def clean_messy_sample(messy_sample: str):
    # if any subsource suffixes found, cleans those from the sample name
    # returns two strings: cleaned_sample_name, suffix
    # suffix often denotes the 'sub source', e.g. forward read vs reverse read and should be included in downstream labeling
    current_sample = messy_sample
    sub_source = ""
    sub_source_suffix = ""  # used to add to sub source

    # first remove adaptor suffixes, harder to catch but hooking onto separater seems fine
    # E.g. 'Atha_Ler_0_sl_FLT_uG_Rep4_R2_raw _ Adapter 1' should move '_R2_raw:Adapter 1' into subsource
    adaptor_tokens = current_sample.split(" - ")
    if len(adaptor_tokens) == 1:  # 1 would indicate no such adaptors
        sub_source_suffix = ""
    elif (
        len(adaptor_tokens) == 2
    ):  # 'Atha_Ler_0_sl_FLT_uG_Rep4_R2_raw - Adapter 1' should trigger this
        logging.debug(f"Removing adaptor token to clean sample name")
        current_sample = adaptor_tokens[0]
        sub_source_suffix = adaptor_tokens[1]
    else:
        raise ValueError(
            f"Can't parse potential adaptors (using ' - ' split) from sample: {current_sample}"
        )

    # remove any substrings that don't indicate sub sources
    for scrub_s in SCRUB_SAMPLES:
        if current_sample.endswith(scrub_s):
            logging.debug(
                f"Removing {scrub_s} suffix to clean sample name, removing completely (redudant)"
            )
            current_sample = current_sample[: -len(scrub_s)]

    # split any true sub sources (e.g. _R1, _R2, _STARpass1)
    for subsource in SUBSOURCES:
        if current_sample.endswith(subsource):
            logging.debug(
                f"Removing {subsource} suffix to clean sample name, using to indicate sub source"
            )
            current_sample = current_sample[: -len(subsource)]
            sub_source += subsource
    # add adaptor suffix, no longer performed as the adapter string is added to the column multi-index
    # sub_source = (
    #    f"{sub_source}:{sub_source_suffix}" if sub_source_suffix else sub_source
    # )
    logging.info(
        f"Cleaned {messy_sample} to derive sample: {current_sample} and subsource: {sub_source}"
    )
    return current_sample, sub_source


def get_reformated_source_dict(source_dict: dict):
    logging.info(
        "Cleaning general data statistics dict and reformating to clean-sample-wise"
    )
    clean_source_dict = defaultdict(dict)
    for messy_sample_k, messy_sample_data in source_dict.items():
        clean_sample, sub_source = clean_messy_sample(messy_sample_k)
        for k, v in messy_sample_data.items():
            # new_k = f"{sub_source}::{k}" if sub_source else k
            entity_key = (clean_sample, sub_source)
            clean_source_dict[entity_key][k] = v
    return clean_source_dict


def get_parsed_data(
    input_f: List[str], modules: List[str] = [], as_dataframe: bool = True
):
    logging.info(f"Using MQC to parse: {input_f}")
    try:
        # a workaround for flushing handlers in MQC version 1.11
        logger = logging.getLogger("multiqc")
        [logger.removeHandler(h) for h in logger.handlers]
        mqc_ret = multiqc.run(
            input_f, no_data_dir=True, module=modules
        )  # note: empty list for modules falls back on all modules
        logging.info(f"Successfully parsed: {input_f}")
    except SystemExit as e:
        # a zero exit code indicates no data extracted (expected if the file isn't covered by a multiqc module yet)
        # catch non-zero exit codes, as raised when multiqc has an actual exception.
        assert (
            e.code == 0
        ), "MQC Module error occured: This is very rare and might indicate an edge case unhandled by Multiqc"
        logging.info(f"MQC does not currently parse: {input_f}, returning None")
        # other
        return None

    if not as_dataframe:
        return mqc_ret
    else:
        mqc_df = format_as_dataframe(mqc_ret["report"])
        return mqc_df


def flatten_raw_data(mqc_rep: dict) -> dict:
    """Generates a raw data flat dictionary from the full nested mqc report

    :param mqc_rep: [description]
    :type mqc_rep: dict
    :return: [description]
    :rtype: dict
    """
    all_raw_data = dict()
    for module_source, raw_data in mqc_rep.saved_raw_data.items():
        logging.info(f"Ingesting {module_source} from multiqc report into dataframe")
        sourced_raw_data = get_reformated_source_dict(raw_data)
        # merge into existing data
        for sample, new_data in sourced_raw_data.items():
            new_data = {
                (module_source, "general_stats", k): v for k, v in new_data.items()
            }
            if existing_dict := all_raw_data.get(sample):
                logging.debug(f"Adding additional data for existing sample: {sample}")
                existing_dict.update(new_data)
            else:
                logging.debug(f"Adding new sample: {sample}")
                all_raw_data[sample] = new_data

    # before final merging, force all sample names to change hypen to underscore
    # This occurs in at least one module in rseqc automatically
    logging.info(
        "Converting '-' in sample names to '_' to accommodate these forced changes in RSeQC"
    )

    def _normalize_dash(s):
        return s.replace("-", "_")

    all_raw_data = {
        tuple(_normalize_dash(sub_k) for sub_k in k): v for k, v in all_raw_data.items()
    }

    return all_raw_data

# A type hint class for the return from multiqc.run
class MQCRunDict(TypedDict):
    report: ModuleType
    config:  ModuleType
    sys_exit_code: int

def get_general_stats(mqc_run_output: MQCRunDict) -> dict[str, dict]:
    returnDict = dict()
    report = mqc_run_output['report']
    mqc_modules = [module.name for module in report.modules_output]
    for mqc_module, single_module_data in zip(mqc_modules, report.general_stats_data):
        returnDict[mqc_module] = single_module_data
    return returnDict

def format_as_dataframe(mqc_rep: MQCRunDict) -> pd.DataFrame:
    logging.info(f"Formatting to dataframe")
    # general stats data to columns
    mqc_modules = [module.name for module in mqc_rep.modules_output]


    flat_dict = flatten_raw_data(mqc_rep)

    # ingest plot data
    flat_plot_dict = format_plot_data(mqc_rep)
    # reformat to flatten list of dicts into single dict
    final_flat_plot_dict = dict()
    for s, list_of_dicts in flat_plot_dict.items():
        final_flat_plot_dict[s] = dict()
        for d in list_of_dicts:
            final_flat_plot_dict[s].update(d)

    # convert to sample wise dataframe and merge
    return pd.DataFrame(flat_dict).T.merge(
        pd.DataFrame(final_flat_plot_dict).T,
        left_index=True,
        right_index=True,
        how="outer",
    )


def parse_bar_graph_to_flat_dict(plot_data):
    # return messy sample:[{key (ylab):value}]
    val = {
        messy_s: [
            {
                (
                    plot_data["config"]["id"],
                    plot_data["config"]["id"],
                    f"{plot_data_subset['name']} ({plot_data['config'].get('ylab')})",
                ): plot_data_subset["data"][i]
            }
            for plot_data_subset in plot_data["datasets"][0]
        ]
        for i, messy_s in enumerate(plot_data["samples"][0])
    }
    return val


def __clean_mapped_data(mapped_data, messy_to_clean_map):
    __clean_mapped_data = defaultdict(lambda: list())
    for messy_s, data in mapped_data.items():
        # print(messy_s, data)
        # clean adapters from messy_s
        clean_s = messy_to_clean_map[messy_s]
        for datum in data:
            # print(datum)
            for k, v in datum.items():  # single item dicts
                __clean_mapped_data[
                    (clean_s["clean_sample"], clean_s["sub_source"])
                ].append({k: v})

                # __clean_mapped_data[clean_s["clean_sample"]].append({new_k: v})
    return __clean_mapped_data


def __parse_xy_line_graph_to_flat_dict(plot_data):
    # return messy sample:[{key (ylab):value}]
    all_flat_dict = dict()
    if "adapt" in plot_data["config"]["id"]:
        print(1)
        pass
    if categories := plot_data["config"].get("categories"):
        logging.debug("Plot has categorical data, extracting by category")
        for line in plot_data["datasets"][0]:
            sample_flat_dict = list()
            messy_s = line["name"]
            sample_flat_dict = [
                {
                    (
                        plot_data["config"]["id"],
                        plot_data["config"]["id"],
                        f"{categories[i_category]} {plot_data['config']['xlab']} ({plot_data['config']['ylab']})",
                    ): val
                }
                for i_category, val in enumerate(line["data"])
            ]
            all_flat_dict[messy_s] = sample_flat_dict
    else:
        logging.debug("Plot does not have categorical data, extracting accordingly")
        for line in plot_data["datasets"][0]:
            sample_flat_dict = list()
            messy_s = line["name"]
            # accomodate adaptor specific content (a feature of the column, not the sample)
            adapter_split_tokens = messy_s.split(" - ")
            if len(adapter_split_tokens) == 2:
                messy_s, adapter_s = adapter_split_tokens
            elif len(adapter_split_tokens) == 1:
                adapter_s = ""
            else:
                raise ValueError(
                    f"Unexpected 'messy' sample name in multiQC failed adapter related tokenization, this sample name raised the issue: {messy_s}"
                )
            sample_flat_dict = [
                {
                    (
                        plot_data["config"]["id"],
                        f"{plot_data['config']['id']}:{adapter_s}"
                        if adapter_s
                        else plot_data["config"]["id"],
                        f"{pos[0]} {plot_data['config']['xlab']} ({plot_data['config']['ylab']})",
                    ): pos[1]
                }
                for pos in line["data"]
            ]
            all_flat_dict[messy_s] = sample_flat_dict
    return all_flat_dict


def format_plot_data(mqc_rep: dict):
    logging.info(f"Attempting to extract data from {len(mqc_rep.plot_data)} plots")
    all_clean_data = dict()
    for plot_key, plot_data in mqc_rep.plot_data.items():
        logging.info(
            f"Attempting to extract data from plot with ID: {plot_data['config']['id']}"
        )
        logging.debug(f"Plot type: {plot_data['plot_type']}")
        # check plot type
        if plot_data["plot_type"] == "bar_graph":
            mapped_data = parse_bar_graph_to_flat_dict(
                plot_data
            )  # {messy_s: [{'sub_source::plot_name::sub_part--(units)':value}]
            messy_to_clean_map = {
                s.split(" - ")[
                    0
                ]: {  # this split effectively cleans adaptors from the string
                    "clean_sample": clean_messy_sample(s)[0],
                    "sub_source": clean_messy_sample(s)[1],
                }
                for s in plot_data["samples"][0]
            }
        elif plot_data["plot_type"] == "xy_line":
            mapped_data = __parse_xy_line_graph_to_flat_dict(
                plot_data
            )  # {messy_s: [{'sub_source::plot_name::sub_part--(units)':value}]
            messy_to_clean_map = {
                s.split(" - ")[
                    0
                ]: {  # this split effectively cleans adaptors from the string
                    "clean_sample": clean_messy_sample(s)[0],
                    "sub_source": clean_messy_sample(s)[1],
                }
                for s in [s["name"] for s in plot_data["datasets"][0]]
            }
        elif plot_data["plot_type"] in ["heatmap"]:
            logging.warning(
                f"Not implemented for dataframe extraction: {plot_data['plot_type']}, skipping this plot with ID: {plot_data['config']['id']}"
            )
        else:
            raise ValueError(
                f"Unexpected plot type encountered: {plot_data['plot_type']} for plot: {plot_key}"
            )

        clean_data = __clean_mapped_data(mapped_data, messy_to_clean_map)

        # before final merging, force all sample names to change hypen to underscore
        # This occurs in at least one module in rseqc automatically
        def _normalize_dash(s):
            return s.replace("-", "_")

        clean_data = {
            tuple(_normalize_dash(sub_k) for sub_k in k): v
            for k, v in clean_data.items()
        }

        # add to existing dict
        for clean_s, data in clean_data.items():
            if existing_data := all_clean_data.get(clean_s):
                data = data + existing_data
            all_clean_data[clean_s] = data
    return all_clean_data
