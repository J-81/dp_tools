# This document explains usage from Gitpod; however, beside installation, these may (untested) also work when running from containers (wrapped in appropriate `singularity` or `docker` invocations)


## Installation

```
cd $REPO_DIRECTORY # e.g. /workspace/dp_tools in gitpod
pip install -e .
```

## Download Relevant MultiQC & ISA archive

> python download_multiqc_from_OSD.py --osd-id <OSD-NNN> --output-dir <OUTPUT_DIR>

### Known Limitations / Issues

* Only supports datasets with `read distribution` MultiQC files (used as a proxy for whether the dataset is actually sequencing  transcriptomics)
** Future: Should rely on parsing metadata from API

## Copy required configuration files

> bash set_up_config_files.sh <OUTPUT_DIR>

This copies template yaml files from the repository code.

## CD into directory

> cd <OUTPUT_DIR>

## Modify configuration files


### isa_config.yaml

1. Initially, no changes
2. If encountering error like: `ValueError: Could not find required column '['Parameter Value[Stranded]', 'Parameter Value[stranded]']' in either ISA sample or assay table.` 
  * Comment out or modify item in `Staging: ->  General: -> Required Metadata: -> From ISA:` section of yaml

### extraction_settings.yaml

1. MUST: change root search directory (line 2) to directory containing multiQC reports generated at start of this document
1. MAY: need to disable section for certain multiQC (not likely useful / will very probably break summarization)

## Run extract & summarize script

> python ../extract_dataset.py --osd-id <OSD_NNN> # You should still be in the directory with the multiQC outputs & yaml files

Outputs:

1. <OSD_NNN>_metrics.csv # Exhaustive metrics as pulled from multiQC reports
2. <OSD_NNN>_summary.csv # Summarization and derived statistics as generated on the exhaustive metrics table


## Overall Known Limitations

* Currently only supports paired end sequencing transcriptomics
   * Updating this will require updating both extraction & summarization code

* Certain ISA archives may not work
   * While most missing or encoded off-spec metadata can be addressed by disabling (commenting out) sections in `extraction_settings.yaml`, certain ones like missing `library layout` (unlikley but an example) will likely require more significant changes to accomodate.
