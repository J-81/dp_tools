""" 
## Description

dp_tools is a collection of data processing tools designed to assist in GeneLab data processing operations

## What This Library Includes

This module currently provides the following major features:

1. [General Validation and Verification (V&V) Framework](dp_tools/core/check_model.html)
2. General Data Model Framework ([base data model](dp_tools/core/entity_model.html), [data model loading functions](dp_tools/core/loaders.html))
3. [Yaml Configuration Interface](dp_tools/config/interface.html)
4. [GLDS API Wrapper Functions](dp_tools/glds_api/commons.html) 

Additionally, the following Assay specific functionality is packaged:

1. bulkRNASeq
  - [configuration files](https://github.com/J-81/dp_tools/tree/development/dp_tools/config)
  - [check functions](dp_tools/bulkRNASeq/checks.html)
  - [validation procotol](dp_tools/bulkRNASeq/vv_protocols.html)

## Installation

#### Using Containers (e.g. Singularity)

This library is available for usage as prebuilt images located at [quay.io](https://quay.io/repository/j_81/dp_tools?tab=tags)
> singularity shell quay.io/repository/j_81/dp_tools

#### Using pip

> pip install git+https://github.com/J-81/dp_tools.git@1.3.4

## CLI Commands

**Note: Most library functionality is only available through using python import.**

The following command line scripts are also available once installed and are defined here:

- dpt-get-isa-archive
``` bash
usage: dpt-get-isa-archive [-h] --accession GLDS-001

Script for downloading latest ISA from GLDS repository

options:
  -h, --help            show this help message and exit
  --accession GLDS-001  GLDS accession number
```

- dpt-isa-to-runsheet
``` bash
usage: dpt-isa-to-runsheet [-h] --accession GLDS-001 --config-type CONFIG_TYPE [--config-version CONFIG_VERSION] --isa-archive ISA_ARCHIVE

Script for downloading latest ISA from GLDS repository

options:
  -h, --help            show this help message and exit
  --accession GLDS-001  GLDS accession number
  --config-type CONFIG_TYPE
                        Packaged config type to use. Currently supports: ['microarray', 'bulkRNASeq']
  --config-version CONFIG_VERSION
                        Packaged config version to use
  --isa-archive ISA_ARCHIVE
                        Local location of ISA archive file. Can be downloaded from the GLDS repository with 'dpt-get-isa-archive'
```

## Examples

### Two step process to create runsheet from GLDS ISA Archive (using Singularity)

> **Note: At this time, No stdout messages are printed for these scripts**

#### Download ISA Archive
``` bash
# First two lines tell Singularity to run the dp_tools container in the current working directory
singularity exec --bind $(pwd):$(pwd) \\
  docker://quay.io/j_81/dp_tools:1.3.4 \\
  dpt-get-isa-archive --accession GLDS-168 # command we want to run
```

#### Convert ISA Archive into bulkRNASeq Runsheet
``` bash
# First two lines tell Singularity to run the dp_tools container in the current working directory
singularity exec --bind $(pwd):$(pwd) \\
  docker://quay.io/j_81/dp_tools:1.3.4 \\
  dpt-isa-to-runsheet --accession GLDS-168 \\
                      --config-type bulkRNASeq \\
                      --config-version Latest \\
                      --isa-archive GLDS-168_metadata_GLDS-168-ISA.zip # command we want to run
```

"""
__version__ = "1.3.4"
