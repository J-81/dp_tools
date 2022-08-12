# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.8rc

### Changed

- Updated GeneLab filename to url mapping to utilize the [GeneLab public API](https://genelab.nasa.gov/genelabAPIs)
  - Addresses removal of prior-used deprecated endpoints

## 1.0.7rc

### Dockerfile

- Added samtools as needed for certain checks

### Checks

#### Fixed

- check_contrasts_table_rows: message no longer introduces extra newlines into log

## Planned

## rc1.0.6

### Added

#### BulkRNASeq V&V Reporting

- A validation protocol that runs on a BulkRNASeq dataset model
- Includes generation of report files

#### BulkRNASeq Data Model From Nextflow RNASeq Concensus Pipeline

- A set of multi-stage loaders to create a data model
- Includes: validation system and multiQC powered data extraction

### Fixed

- Tilde characters are not converted to periods in contrasts: https://tower.nf/orgs/GL_Testing_Nextflow/workspaces/Nextflow_RCP_Testing/watch/1t8TfGbpDmCVNK
  - this should emulate R make.names behaviour completely

#### BulkRNASeq Reporter File Generation

- Data assets tagged with file categories for reporter file export including:
  - md5sum table
  - curation tables [GeneLab internal use]
