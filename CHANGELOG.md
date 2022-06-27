# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
