# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.4]

### Changed

- Table updates (associated with updating ISA archive files) now separates multiple files in a field with ',' instead of ', '

## [1.3.3]

### Added

- Support for data asset key sets and run components in updated validation interface (i.e. by 'dpt validation')

## [1.3.2]

### Fixed

- Refactored ISA archive parsing functions as prior the fallback wasn't being used in all calls (specifically the plug in based ones)

## [1.3.1]

### Fixed

- Parsing for ISA Archives met 'ISO-8859-1' encoding but not 'utf-8'
  - Specifically, 'utf-8' is attempted and 'ISO-8859-1' is used as a fallback

## [1.3.0]

### Added

- Improved V&V interface
  - Plugin support for protocols & checks
  - Specification generation via `dpt validation spec`
  - Generalized validation runner via `dpt validation run`
  - Manual check interface via `dpt validation manual-checks`
  - CLI interface implemented using [click](https://click.palletsprojects.com/en/8.1.x/)

- Started to use logging via [loguru](https://github.com/Delgan/loguru)

### Removed

- Microarray specific protocols/checks (now implemented as plugins outside dp_tools)

## [1.2.1]

### Added

- Ability to inject columns during runsheet generation

## [1.2.0]

### Added

- Microarray (Agilent 1 Channel) V&V protocol
- Pandera as dependency for better validation tooling

### Changed

- BulkRNASeq runsheet validation enhanced
  - Upgraded from Schema to Pandera
  - Added checks for dataset metadata columns like 'paired_end'
  - Added sanity check for 'read2_path' column optional nature

### Fixed

- [#19](https://github.com/J-81/dp_tools/issues/19)

## [1.1.9]

### Added
- Runsheet generation for methlySeq ISA archives

## [1.1.8]

### Changed
- GLDS API usage now considers the 'OSD' accession ID as the study ID instead of 'GLDS'.  This is consistent with the recent release of the [OSDR](https://osdr.nasa.gov/bio/)
### Fixed
- Fixes incorrect numeric inferrence for strings (commit: 3b0d953)[https://github.com/J-81/dp_tools/commit/3b0d9537de73363aaa78979b78b3a209c69ccd45]

## [1.1.7]

### Fixed
- Fixes incorrect unit detection for runsheet generation [#14](https://github.com/J-81/dp_tools/issues/14)
## [1.1.6]

### Added
- Stdout logging for scripts, this better explains what is happening during the script

### Fixed
- Missing Microarray technology valid combination and handling of multiple valid combinations
## [1.1.5]

### Fixed
- Staging runsheets failing to extract unit columns 
- V&V crash related to factor columns being inferred as numeric. Now correctly inferring as string values.

## [1.1.4]
### Added
- Integrity check for gzipped files to bulkRNASeq checks and protocol

### Changed
- Pinned Pandas version to 1.4.4 (prior: no pin, most recent version installed)
  - Version 1.5 causes changes to checksum for pandas objects and would require updating all tests that include a checksum (planned for future)
## [1.1.3]

### Fixed
- Fixing false V&V halt flagging: Add in micro sign as whitelisted (better in sync with r make.names function)
- Expected location of SampleTable.csv and ERCC_SampleTable.csv in 
## [1.1.2]

### Fixed
- Fixing false V&V halt flagging: Add in greek characters as whitelisted (better in sync with r make.names function)

## [1.1.1]
### Fixed
- Incorrect detection of has_ERCC from ISA Archives
  - Example Impacted GLDS: 161,162,163,173
- Runsheet generation failing for different variations of raw reads data column names
  - Example Impacted GLDS: 105,138

## 1.1.0

### First Production Release
- Prior 1.0.0 tagged versions were actually develop style releases
- Moving ahead only production releases will have tags without 'rc' (release candidate) in the name

### Quality Updates
- Various flag messages improved
- Documentation updated

### Added
- Check related to multiQC samples inclusion

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

[1.1.1]: https://github.com/j-81/dp_tools/compare/1.1.0...1.1.1
[1.1.2]: https://github.com/j-81/dp_tools/compare/1.1.1...1.1.2
[1.1.3]: https://github.com/j-81/dp_tools/compare/1.1.2...1.1.3
[1.1.4]: https://github.com/j-81/dp_tools/compare/1.1.3...1.1.4
[1.1.5]: https://github.com/j-81/dp_tools/compare/1.1.4...1.1.5
[1.1.6]: https://github.com/j-81/dp_tools/compare/1.1.5...1.1.6
[1.1.7]: https://github.com/j-81/dp_tools/compare/1.1.6...1.1.7
[1.1.8]: https://github.com/j-81/dp_tools/compare/1.1.7...1.1.8
[1.1.9]: https://github.com/j-81/dp_tools/compare/1.1.8...1.1.9
[1.2.0]: https://github.com/j-81/dp_tools/compare/1.1.9...1.2.0
[1.2.1]: https://github.com/j-81/dp_tools/compare/1.2.0...1.2.1
[1.3.0]: https://github.com/j-81/dp_tools/compare/1.2.1...1.3.0
[1.3.1]: https://github.com/j-81/dp_tools/compare/1.3.0...1.3.1
[1.3.2]: https://github.com/j-81/dp_tools/compare/1.3.1...1.3.2
[1.3.3]: https://github.com/j-81/dp_tools/compare/1.3.2...1.3.3
[1.3.4]: https://github.com/j-81/dp_tools/compare/1.3.3...1.3.4
