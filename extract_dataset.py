from pathlib import Path

from dp_tools.core.utilites.metrics_extractor import (
    generate_extractor_from_yaml_config,
    AssayType,
)


CONFIG_YAML = "extraction_conf.yaml"
ACCESSION = "OSD-511"
ISA_PATH = "OSD-511_metadata_OSD-511-ISA.zip"
ISA_PARSE_PATH = "isa_config.yaml"

metricsExtractor = generate_extractor_from_yaml_config(config=CONFIG_YAML)

metricsExtractor.extract_data_from_isa(
    accession=ACCESSION,
    isa_archive=ISA_PATH,
    config=Path(ISA_PARSE_PATH),
)

# metricsExtractor.append_manual_yaml_data(target_yaml=test_yaml)

metricsExtractor.extract_sections()

metricsExtractor.metrics.to_csv("metrics.csv")

metricsExtractor.process_metrics(assay_type=AssayType.bulkRNASeq).to_csv("summary.csv")
