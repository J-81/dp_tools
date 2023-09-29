from pathlib import Path

import click

from dp_tools.core.utilites.metrics_extractor import (
    generate_extractor_from_yaml_config,
    AssayType,
)


CONFIG_YAML = "extraction_conf.yaml"
ISA_PARSE_PATH = "isa_config.yaml"

@click.command()
@click.option("--osd-id", help='OSD Accession ID. e.g. "OSD-194"', required=True)
def main(osd_id):
    ISA_PATH = list(Path.cwd().glob("*ISA*.zip"))[0]

    metricsExtractor = generate_extractor_from_yaml_config(config=CONFIG_YAML)

    metricsExtractor.extract_data_from_isa(
        accession=osd_id,
        isa_archive=ISA_PATH,
        config=Path(ISA_PARSE_PATH),
    )

    # metricsExtractor.append_manual_yaml_data(target_yaml=test_yaml)

    metricsExtractor.extract_sections()

    metricsExtractor.metrics.to_csv(f"{osd_id}_metrics.csv")

    metricsExtractor.process_metrics(assay_type=AssayType.bulkRNASeq).to_csv(f"{osd_id}_summary.csv")

if __name__ == '__main__':
    main()