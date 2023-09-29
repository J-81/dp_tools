from pathlib import Path

import pytest

from dp_tools.core.utilites.metrics_extractor import (
    generate_extractor_from_yaml_config,
    AssayType,
)


@pytest.fixture
def test_yaml():
    # Make the path relative to this file
    TEST_DIR = Path(__file__).parent
    return TEST_DIR / "assets/test.yaml"


@pytest.fixture
def configuration_yaml():
    # Make the path relative to this file
    TEST_DIR = Path(__file__).parent
    return TEST_DIR / "assets/config.yaml"


@pytest.fixture
def OSD_576_metrics_csv():
    TEST_DIR = Path(__file__).parent
    return TEST_DIR / "assets/OSD-576_on_cluster_metrics.csv"


@pytest.fixture
def OSD_281_metrics_csv():
    TEST_DIR = Path(__file__).parent
    return TEST_DIR / "GLDS-281_on_cluster_metrics.csv"


def test_extract_general_information(test_yaml):
    MetricsExtractor().extract_general_information(assay_type=1, yaml_file=test_yaml)
    pass


def test_isa_to_yaml(glds194_test_dir, test_yaml, configuration_yaml):
    metricsExtractor = generate_extractor_from_yaml_config(config=configuration_yaml)

    metricsExtractor.extract_data_from_isa(
        accession="GLDS-194",
        isa_archive=glds194_test_dir / "Metadata/GLDS-194_metadata_GLDS-194-ISA.zip",
        config=Path("/workspace/metrics_bulkRNASeq.yaml"),
    )

    metricsExtractor.append_manual_yaml_data(target_yaml=test_yaml)

    metricsExtractor.extract_sections()

    assert (
        set(
            [
                "has_ERCC",
                "organism",
                "Tissue Type",
                "Library Prep Method",
                "PE or SE",
                "Stranded or Unstranded",
                "% rRNA contamination",
                "Original Sample Name",
                "OSD-#",
                "GLDS-#",
                "Data Source",
            ]
        ).difference(set(metricsExtractor.metrics))
        == set()
    )

    metricsExtractor.metrics.to_csv(glds194_test_dir / "test.csv")

    metricsExtractor.process_metrics(assay_type=AssayType.bulkRNASeq)


def test_load_and_process_metrics_table(configuration_yaml, OSD_576_metrics_csv):
    metricsExtractor = generate_extractor_from_yaml_config(config=configuration_yaml)

    metricsExtractor.load_metrics_csv(metrics_csv=OSD_576_metrics_csv)

    metricsExtractor.process_metrics(assay_type=AssayType.bulkRNASeq).to_csv(
        "/workspace/OSD_576_metrics_summary.csv"
    )
