import pandas as pd


def test_dataset_accessors(glds194_dataSystem):
    ds = glds194_dataSystem.dataset

    ds.get_assets()

    ds.get_assets(filter_to=["runsheet"])


def test_bulkRNASeq_mqc_api(glds194_dataSystem):
    # check plot name retrieval
    ds = glds194_dataSystem.dataset
    data = ds.compile_multiqc_data(
        data_asset_keys=[
            "raw reads fastQC ZIP",
            "raw forward reads fastQC ZIP",
            "raw reverse reads fastQC ZIP",
        ]
    )
    data = ds.compile_multiqc_data()
    plots = data["plots"]["FastQC"]
    assert set(plots) == set(
        [
            "Overrepresented sequences",
            "Per Sequence GC Content",
            "Per Base N Content",
            "Adapter Content:Subplot::illumina_universal_adapter",
            "Sequence Counts",
            "Mean Quality Scores",
            "Per Sequence Quality Scores",
            "Sequence Length Distribution",
            "Sequence Duplication Levels",
        ]
    )

    assert isinstance(data["general_stats"]["FastQC"], pd.DataFrame)
