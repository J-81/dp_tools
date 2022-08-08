""" Test CLI scripts """


import os


def test_dpt_get_isa_archive(script_runner, tmpdir):
    os.chdir(tmpdir)

    ret = script_runner.run("dpt-get-isa-archive", "--accession", "GLDS-194")
    assert ret.success


def test_dpt_isa_to_runsheet(script_runner, tmpdir, glds194_test_dir):
    os.chdir(tmpdir)
    isaPath = glds194_test_dir / "Metadata" / "GLDS-194_metadata_GLDS-194-ISA.zip"

    ret = script_runner.run(
        "dpt-isa-to-runsheet",
        "--accession",
        "GLDS-194",
        "--config-type",
        "bulkRNASeq",
        "--config-version",
        "Latest",
        "--isa-archive",
        str(isaPath),
    )
    assert ret.success
