import copy
import gzip
from io import StringIO
import io
import os
from pathlib import Path
from unittest.mock import MagicMock
from dp_tools.bulkRNASeq.checks import COMPONENT_RAWREADS_0001, COMPONENT_TRIMREADS_0001
from dp_tools.bulkRNASeq.loaders import load_BulkRNASeq_STAGE_00
from dp_tools.core import check_model
from dp_tools.core.entity_model import DataFile

from pytest import MonkeyPatch

# set for testing
TEST_DIR = Path(os.environ["TEST_ASSETS_DIR"])

    

def test_COMPONENT_READS_0001():
    target_data_dir = TEST_DIR / "GLDS-194"
    ds = load_BulkRNASeq_STAGE_00(target_data_dir)
    
    test_component = list(ds.dataset.samples.values())[0].rawForwardReads

    # expected GREEN
    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    assert flag.code.name == 'GREEN'

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        m.setattr(test_component, "fastqGZ", None)
        m.setattr(test_component, "fastqcReportZIP", None)
        flag = COMPONENT_RAWREADS_0001.validate(test_component)
        assert flag.code.name == 'HALT1'
        assert flag.message == "Missing expected files: ['fastqGZ', 'fastqcReportZIP']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    #with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.code.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mock_path = io.BytesIO(b'bad file contents')
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = COMPONENT_RAWREADS_0001.validate(test_component_monkey_patch)
        assert flag.code.name == 'HALT3'
        assert flag.message == "Corrupted Fastq.gz file suspected, last line number encountered: 0"

def test_COMPONENT_TRIMREADS_0001():
    target_data_dir = TEST_DIR / "GLDS-194"
    ds = load_BulkRNASeq_STAGE_00(target_data_dir)
    
    test_component = list(ds.dataset.samples.values())[0].rawForwardReads

    # TODO: remove after implementing trimReads
    test_component.trimmingReportTXT = MagicMock(spec=DataFile)
    test_component.trimmingReportTXT.path = MagicMock(spec=Path)
    test_component.trimmingReportTXT.path.exists = lambda: True

    # expected GREEN
    flag = COMPONENT_TRIMREADS_0001.validate(test_component)
    assert flag.code.name == 'GREEN'

    # expected HALT1
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        m.setattr(test_component, "trimmingReportTXT", None)
        flag = COMPONENT_TRIMREADS_0001.validate(test_component)
        assert flag.code.name == 'HALT1'
        assert flag.message == "Missing expected files: ['trimmingReportTXT']"

    # expected HALT2
    # TODO: requires a true truncated test file or something similar
    #with MonkeyPatch.context() as m:
    #    m.setattr(gzip, "open", lambda: raise EOFError)
    #    flag = COMPONENT_RAWREADS_0001.validate(test_component)
    #    assert flag.code.name == 'HALT2'

    # expected HALT3
    with MonkeyPatch.context() as m:
        m.setattr(check_model, "ALLOWED_DEV_EXCEPTIONS", (SystemExit))
        mock_path = io.BytesIO(b'bad file contents')
        mock_path.exists = lambda: True
        test_component_monkey_patch = copy.deepcopy(test_component)
        test_component_monkey_patch.fastqGZ.path = mock_path
        flag = COMPONENT_TRIMREADS_0001.validate(test_component_monkey_patch)
        assert flag.code.name == 'HALT3'
        assert flag.message == "Corrupted Fastq.gz file suspected, last line number encountered: 0"
