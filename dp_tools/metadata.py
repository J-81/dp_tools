""" Protocols between metadata files """
from pathlib import Path
from typing import List, Protocol, runtime_checkable

import pandas as pd


@runtime_checkable
class BulkRNASeqMetadata(Protocol):
    paired_end: bool

    def __init__(self):
        pass

    def samples(self) -> List[str]:
        pass


class Runsheet:
    def __init__(self, path: Path):
        self.df = pd.read_csv(path)
        self.paired_end = self.df["paired_end"].unique()[0]

    def samples(self):
        return self.df["sample_name"]


class DummyMetadata:
    def samples(self):
        return ["testS1", "testS2", "testS3"]

