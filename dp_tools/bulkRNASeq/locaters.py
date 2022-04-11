""" Tools for finding specific paths 

These should rely on configuratble files to accommodate changing directory structures.
For now, the regexes will be hardcoded.

Such Locators SHOULD:
  - have find function that returns the path
  - search starting at root data directory
  - any RE patterns should be relative to the search_path
"""
from dataclasses import dataclass, field
from functools import cache
import os
from pathlib import Path
from typing import ClassVar, List, Protocol, Set, Tuple
import logging

log = logging.getLogger(__name__)

'''
TODO: implement as a class to achieve the following functionality:
- Track all unlocated data assets to report likely data assets missing from data-model
@dataclass
class Locator():
    all_data_assets: ClassVar[Set] = field(init=False, repr=False)
    root_dir: Path
    data_asset_config: dict
'''

def find_data_asset_path(
    data_asset_config: dict,
    root_dir: Path,
    search: bool = True,
    **template_kwargs
    ): # type: ignore
    """ Find a data asset 
    **template_kwargs include substrings like sample and dataset.
    This can be 'overloaded' but will throw an exception if a requested template value
    is not provided.

    Uses the data assets yaml format
    """
    @cache
    def _rIterDir(p: Path) -> Set[Path]:
        files = {f for f in p.rglob("*")}
        log.debug(f"Caching {len(files)} files found in {p}")
        return files        
    # format search pattern
    search_pattern = Path(root_dir) / Path(os.path.join(*data_asset_config["processed location"]).format(**template_kwargs))
    if not search:
        return search_pattern

    # get list of a files to check
    targets = _rIterDir(root_dir)

    # find if matches, else return None
    if search_pattern in targets:
        return search_pattern
    else:
        raise FileNotFoundError(f"Did not find {search_pattern} anywhere in {root_dir}")
