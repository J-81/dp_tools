""" Tools for finding specific paths 

These should rely on configuratble files to accommodate changing directory structures.
For now, the regexes will be hardcoded.

Such Locators SHOULD:
  - have find function that returns the path
  - search starting at root data directory
  - any RE patterns should be relative to the search_path
"""
from dataclasses import dataclass, field
import fnmatch
from functools import cache
import os
from pathlib import Path
from typing import ClassVar, List, Protocol, Set, Tuple, Union
import logging

log = logging.getLogger(__name__)

"""
TODO: implement as a class to achieve the following functionality:
- Track all unlocated data assets to report likely data assets missing from data-model
@dataclass
class Locator():
    all_data_assets: ClassVar[Set] = field(init=False, repr=False)
    root_dir: Path
    data_asset_config: dict
"""


def find_data_asset_path(
    data_asset_config: dict,
    root_dir: Path,
    search: bool = True,
    return_parsed_config_as_metadata: bool = True,
    glob: bool = False,
    **template_kwargs,
):  # type: ignore
    """Find a data asset
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
    search_pattern = Path(root_dir) / Path(
        os.path.join(*data_asset_config["processed location"]).format(**template_kwargs)
    )
    if not search:
        return search_pattern

    # get list of a files to check
    targets = _rIterDir(root_dir)

    # track if data asset is found
    found: Union[Path, bool]

    if glob:
        # glob search and check if one unique match found
        glob_found = fnmatch.filter([str(p) for p in targets], search_pattern)
        if len(glob_found) == 1:
            found = Path(glob_found[0])  # make sure to convert back into path
        elif len(glob_found) > 1:
            raise ValueError(
                f"Glob Pattern: {search_pattern} matched multiple targets: {glob_found}. One and only one must match."
            )
        else:
            found = False
    else:
        found = search_pattern if search_pattern in targets else False

    # return found or raise Exception
    if found != False:
        if return_parsed_config_as_metadata:
            return {
                "metadata": data_asset_config
                | {"template_kwargs": {**template_kwargs}},
                "path": found,
            }
        else:
            return found
    else:
        raise FileNotFoundError(f"Did not find {search_pattern} anywhere in {root_dir}")
