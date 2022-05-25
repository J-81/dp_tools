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


@cache
def _rIterDir(p: Path) -> Set[Path]:
    files = {f for f in p.rglob("*")}
    log.debug(f"Caching {len(files)} files found in {p}")
    return files


@dataclass
class Locator:

    root_dir: Path
    data_asset_config: dict
    validation_enabled: bool = field(default=True)
    return_parsed_config_as_metadata: bool = field(default=True)

    def find_data_asset_path(
        self,
        config_key: str,
        search: bool = True,
        glob: bool = False,
        **template_kwargs,
    ):  # type: ignore
        """Find a data asset
        **template_kwargs include substrings like sample and dataset.
        This can be 'overloaded' but will throw an exception if a requested template value
        is not provided.

        Uses the data assets yaml format
        """
        this_data_asset_config = self.data_asset_config[config_key]

        # format search pattern
        search_pattern = Path(self.root_dir) / Path(
            os.path.join(*this_data_asset_config["processed location"]).format(
                **template_kwargs
            )
        )

        # check a number of conditions to see if validation is disabled
        if (
            any(
                [
                    (
                        not this_data_asset_config.get("validate exists", True)
                    ),  # specific data asset config
                    (
                        not self.validation_enabled
                    ),  # global control, runtime initialized
                ]
            )
            and not search
        ):  # runtime controlled, overrides other config
            log.debug(f"Skipping validation as configured for {search_pattern}")
            if self.return_parsed_config_as_metadata:
                return {
                    "metadata": this_data_asset_config
                    | {"template_kwargs": {**template_kwargs}}
                    | {"data asset key": config_key},
                    "path": search_pattern,
                }
            else:
                return search_pattern

        # get list of a files to check
        targets = _rIterDir(self.root_dir)

        # track if data asset is found
        found: Union[Path, bool]

        if glob:
            # glob search and check if one unique match found
            glob_found = fnmatch.filter([str(p) for p in targets], search_pattern)  # type: ignore
            if len(glob_found) == 1:
                found = Path(
                    glob_found[0]  # type: ignore
                )  # make sure to convert back into path
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
            if self.return_parsed_config_as_metadata:
                return {
                    "metadata": this_data_asset_config
                    | {"template_kwargs": {**template_kwargs}}
                    | {"data asset key": config_key},
                    "path": found,
                }
            else:
                return found
        else:
            raise FileNotFoundError(
                f"Did not find {search_pattern} anywhere in {self.root_dir}"
            )
