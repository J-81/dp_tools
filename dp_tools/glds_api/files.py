""" Module for finding file urls as hosted on the GLDS repository """
import functools
from typing import Union
from urllib.parse import quote
from dp_tools.core.entity_model import DataFile, TemplateDataset, TemplateSample
from dp_tools.glds_api.commons import GENELAB_ROOT, get_glds_filelisting_json
import logging

log = logging.getLogger(__name__)


def get_valid_url(fn: str, version: str) -> str:
    return f"{GENELAB_ROOT}/genelab/static/media/dataset/{quote(fn)}?version={version}"


def reformat_filelisting(original_json: dict) -> dict:
    """Reformats the original file listing

    :param original_json: _description_
    :type original_json: dict
    :return: _description_
    :rtype: dict
    """
    return {
        data["file_name"]: {"url": get_valid_url(data["file_name"], data["version"])}
        for data in original_json
    }


@functools.cache
def get_urls(accession: str) -> dict[str, str]:
    file_listing_json = get_glds_filelisting_json(accession)
    return reformat_filelisting(file_listing_json)


def associate_url_for_data_file(
    entity: Union[TemplateDataset, TemplateSample], target_asset: DataFile
) -> str:
    match entity:
        case TemplateDataset():
            dataset = entity
        case TemplateSample():
            dataset = entity.dataset
    urls = get_urls(accession=dataset.dataSystem.name)
    glds_name = f"{dataset.config['ISA Meta']['Global file prefix'].format(datasystem=dataset.dataSystem.name)}{target_asset.path.name}"
    target_asset.glds_url = urls[glds_name]["url"]
    return urls[glds_name]["url"]
