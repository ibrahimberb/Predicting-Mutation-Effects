import os
import pandas as pd
from pandas import DataFrame
import os.path
from datetime import datetime

from .mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def print_annotation(s):
    print(f"\n{s}\n{'-' * len(s)}")


def get_file_path(project_common_file_dir, filename):
    return os.path.join(project_common_file_dir, filename)


def export_data(
        data: DataFrame,
        file_name: str,
        file_extension='csv',
        overwrite=False
):
    """
    :param data:
    :param file_name:
    :param file_extension:
    :param overwrite:
    :return:
    """

    log.debug(f"Exporting data {file_name} ..")

    file_date = datetime.today().strftime('%Y-%m-%d')
    file_name = f'{file_name}_{file_date}.{file_extension}'

    # Ensure the file is not exists before creating to prevent overwriting.
    if os.path.isfile(file_name) and not overwrite:  ## fixme: i am not sure if logic is correct
        print(f'File {file_name} is already exist.')

    else:
        # Export
        data.to_csv(file_name, index=False)
        log.debug(f'{file_name} is exported successfully.')

