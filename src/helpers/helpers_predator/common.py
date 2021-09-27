import os
from pandas import DataFrame
import os.path
from datetime import datetime
from pathlib import Path

from ..mylogger import get_handler
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
        folder_path: Path,
        data: DataFrame,
        file_name: str,
        file_extension='csv',
        overwrite=False
) -> None:
    """
    A helper function to export given data with specified name and extension.
    """

    log.debug(f"Exporting data {file_name} at location {folder_path} ..")
    file_name = os.path.join(folder_path, file_name)
    file_date = datetime.today().strftime('%Y-%m-%d')
    file_name = f'{file_name}_{file_date}.{file_extension}'

    # Ensure the file is not exists before creating to prevent overwriting.
    if os.path.isfile(file_name) and not overwrite:
        log.error(f"File {file_name} is already exist.\n"
                  "To overwrite existing file, use `overwrite=True`.")

    else:
        # Export
        data.to_csv(file_name, index=False)
        log.info(f'{file_name} is exported successfully.')


def get_random_id(path):
    import uuid

    while True:
        random_id = uuid.uuid4().hex[:8]
        id_path = os.path.join(path, random_id)
        is_valid = not os.path.isdir(id_path)

        if is_valid:
            # Path(path).mkdir(parents=True)
            log.debug(f"Results folder with ID {random_id} is created.")
            return random_id


def export_prediction_data(
        tcga,
        data,
        file_name,
        folder_path,
        config,
        overwrite,
        voting,
        file_extension='csv'
) -> None:
    """
    A helper function to export given data with specified name and extension along with metadata.
    """
    current_date = datetime.today().strftime('%Y-%m-%d')
    current_time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

    random_id = get_random_id(folder_path)
    folder_name = os.path.join(f"{tcga}_prediction_{current_date}", random_id)
    log.debug(f"Exporting data {file_name} at location {folder_path} in folder {folder_name}..")

    folder_path = os.path.join(folder_path, folder_name)
    file_name = os.path.join(folder_path, file_name)
    file_name = f"{file_name}_{voting}_{current_date}.{file_extension}"

    Path(f"{folder_path}").mkdir(parents=True, exist_ok=True)

    # Ensure the file is not exists before creating to prevent overwriting.
    if os.path.isfile(file_name) and not overwrite:
        log.error(f"File {file_name} is already exist.\n"
                  "To overwrite existing file, use `overwrite=True`.")

    else:
        # Export data
        data.to_csv(file_name, index=False)
        log.info(f'{file_name} is exported successfully.')

        # Export metadata
        config_path = os.path.join(folder_path, "config.txt")
        with open(config_path, "w") as file:
            for config_name, config_dict in config.items():
                file.write(f"{config_name.upper()}\n")
                for key, values in config_dict.items():
                    file.write(
                        f"    {key.upper()}: {config[f'{config_name}'][f'{key}']}\n"
                    )

        log.debug(f"Config is exported.")
