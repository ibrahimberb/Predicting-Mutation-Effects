from datetime import datetime

from ..mylogger import get_handler
import logging


import os
import pathlib

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def save_auc_scores(prediction_file_path, file_name, auc_scores, overwrite):
    folder_path = pathlib.Path(prediction_file_path).parent
    file_name = os.path.join(folder_path, file_name)
    file_date = datetime.today().strftime('%Y-%m-%d')
    file_name = f'{file_name}_{file_date}.txt'

    # Ensure the file is not exists before creating to prevent overwriting.
    if os.path.isfile(file_name) and not overwrite:
        log.error(f"File {file_name} is already exist.\n"
                  "To overwrite existing file, use `overwrite=True`.")

    else:
        data_entries = []
        with open(file_name, 'w') as file:
            for analysis_type, auc_scores in auc_scores.items():
                log.debug(f"{analysis_type}: {auc_scores}")
                file.write(f"{analysis_type}: {auc_scores}\n")
                data_entries.append((analysis_type, auc_scores))

        log.info(f"AUC scores are saved into file {file_name}")

        return data_entries
