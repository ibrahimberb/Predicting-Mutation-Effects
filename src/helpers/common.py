import os
import pandas as pd
from pandas import DataFrame


def print_annotation(s):
    print(f"\n{s}\n{'-' * len(s)}")


def get_file_path(project_common_file_dir, filename):
    return os.path.join(project_common_file_dir, filename)


def extract_data(filename: str, data: DataFrame, extension: str, overwrite=False):
    """
    Extract the given dataframe
    # todo add time and date.
    :param overwrite:
    :param extension:
    :param filename:
    :param data:
    :return:
    """
