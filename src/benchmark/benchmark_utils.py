from typing import Tuple, Set

import pandas as pd
from pandas import DataFrame, Series
from datetime import datetime
import os.path as op
import re

from src.helpers.mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class ReadyResults:

    def __init__(self, downloaded_results_path):
        self.downloaded_results = pd.read_csv(downloaded_results_path, sep='\t', low_memory=False)

    def get_interface_ready_results(self):
        _, interface_data = self.get_core_and_interface_results()

        return interface_data

    def get_core_and_interface_results(self):
        # Drop entries where status is not 'done' (either 'err' or 'running'), i.e. keep only "done" entries.
        results_done = self.downloaded_results[self.downloaded_results["Status"] == "done"]

        core_data, interface_data = self.separate_core_interface(results_done)

        # Handle core mutations
        core_data_dropped = self.drop_duplicates(core_data)

        # Handle interface mutations
        interface_data_no_self_interactions = self.drop_self_interactions(interface_data)
        interface_data_dropped = self.drop_duplicates(interface_data_no_self_interactions)

        return core_data_dropped, interface_data_dropped

    @staticmethod
    def get_entries_core(data):
        core_data = data[data["Type"] == 'core'].reset_index(drop=True)
        return core_data

    @staticmethod
    def get_entries_interface(data):
        interface_data = data[data["Type"] == 'interface'].reset_index(drop=True)
        return interface_data

    def separate_core_interface(self, data) -> (DataFrame, DataFrame):
        """Returns Core Data and Interface Data"""
        print("Separating Core and Interface entries ..")
        core_data = self.get_entries_core(data)
        interface_data = self.get_entries_interface(data)

        print(f"Core data dimensions: {core_data.shape}")
        print(
            "Core data preview: \n{}\n".format(
                core_data[["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]].head(3)
            )
        )

        print(f"Interface data dimensions: {interface_data.shape}")
        print(
            "Interface data preview: \n{}\n".format(
                interface_data[["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]].head(3)
            )
        )

        return core_data, interface_data

    @staticmethod
    def drop_self_interactions(data):
        """
        Entries whose UniProt_ID protein is the same (or isoform) as Interactor_UniProt_ID will be removed.
        """

        print("Dropping self interactions ..")

        # Take the entries where UniProt_ID different than Interactor_UniProt_ID (dropping self and isoform)
        data_dropped = data[
            data["UniProt_ID"].apply(lambda x: x.split('-')[0]) != data["Interactor_UniProt_ID"].apply(
                lambda x: x.split('-')[0])
            ]

        # Reset index of the dataframe to avoid any possible errors
        data_dropped = data_dropped.reset_index(drop=True)

        return data_dropped

    @staticmethod
    def drop_duplicates(data):
        """
        Remove duplicated entries in the given dataframe.
        """

        print("Dropping duplicated entries ..")

        # Size of dataframe before dropping duplicated entries.
        print(f"Size of dataframe before dropping duplicated entries: {data.shape}")

        # Drop duplicates by keeping the 'first' one.
        data = data.drop_duplicates(keep="first")

        # Size of dataframe after dropping duplicated entries.
        print(f"Size of dataframe  after dropping duplicated entries: {data.shape}")

        # Reset index of the dataframe to avoid any possible errors
        data.reset_index(drop=True, inplace=True)

        return data


class SupplementaryExcelParser:
    """
    A class to handle io and parse operations in a convenient way.
    """

    def __init__(self, supp_excel_path):
        self.supp_excel_path = supp_excel_path
        self.pairs = None

    @staticmethod
    def get_input_pairs_via_series(
            protein_series: Series, mutation_series: Series
    ) -> Set[Tuple[str, str]]:
        """
        Get (protein,mutation) pairs using two series.
        """

        # Ensure both series provided correctly.
        assert len(protein_series) == len(mutation_series)
        assert (protein_series.index == mutation_series.index).all()

        pairs = set(zip(protein_series, mutation_series))

        return pairs

    @staticmethod
    def get_input_pairs_single_data(
            data, protein_column_name, mutation_column_name
    ) -> Set[Tuple[str, str]]:
        """
        Get (protein,mutation) pairs from single dataset.
        """
        data_protein_mutation = data[[protein_column_name, mutation_column_name]]
        pairs = set()
        for index, row in data_protein_mutation.iterrows():
            pairs.add((row[protein_column_name], row[mutation_column_name]))

        return pairs

    def extract_pairs(self, filename="ELASPIC_input"):

        file_date = datetime.today().strftime('%Y-%m-%d')
        filename = f"{filename}_{file_date}.txt"

        if self.pairs is None:
            raise NotImplementedError("Pairs are not filled.")

        if op.isfile(filename):
            raise FileExistsError(f"You already have the file {filename}")

        with open(filename, "w") as file_out:
            for pair in sorted(self.pairs):
                file_out.write(f"{pair[0]}.{pair[1]}\n")

        print("Input pairs are extracted.")

    @staticmethod
    def load_excel(excel_path, sheet_name):
        log.debug(f"Loading sheet {sheet_name} ..")
        data = pd.read_excel(excel_path, sheet_name=sheet_name)
        log.debug(f"Excel sheet {sheet_name} is loaded.")
        return data


def get_mutation_corrected(mutation):
    """Given mutation may have lower case letters. Here, convert it to proper form."""
    mutation_match = re.match(r"^([A-Za-z][0-9]+[A-Za-z]_?[0-9]*)$", mutation)
    if not mutation_match:
        # I could also change here to: # return None
        raise ValueError(f"Unexpected mutation countered: {mutation}")
    return mutation_match.group(1).upper()


def split_pairs_tuple(pairs):
    """
    Given a tuple `pairs` containing protein and mutation, return the
    :param pairs:
    :return:
    """
    return pairs[0], get_mutation_corrected(pairs[1])
