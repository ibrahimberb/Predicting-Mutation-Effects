import pickle
import re
from pathlib import Path
from typing import List

from IPython.display import display

import os.path as op
from datetime import datetime

import pandas as pd
from tqdm.notebook import tqdm

from matplotlib import pyplot as plt

from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
)

from src.helpers.mylogger import get_handler

import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class EntryNotFoundError(Exception):
    pass


class BidirectionalEntriesFoundError(Exception):
    pass


class MultipleEntryFoundError(Exception):
    pass


class ClassLabels:
    DISRUPTIVE = "Disruptive"
    NON_DISRUPTIVE = "Non-disruptive"
    NOT_AVAILABLE = "N/A"


def unzip_res_range(res_range):
    """
    Converts ranges in the form: [2-210] or [3-45,47A,47B,51-67] into lists of strings including all numbers in
    these ranges in order .
    Items are of type <str>.
    """
    res_ranges = res_range.strip()[1:-1].split(',')
    index_list = []
    for r in res_ranges:
        if re.match('.+-.+', r):
            a, b = r.split('-')
            index_list += [str(n) for n in range(int(a), int(b) + 1)]
        else:
            index_list.append(r)

    return index_list


# print(unzip_res_range("[95-96,98-100,102-103,122,262,266-267,270,273,294]"))


def get_mutation_position(mutation):
    return mutation[1:-1]


class CancerValidation:
    def __init__(self, interfaces_data_path):
        self.interfaces_data = self.load_data(interfaces_data_path)

    @staticmethod
    def load_data(data_path):
        interfaces_data = pd.read_csv(
            data_path, sep="\t", usecols=["P1", "P2", "P1_IRES", "P2_IRES"]
        )
        interfaces_data["P1_IRES"] = interfaces_data["P1_IRES"].apply(lambda x: unzip_res_range(x))
        interfaces_data["P2_IRES"] = interfaces_data["P2_IRES"].apply(lambda x: unzip_res_range(x))
        return interfaces_data

    def check(self, protein: str, mutation: str, interactor: str):
        # print(f"Checking ..\n"
        #       f"> PROTEIN: {protein} \n"
        #       f"> MUTATION: {mutation} \n"
        #       f"> INTERACTOR: {interactor}")
        try:
            data, res = self._get_entry(protein, interactor)
        except EntryNotFoundError:
            return ClassLabels.NOT_AVAILABLE

        mut_pos = get_mutation_position(mutation)
        if mut_pos in res:
            return ClassLabels.DISRUPTIVE
        else:
            return ClassLabels.NON_DISRUPTIVE

    @staticmethod
    def _handle_check_duplicated_entries(data, p_ires) -> List[int]:
        """
        Checks if all entries of given data duplicated. Each cell contains list in it.
        If all entries are duplicated entries, then we have no problem, just get the res.
        """
        # display(data)

        # In order for us to check if the entries are duplicated, we'll have to
        # convert list item in the cells to tuple. Otherwise, we get the following error:
        # TypeError: unhashable type: 'list'
        data_tuple = data[["P1_IRES", "P2_IRES"]].applymap(lambda x: tuple(x))
        data_tuple.duplicated(keep=False).all()

        # no problem, then.
        if p_ires == "P1_IRES":
            [p1] = data["P1"].unique()
            [p2] = data["P2"].unique()
        elif p_ires == "P2_IRES":
            [p2] = data["P1"].unique()
            [p1] = data["P2"].unique()
        else:
            raise ValueError(f"Illegal argument provided for parameter `p_ires`: {p_ires}")

        # check if all entries are duplicated
        if data_tuple.duplicated(keep=False).all():

            log.warning(
                f"Multiple entries but they were duplicated. PROTEIN: {p1}, INTERACTOR: {p2}"
            )
            [p_res] = data_tuple[p_ires].unique()
            p_res = list(p_res)

            return p_res

        else:
            log.error("MultipleEntryError with following data: ")
            display(data)
            p_res_list = data[p_ires].tolist()
            p_res = sorted(
                set([item for sublist in p_res_list for item in sublist])
            )
            log.error(F"Returned RES: {p_res}")

            return p_res
            # data.to_csv("ERROR_data.csv", index=False)
            # raise MultipleEntryFoundError

    def _get_entry(self, protein, interactor):
        a_b_interface_data = self.interfaces_data[
            (self.interfaces_data["P1"] == protein) &
            (self.interfaces_data["P2"] == interactor)
            ]

        b_a_interface_data = self.interfaces_data[
            (self.interfaces_data["P1"] == interactor) &
            (self.interfaces_data["P2"] == protein)
            ]

        # Both of them contains entry -- this is an unlikely situation, unless there is problem with the text file..
        if len(a_b_interface_data) != 0 and len(b_a_interface_data) != 0:
            raise BidirectionalEntriesFoundError

        # First data contains entry and the second one is empty
        elif len(a_b_interface_data) != 0 and len(b_a_interface_data) == 0:
            if len(a_b_interface_data) != 1:
                p1_res = self._handle_check_duplicated_entries(a_b_interface_data, "P1_IRES")

            else:
                [p1_res] = a_b_interface_data["P1_IRES"]

            return a_b_interface_data, p1_res

        # First data is empty and the second one contains entry
        elif len(a_b_interface_data) == 0 and len(b_a_interface_data) != 0:
            if len(b_a_interface_data) != 1:
                p2_res = self._handle_check_duplicated_entries(b_a_interface_data, "P2_IRES")

            else:
                [p2_res] = b_a_interface_data["P2_IRES"]

            return b_a_interface_data, p2_res

        # Both of them are empty
        else:
            raise EntryNotFoundError

    def validate(
            self,
            tcga_type: str,
            tcga_data: pd.DataFrame,
    ):
        tcga_data_validation = tcga_data.copy()
        validation_results = []
        for index, row in tqdm(
                tcga_data_validation.iterrows(),
                total=len(tcga_data_validation)
        ):
            protein = row["UniProt_ID"]
            mutation = row["Mutation"]
            interactor = row["Interactor_UniProt_ID"]
            valid_label = self.check(
                protein=protein, mutation=mutation, interactor=interactor
            )
            # print(f">> RESULT: {valid_label}")
            validation_results.append(valid_label)

        tcga_data_validation["Validation"] = validation_results
        tcga_data_validation["Validation"].value_counts().plot(
            kind="bar", title=f"{tcga_type} Validation Results"
        )
        plt.show()

        tcga_data_validation_processed = process_validation_data(tcga_data_validation)
        tcga_data_validation_processed["Validation"].value_counts().plot(
            kind="bar", title=f"{tcga_type} Validation Processed Results"
        )
        plt.show()

        metrics_data = get_scoring_metrics(tcga_data_validation_processed)
        num_entries = len(tcga_data_validation)
        counts = tcga_data_validation_processed["Validation"].value_counts().to_dict()
        num_disruptive = counts[0]
        num_non_disruptive = counts[1]
        metrics_data.insert(0, "TCGA", tcga_type)
        metrics_data.insert(1, "#_Entries", num_entries)
        metrics_data.insert(2, "#_Disruptive", num_disruptive)
        metrics_data.insert(3, "#_Non_disruptive", num_non_disruptive)

        return {
            "data_validation": tcga_data_validation,
            "data_validation_processed": tcga_data_validation_processed,
            "metrics_data": metrics_data,
        }

    @staticmethod
    def validate_single_class(
            tcga_type: str,
            output_already_calculated: dict,
            single_class: int,
    ):
        """
        Requires the positions to be already calculated.
        """

        tcga_data_validation = output_already_calculated["data_validation"]

        print(f"Using the class {single_class} only.")
        tcga_data_validation = tcga_data_validation[
            tcga_data_validation["Prediction"] == single_class
            ].copy()

        tcga_data_validation_processed = process_validation_data(tcga_data_validation)
        metrics_data = get_scoring_metrics(tcga_data_validation_processed)
        num_entries = len(tcga_data_validation)
        counts = tcga_data_validation_processed["Validation"].value_counts().to_dict()
        num_disruptive = counts[0]
        num_non_disruptive = counts[1]
        metrics_data.insert(0, "TCGA", tcga_type)
        metrics_data.insert(1, "#_Entries", num_entries)
        metrics_data.insert(2, "#_Disruptive", num_disruptive)
        metrics_data.insert(3, "#_Non_disruptive", num_non_disruptive)

        return {
            "data_validation": tcga_data_validation,
            "data_validation_processed": tcga_data_validation_processed,
            "metrics_data": metrics_data,
        }

    @staticmethod
    def extract_output_dict(name, dict_obj):
        folder_path = "outputs"
        Path(f"{folder_path}").mkdir(parents=True, exist_ok=True)
        current_date = datetime.now().strftime("%Y-%m-%d")
        file_name = f"{name}_{current_date}.pickle"
        file_path = op.join(folder_path, file_name)
        if op.exists(file_path):
            raise FileExistsError("File already exists")

        pickle.dump(dict_obj, open(file_path, "wb"))
        print("Object extracted successfully.")

    @staticmethod
    def load_output_dict(pickle_path):
        obj_loaded = pickle.load(open(pickle_path, "rb"))
        return obj_loaded


def test_entry_not_found(df, p, i):
    a_b = df[
        (df["P1"] == p) &
        (df["P2"] == i)
        ]

    b_a = df[
        (df["P1"] == i) &
        (df["P2"] == p)
        ]

    assert len(a_b) == len(b_a) == 0


def get_scoring_metrics(tcga_validation_data):
    y_true = tcga_validation_data["Validation"]
    y_pred = tcga_validation_data["Prediction"]
    metrics_data = pd.DataFrame(
        [
            accuracy_score(y_true, y_pred),
            balanced_accuracy_score(y_true, y_pred),
            f1_score(y_true, y_pred),
            precision_score(y_true, y_pred),
            recall_score(y_true, y_pred),
            matthews_corrcoef(y_true, y_pred),
        ],
        index=["ACCURACY", "BALANCED_ACCURACY", "F1", "PRECISION", "RECALL", "MATTHEWS_COR"]
    ).T

    return metrics_data


def process_validation_data(tcga_data: pd.DataFrame):
    """
    Process the validation data.
    1. Drop N/A entries
    2. Convert Labels as follows:
        DISRUPTIVE → 0
        NON_DISRUPTIVE → 1
    3. Convert its type to int.
    """

    tcga_processed = tcga_data[tcga_data["Validation"] != "N/A"].copy()
    tcga_processed["Validation"] = tcga_processed["Validation"].replace(
        {
            ClassLabels.DISRUPTIVE: 0,
            ClassLabels.NON_DISRUPTIVE: 1,
        }
    )
    tcga_processed = tcga_processed.astype({"Validation": "int"})

    return tcga_processed
