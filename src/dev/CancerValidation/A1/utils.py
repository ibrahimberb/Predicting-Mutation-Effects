import re
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


class EntryNotFoundError(Exception):
    pass


class BidirectionalEntriesFoundError(Exception):
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
            [p1_res] = a_b_interface_data["P1_IRES"]
            return a_b_interface_data, p1_res

        # First data is empty and the second one contains entry
        elif len(a_b_interface_data) == 0 and len(b_a_interface_data) != 0:
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
        num_entries = len(tcga_data)
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

