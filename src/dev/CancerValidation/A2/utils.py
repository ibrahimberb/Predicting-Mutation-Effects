from typing import List
import logging

from tqdm.notebook import tqdm

from src.helpers.mylogger import get_handler

from IPython.display import display

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class EntryNotFoundError(Exception):
    """Entry does not exists"""


class BidirectionalEntriesFoundError(Exception):
    """Both ways"""


class MultipleEntryFoundError(Exception):
    """There are multiple entries."""


def get_mutation_position(mutation):
    return mutation[1:-1]


class Direction:
    A_B = "A_to_B"
    B_A = "B_to_A"


class ClassLabels:
    DISRUPTIVE = "Disruptive"
    NON_DISRUPTIVE = "Non-disruptive"
    NOT_AVAILABLE = "N/A"


def _get_tcga_entry(tcga_data, protein, interactor):
    a_b_tcga_data = tcga_data[
        (tcga_data["UniProt_ID"] == protein) &
        (tcga_data["Interactor_UniProt_ID"] == interactor)
        ]

    b_a_tcga_data = tcga_data[
        (tcga_data["UniProt_ID"] == interactor) &
        (tcga_data["Interactor_UniProt_ID"] == protein)
        ]

    if len(a_b_tcga_data) != 0 and len(b_a_tcga_data) != 0:
        log.critical("Oh no.")
        display(a_b_tcga_data)
        raise BidirectionalEntriesFoundError

    # First data contains entry and the second one is empty
    elif len(a_b_tcga_data) != 0 and len(b_a_tcga_data) == 0:
        if len(a_b_tcga_data) != 1:
            log.critical("Oh no.")
            display(a_b_tcga_data)
            raise MultipleEntryFoundError

        [mut] = a_b_tcga_data["Mutation"]
        direction = Direction.A_B

        return a_b_tcga_data, mut, direction

    # First data is empty and the second one contains entry
    elif len(a_b_tcga_data) == 0 and len(b_a_tcga_data) != 0:
        if len(b_a_tcga_data) != 1:
            raise

        [mut] = b_a_tcga_data["Mutation"]
        direction = Direction.B_A

        return b_a_tcga_data, mut, direction

    # Both of them are empty
    else:
        raise EntryNotFoundError


def check_tcga(tcga_data, protein, interactor, p1_ires: List[str], p2_ires: List[str]):
    try:
        filtered_data, mutation, direction = _get_tcga_entry(
            tcga_data=tcga_data, protein=protein, interactor=interactor
        )
    except EntryNotFoundError:
        # print(f"{protein} {interactor} NOT FOUND IN TCGA.")
        return ClassLabels.NOT_AVAILABLE

    except MultipleEntryFoundError:
        return "i will handle these later."

    mut_position = get_mutation_position(mutation)

    if direction == Direction.A_B:
        if mut_position in p1_ires:
            return ClassLabels.DISRUPTIVE
        else:
            return ClassLabels.NON_DISRUPTIVE

    elif direction == Direction.B_A:
        if mut_position in p2_ires:
            return ClassLabels.DISRUPTIVE
        else:
            return ClassLabels.NON_DISRUPTIVE

    else:
        raise


def validate(validation_data, tcga_prediction_data):
    validation_results = []
    for index, row in tqdm(
        validation_data.iterrows(),
        total=len(validation_data)
    ):
        txt_protein = row["P1"]
        txt_interactor = row["P2"]
        p1_ires = row["P1_IRES"]
        p2_ires = row["P2_IRES"]

        result = check_tcga(
            tcga_data=tcga_prediction_data,
            protein=txt_protein,
            interactor=txt_interactor,
            p1_ires=p1_ires,
            p2_ires=p2_ires,
        )

        validation_results.append(result)

    return validation_results
