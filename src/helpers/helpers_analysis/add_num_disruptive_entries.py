import pandas as pd


def add_num_disruptive_entries(preliminary_data: pd.DataFrame, prediction_data: pd.DataFrame):
    """
    Given a tcga preliminary_data data, find the number of disruptive predicted entries for each protein.
    E.g.
        P04637 (TP53) → 76 → 74
        
    Note: For some entries, it is possible that we may not have any predictions, either because it is 
        not an interface mutation or it was causing "invalid" predictions with its triplet (0 and 1 at the same time).

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `NUM_DISRUPTIVE_ENTRIES` column will be to be added on.

        prediction_data : <DataFrame>
            The dataframe which contains prediction column, along with protein, mutation, interactor columns.

    Returns
    -------
        None. Modifies the dataframe.
    """

    num_disruptive_entries = []
    for protein in preliminary_data["PROTEIN"]:
        # Get the disruptive entries, i.e. prediction is 0.
        disruptive_entries_data = prediction_data[(prediction_data['UniProt_ID'] == protein) &
                                                  (prediction_data['Predictions'] == 0)]
        num_disruptive_entries.append(len(disruptive_entries_data))

    preliminary_data['NUM_DISRUPTIVE_ENTRIES'] = num_disruptive_entries
