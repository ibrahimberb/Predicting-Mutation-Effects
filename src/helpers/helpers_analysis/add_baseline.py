def add_baseline(preliminary_data, baseline_dict):
    """
    Given a tcga preliminary_data data, add baseline counts.

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `BASELINE` column will be to be added on.

        baseline_dict : <dict>
            A dictionary that maps proteins to baseline counts.

    Returns
    -------
        None. Modifies the dataframe.
    """

    baseline_counts = []
    for protein in preliminary_data["PROTEIN"]:
        baseline_counts.append(baseline_dict[protein])

    preliminary_data['BASELINE'] = baseline_counts
