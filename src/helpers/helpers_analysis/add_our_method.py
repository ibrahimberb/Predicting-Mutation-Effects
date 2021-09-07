def add_our_method(preliminary_data, our_method_dict):
    """
    Given a tcga preliminary_data data, add our_method counts.

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `OUR_METHOD` column will be to be added on.

        our_method_dict : <dict>
            A dictionary that maps proteins to our_method counts.

    Returns
    -------
        None. Modifies the dataframe.
    """

    our_method_counts = []
    for protein in preliminary_data["PROTEIN"]:
        our_method_counts.append(our_method_dict[protein])

    preliminary_data['OUR_METHOD'] = our_method_counts
