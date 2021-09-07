def add_num_elaspic_interface_entries(preliminary_data, protein_to_num_elaspic_interface_entries):
    """
    Given a tcga preliminary_data data, add corresponding number of ELASPIC interface entries for each protein.
    E.g.
        P04637 (TP53) → 76

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `NUM_ELASPIC_INTERFACE_ENTRIES` column will be to be added on.

        protein_to_num_elaspic_interface_entries : <defaultdict>
            A (default) dictionary that maps each protein to its number of ELASPIC interface entries.

    Returns
    -------
        None. Modifies the dataframe.
    """

    elaspic_interface_entries = []
    for protein in preliminary_data["PROTEIN"]:
        elaspic_interface_entries.append(protein_to_num_elaspic_interface_entries[protein])

    preliminary_data['NUM_ELASPIC_INTERFACE_ENTRIES'] = elaspic_interface_entries
