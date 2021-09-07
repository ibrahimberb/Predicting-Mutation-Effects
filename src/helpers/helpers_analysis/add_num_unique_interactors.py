def add_num_unique_interactors(preliminary_data, protein_to_num_unique_interactors):
    """                                          
    Given a tcga preliminary_data data, add corresponding number of unique ELASPIC interface interactors 
    for each protein.
    E.g.
        P04637 (TP53) → 15

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `NUM_UNIQUE_INTERACTORS` column will be to be added on.

        protein_to_num_unique_interactors : <defaultdict>
            A (default) dictionary that maps each protein to its number of unique interactors in ELASPIC interface data.


    Returns
    -------
        None. Modifies the dataframe.
    """

    num_unique_interactors = []
    for protein in preliminary_data["PROTEIN"]:
        num_unique_interactors.append(protein_to_num_unique_interactors[protein])

    preliminary_data['NUM_UNIQUE_INTERACTORS'] = num_unique_interactors
