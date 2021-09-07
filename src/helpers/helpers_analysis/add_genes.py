def add_genes(preliminary_data, protein_to_gene_dict):
    """
    Given a tcga preliminary_data data, add corresponding genes.
    E.g.
        P62837 → UBE2D2

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `GENES` column will be to be added on.

        protein_to_gene_dict : <dict>
            A dictionary that mappes protein name to gene name. 
            One of the followings:
            - brca_protein_to_gene_dict
            - coad_protein_to_gene_dict
            - ov_protein_to_gene_dict

    Returns
    -------
        None. Modifies the dataframe.
    """

    genes = []
    for protein in preliminary_data["PROTEIN"]:
        genes.append(protein_to_gene_dict[protein])

    preliminary_data['GENE'] = genes
