def add_cancermine_status(preliminary_data, cancermine_genes, cancermine_cohort_genes, tcga_type):
    """
    Given a tcga preliminary_data data, add the CancerMine and CancerMine_Cohort status. In other words, find
    if given gene is CancerMine and CancerMine_Cohort gene lists.
    E.g.
        ----+------------+-------------------+
        ..  | CancerMine | CancerMine_Cohort |
        ----+------------+-------------------+
        ..  |    YES     |         NO        |
        ..  |     NO     |         NO        |
        ..  |     ..     |         ..        |
        ----+------------+-------------------+


    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `CancerMine_STATUS` and `CancerMine_Cohort_STATUS` column will be to be added on.

        cancermine_genes : <list>
            CancerMine gene list.

        cancermine_cohort_genes : <list>
            CancerMine cohort gene list. Must be one of BRCA, COAD, or OV, depending on `preliminary_data`.

        tcga_type : <str>
            The TCGA abbreviation, used for column name. Must be one of BRCA, COAD, or OV.

    Returns
    -------
        None. Modifies the dataframe.
    """
    cancermine_statuses = []
    cancermine_cohort_statuses = []

    for gene in preliminary_data["GENE"]:
        cancermine_status = "+" if gene in cancermine_genes else "-"
        cancermine_cohort_status = "+" if gene in cancermine_cohort_genes else "-"
        cancermine_statuses.append(cancermine_status)
        cancermine_cohort_statuses.append(cancermine_cohort_status)

    # Add the columns
    preliminary_data['CancerMine_STATUS'] = cancermine_statuses
    preliminary_data[f'CancerMine_STATUS ({tcga_type})'] = cancermine_cohort_statuses
