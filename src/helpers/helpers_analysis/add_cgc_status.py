def add_cgc_status(preliminary_data, cgc_genes, cgc_cohort_genes, tcga_type):
    """
    Given a tcga preliminary_data data, add the CGC and CGC_Cohort status. In other words, find 
    if given gene is CGC and CGC_Cohort gene lists.
    E.g.
        ----+-----+------------+
        ..  | CGC | CGC_Cohort |
        ----+-----+------------+
        ..  | YES |     NO     |
        ..  | NO  |     NO     |
        ..  | ..  |     ..     |
        ----+-----+------------+


    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `CGC_STATUS` and `CGC_Cohort_STATUS` column will be to be added on.

        cgc_genes : <list>
            CGC gene list.

        cgc_cohort_genes : <list>
            CGC cohort gene list. Must be one of BRCA, COAD, or OV, depending on `preliminary_data`.

        tcga_type : <str>
            The TCGA abbreviation, used for column name. Must be one of BRCA, COAD, or OV.

    Returns
    -------
        None. Modifies the dataframe.
    """
    cgc_statuses = []
    cgc_cohort_statuses = []

    for gene in preliminary_data["GENE"]:
        cgc_status = "+" if gene in cgc_genes else "-"
        cgc_cohort_status = "+" if gene in cgc_cohort_genes else "-"
        cgc_statuses.append(cgc_status)
        cgc_cohort_statuses.append(cgc_cohort_status)

    # Add the columns
    preliminary_data['CGC_STATUS'] = cgc_statuses
    preliminary_data[f'CGC_STATUS ({tcga_type})'] = cgc_cohort_statuses
