def add_num_interface_patients_disruptive_interactor(
        preliminary_data, patient_interaction_data
):
    """
    Given a tcga preliminary_data and patient interaction data for that cohort,
    add number of interface patients where *at least one* interaction is disrupted
    for that gene.

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `NUM_INTERFACE_PATIENTS_DISRUPTIVE_INTERACTOR` column
            will be added on.

        patient_interaction_data : <DataFrame>
            The patient interaction data with disruptive interactors for each patient. Also, there is
            CORE / INTERFACE categorization for each patient, shown with "C" or "I".

    Returns
    -------
        None. Modifies the dataframe.
    """

    num_interface_patients_disruptive_interactor = []
    for protein in preliminary_data["PROTEIN"]:
        num_patients = patient_interaction_data[
            # Looking for Interface patients only -- with "I".
            (patient_interaction_data["CORE_INTERFACE_VS_INTERFACE_STATUS"] == "I") &
            # Looking for at least one disruptive interactors.
            (patient_interaction_data["NUM_DISRUPTIVE_INTERACTORS"] > 0) &
            # Looking for current protein.
            (patient_interaction_data["PROTEIN_GENE"].apply(lambda x: x.split(":")[0] == protein))
        ]["PATIENT"].nunique()

        num_interface_patients_disruptive_interactor.append(num_patients)

    preliminary_data['NUM_INTERFACE_PATIENTS_DISRUPTIVE_INTERACTOR'] = num_interface_patients_disruptive_interactor
