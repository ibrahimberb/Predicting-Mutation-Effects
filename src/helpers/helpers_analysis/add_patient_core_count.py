from tqdm.notebook import tqdm
from .get_patient_protein_to_mutations_dict import get_patient_protein_to_mutations_dict
from .is_core import is_core
from .is_in_elaspic import is_in_elaspic


def add_patient_core_count(preliminary_data, snv_data, elaspic_core_data, elaspic_interface_data):
    """
    Given a tcga preliminary_data data, add number of patients that core mutation occurred
    for each protein as a new column.

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `PATIENT_CORE_COUNT` column will be to be added on.

        snv_data : <DataFrame>
            An SNV dataframe, we use the processed version of SNV.

        elaspic_core_data : <DataFrame>
            The ELASPIC results file that contains only the `core` type entries.

        elaspic_interface_data : <DataFrame>
            The ELASPIC results file that contains only the `interface` type entries. It will be used to
            check if a specific (protein, mutation) pair is an interface via `is_interface` function.

    Returns
    -------
        None. Modifies the dataframe.
    """

    # Initialize the `proteins_to_patient_core_counts_dict` dictionary with proteins in preliminary_data in
    # correct order, and the default values of 0.
    proteins_to_patient_core_counts_dict = dict.fromkeys(list(preliminary_data['PROTEIN']), 0)

    # Get the patient IDs.
    patients = list(snv_data['Tumor_Sample_Barcode'].unique())

    for patient in tqdm(patients):

        # Patient filtered dataframe: Filter SNV file for current patient.
        patient_snv_data = snv_data[snv_data["Tumor_Sample_Barcode"] == patient]

        for protein, mutations in get_patient_protein_to_mutations_dict(patient_snv_data).items():

            core_flag = 'N/A'
            # print(protein, mutations)
            for mutation in mutations:

                # Check if (protein.mutation) is in ELASPIC.
                if is_in_elaspic(protein, mutation, elaspic_core_data, elaspic_interface_data):
                    # print(f'{protein}.{mutation} IS IN ELASPIC.')

                    if is_core(protein, mutation, elaspic_core_data):
                        # print(' → core found!')
                        core_flag = 1
                        break

                    else:
                        # print(' → interface found!')
                        core_flag = 0

                else:
                    # print(f'{protein}.{mutation} IS NOT IN ELASPIC.')
                    # print(f'CORE_FLAG = {core_flag}')
                    continue

            if core_flag == 1:
                # Adding the corresponding gene counter
                proteins_to_patient_core_counts_dict[protein] += 1

    # Add the column
    preliminary_data['PATIENT_CORE_COUNT'] = list(proteins_to_patient_core_counts_dict.values())
