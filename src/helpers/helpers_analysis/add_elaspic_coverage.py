from tqdm.notebook import tqdm
from .get_patient_protein_to_mutations_dict import get_patient_protein_to_mutations_dict


def add_elaspic_coverage(preliminary_data, elaspic_combined_data, snv_data):
    """
    Given preliminary_data, elaspic_combined_data and snv_data;
    
        Go over each patients in SNV data.
        List all specific protein.mutation pairs for current patient.
        If *any* of these mutations are in ELASPIC, increase the count (without considering CORE or INTERFACE)
        Then, go to the next patient.
        
    This can be thought of as another kind of baseline.
        
    Note: We count only once if there are multiple occurrance of the same protein's mutations.
    
    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `ELASPIC_COVERAGE` column will be to be added on.
            
        elaspic_combined_data : <dataframe>
            The ELASPIC combined results file which contains both CORE and INTERFACE.
        
        snv_data : <DataFrame>
            An SNV dataframe, we use the processed version of SNV.
    
    Returns
    -------
        None. Modifies the dataframe.
    """
    
    # Initialize the `proteins_to_counts_dict` dictionary with proteins in preliminary_data in corrent order, 
    # and the default values of 0.
    proteins_to_counts_dict = dict.fromkeys(list(preliminary_data['PROTEIN']), 0)
    
    # Get the patient IDs.
    patients = list(snv_data['Tumor_Sample_Barcode'].unique())
    
    for patient in tqdm(patients):
        
        # Patient filtered dataframe: Filter SNV file for current patient.
        patient_snv_data = snv_data[snv_data["Tumor_Sample_Barcode"] == patient]
        
        # Get the dictionary that maps protein <str> to mutations <list> for *current patient*.
        patient_protein_to_mutations = get_patient_protein_to_mutations_dict(patient_snv_data)      
        
        # Fill the counts dictionary `proteins_to_counts_dict`
        for protein, mutations in patient_protein_to_mutations.items():
            for mutation in mutations:
                
                # Search in ELASPIC.
                elaspic_search_data = elaspic_combined_data[(elaspic_combined_data["UniProt_ID"] == protein) & 
                                                            (elaspic_combined_data["Mutation"] == mutation)]
                
                # If search data is not empty, i.e. A specific Mutation of current Protein
                # exits in ELASPIC data, increase the count and skip to the next Protein.
                if not elaspic_search_data.empty:
                    proteins_to_counts_dict[protein] += 1
                    break

                    
    # Add the column
    preliminary_data['ELASPIC_COVERAGE'] = list(proteins_to_counts_dict.values())
    