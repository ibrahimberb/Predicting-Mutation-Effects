import pickle
import pandas as pd
from IPython.display import display

with open("brca_patient_protein_mutation_to_C_I_status.pickle", 'rb') as inp:
    brca_patient_protein_mutation_to_C_I_status = pickle.load(inp)

# print(brca_patient_protein_mutation_to_C_I_status)

data = pd.DataFrame(
    brca_patient_protein_mutation_to_C_I_status.keys(),
    columns=["PROTEIN", "PATIENT"]
)
data["C_I_STATUS"] = brca_patient_protein_mutation_to_C_I_status.values()
display(data.head())

