from collections import Counter

from .gene_id_retrieval import GeneIDFetcher
from .loaders import load_snv_datasets
from ..mylogger import get_handler
import logging
from tqdm.auto import tqdm
import pandas as pd
from IPython.display import display

import os.path as op
from datetime import datetime

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class MutualExclusivity:
    def __init__(
            self,
            tcga,
            tcga_snv_path,
            patient_interaction_data_path,
    ):
        self.tcga = tcga.upper()
        self.tcga_snv_path = tcga_snv_path
        self.patient_interaction_data_path = patient_interaction_data_path
        self.snv_data = None
        self.patients = None
        self.patient_to_snv_data = None
        self.patient_interaction_data = None

        self.load_snv_data()
        self.load_patient_ids()
        self.load_patient_to_snv_data()
        self.load_patient_interaction_data()

    def load_snv_data(self):
        log.debug("Loading SNV data simplified ..")
        data_materials = {}
        load_snv_datasets(self.tcga, self.tcga_snv_path, data_materials)
        self.snv_data = data_materials[f"{self.tcga}_snv_data_simplified"]

    def load_patient_ids(self):
        log.debug("Loading patient ids ..")
        patients = list(self.snv_data["Tumor_Sample_Barcode"].unique())
        self.patients = patients

    def get_patient_snv_data(self, patient_id):
        patient_snv_data = self.snv_data[
            self.snv_data["Tumor_Sample_Barcode"] == patient_id
            ]
        return patient_snv_data

    def load_patient_to_snv_data(self):
        log.debug("Loading patient to snv_data ..")
        patient_to_snv_data = {}
        for patient in tqdm(self.patients):
            patient_snv_data = self.get_patient_snv_data(patient)
            patient_to_snv_data[patient] = patient_snv_data

        self.patient_to_snv_data = patient_to_snv_data

    def load_patient_interaction_data(self):
        log.debug("patient interaction data patient data ..")
        patient_interaction_data = pd.read_excel(self.patient_interaction_data_path)
        self.patient_interaction_data = patient_interaction_data

    def get_patients_with(self, identifier, identifier_type):

        if identifier_type == "protein":
            patients = self.snv_data[
                self.snv_data["SWISSPROT"] == identifier
                ]["Tumor_Sample_Barcode"].unique()

            return list(patients)

        elif identifier_type == "gene":
            raise NotImplementedError

        else:
            raise ValueError(f"Invalid argument `{identifier_type}` for `identifier_type`.")

    def calculate_mutual_exclusivity(self, identifier_1, identifier_2, verbose=False, return_num_patients=False):
        """
                             | S1 union S2 |
        Mutual exclusivity = ---------------
                             |S1|  +  |S2|

        """
        # todo: using genes to SNV data may not be a good idea
        #  an alternative is using Protein, and retrieving gene id of that protein, and then querying in SNV data.
        s1 = set(self.get_patients_with(identifier_1, "protein"))
        s2 = set(self.get_patients_with(identifier_2, "protein"))

        if verbose:
            print(f"S1: {s1} ({len(s1)} patients)")
            print(f"S2: {s2} ({len(s2)} patients)")

        mutual_exclusivity_value = len(s1.union(s2)) / (len(s1) + len(s2))

        if return_num_patients:
            return mutual_exclusivity_value, len(s1), len(s2)

        else:
            return mutual_exclusivity_value

    def get_disrupted_interactors(self, identifier, identifier_type, return_counter=False):
        if identifier_type == "protein":
            identifier_pos = 0
        elif identifier_type == "gene":
            identifier_pos = 1
        else:
            raise ValueError(f"Invalid argument `{identifier_type}` for `identifier_type`.")

        patient_data = self.patient_interaction_data.copy()

        identifier_disruptive_interactors_series = patient_data[
            (patient_data["PROTEIN_GENE"].apply(lambda x: x.split(':')[identifier_pos]) == identifier)
        ]["DISRUPTIVE_INTERACTORS"]

        # Drop nan entries
        identifier_disruptive_interactors_series = identifier_disruptive_interactors_series.dropna()

        # Explode values
        identifier_disruptive_interactors_series = identifier_disruptive_interactors_series.apply(
            lambda x: x.split(",")
        ).explode()

        # get rid of probability value.
        identifier_disruptive_interactors = list(
            identifier_disruptive_interactors_series.apply(lambda x: ":".join(x.split(":")[:-1]))
        )

        identifier_disruptive_interactors_counter = Counter(identifier_disruptive_interactors)

        identifier_disruptive_interactors_unique = [
            f"{interactor.split(':')[0]}:{interactor.split(':')[1]}"
            for interactor, _ in identifier_disruptive_interactors_counter.most_common()
        ]

        # # get rid of probability value.
        # identifier_disruptive_interactors_series = identifier_disruptive_interactors_series.apply(
        #     lambda x: f"{x.split(':')[0]}:{x.split(':')[1]}"
        # )

        # identifier_disruptive_interactors = list(
        #     identifier_disruptive_interactors_series.apply(lambda x: x.split(',')).explode())
        # identifier_disruptive_interactors_counter = Counter(identifier_disruptive_interactors)

        # identifier_disruptive_interactors_unique = [
        #     f"{interactor.split(':')[0]}:{interactor.split(':')[1]}"
        #     for interactor, _ in identifier_disruptive_interactors_counter.most_common()
        # ]

        if return_counter:
            return identifier_disruptive_interactors_counter

        return identifier_disruptive_interactors_unique

    def get_disruptive_mutual_exclusivity_data(self, protein):
        gene = GeneIDFetcher().fetch(protein=protein)
        identifier_1_disruptive_interactors = self.get_disrupted_interactors(protein, identifier_type="protein")
        log.info(f"Calculating Mutual Exclusivity over {protein}'s interactors ..")
        log.debug(f"{protein} have {len(identifier_1_disruptive_interactors)} interactors:\n"
                  f"{identifier_1_disruptive_interactors}")

        mutual_exclusivity_entries = []
        for interactor in identifier_1_disruptive_interactors:
            interactor_protein, interactor_gene = interactor.split(":")
            mut_ex_value, len_s1, len_s2 = self.calculate_mutual_exclusivity(
                protein, interactor_protein, return_num_patients=True
            )
            # log.debug(f"{interactor} \t {mut_ex_value}")
            mutual_exclusivity_entries.append(
                (f"{protein}:{gene}", len_s1, interactor, len_s2, round(mut_ex_value, 4))
            )

        mutual_exclusivity_data = pd.DataFrame(
            mutual_exclusivity_entries,
            columns=["PROTEIN:GENE", "NUM_PATIENTS", "INTERACTOR", "NUM_PATIENTS_INTERACTOR", "MUTUAL_EXCLUSIVITY"]
        )

        return mutual_exclusivity_data

    def export_disruptive_mutual_exclusivity_data(self, folder_path, protein, file_extension="csv"):
        # Get Mut Ex data for given identifier.
        mutual_exclusivity_data = self.get_disruptive_mutual_exclusivity_data(protein)

        gene = GeneIDFetcher().fetch(protein=protein)

        log.debug(f"Exporting Mutual Exclusivity {self.tcga} {protein} ..")
        file_name = f"{self.tcga}_{protein}_{gene}"
        file_name = op.join(folder_path, file_name)
        file_date = datetime.today().strftime('%Y-%m-%d')
        file_name = f'{file_name}_{file_date}.{file_extension}'

        # Ensure the file is not exists before creating to prevent overwriting.
        if op.isfile(file_name):
            log.warning(f"File {file_name} is already exist.\n"
                        "To overwrite existing file, use `overwrite=True`.")
        else:
            # Export
            mutual_exclusivity_data.to_csv(file_name, index=False)
            log.info(f'{file_name} is exported successfully.')
