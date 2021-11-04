# Disruptive Interactions for Each Patients.

# Let's do it for BRCA

# read the prediction data

# read the SNV data and get `patient_to_mutations` dictionary.

# For each patient, obtain the disruptive interactions located in prediction.
# We look for protein.mutation only here.


# import os

# os.chdir('../../')
# import pandas as pd
# from pandas import DataFrame
# import re
# import time
#
# import requests
from collections import defaultdict
from datetime import datetime

from pandas import read_csv
import pandas as pd
from pathlib import Path
import os.path as op

from .loaders import load_snv_datasets

from tqdm.auto import tqdm
# import os.path as op

from ..mylogger import get_handler
import logging

from .gene_id_retrieval import GeneIDRetriever

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class InteractionsPerPatient:
    def __init__(
            self,
            tcga: str,
            prediction_data_path: Path,
            tcga_snv_path: Path,
            identifier: str,
            verbose: bool = False
    ):
        """
        Interactions per Patient.

        Parameters
        ----------
            tcga : <str>
                TCGA Cohort name, e.g. `BRCA`.

            prediction_data_path : <Path>
                The path of prediction data.

            tcga_snv_path : <Path>
                The path of original TCGA SNV data.

            identifier : <str>
                The identifier to be used in protein/gene representation.
                - `uniprot` for UniProt_ID protein representation.
                - `hugo` for Hugo ID Gene representation.

            verbose : <bool>
                Controls the detail progress to be displayed. E.g. decides whether patient in
                progress to be displayed when finding that patient's disruptive interactions.

        """
        self.tcga = tcga.upper()
        self.prediction_data_path = prediction_data_path
        self.tcga_snv_path = tcga_snv_path
        self.identifier = identifier
        self.verbose = verbose

        self.snv_data_simplified = None
        self.patients = None
        self.patient_to_snv_data = None
        self.disruptive_prediction_data = None
        self.prediction_data = None
        self.analysis_table = None

        self.patient_to_interactions = None
        self.patient_to_disruptive_interactions = None

        self.uniprot_to_gene_id = None
        self.gene_retriever = GeneIDRetriever()
        self.load_materials()

        self.find_all_interactions()
        self.find_disruptive_interactions()

        log.info("Disruptive interactions per patient completed.")

    def load_materials(self):
        log.info("Loading materials ..")
        self.load_snv_data_simplified()
        self.load_patient_ids()
        self.load_patient_to_snv_data()
        self.load_prediction_data()
        self.load_disruptive_prediction_data()
        self.load_uniprot_to_gene_id()
        log.info("Materials loaded.")
        log.info(f"Number of {self.tcga} patients: {len(self.patients)}.")

    def load_snv_data_simplified(self):
        log.debug("Loading SNV data simplified ..")
        data_materials = {}
        load_snv_datasets(self.tcga, self.tcga_snv_path, data_materials)
        self.snv_data_simplified = data_materials[f"{self.tcga}_snv_data_simplified"]

    def load_patient_ids(self):
        log.debug("Loading patient ids ..")
        patients = list(self.snv_data_simplified["Tumor_Sample_Barcode"].unique())
        self.patients = patients

    def get_patient_snv_data(self, patient_id):
        patient_snv_data = self.snv_data_simplified[
            self.snv_data_simplified["Tumor_Sample_Barcode"] == patient_id
            ]
        return patient_snv_data

    def load_patient_to_snv_data(self):
        log.debug("Loading patient to snv_data ..")
        patient_to_snv_data = {}
        for patient in tqdm(self.patients):
            patient_snv_data = self.get_patient_snv_data(patient)
            patient_to_snv_data[patient] = patient_snv_data

        self.patient_to_snv_data = patient_to_snv_data

    def load_prediction_data(self):
        log.debug("Loading the prediction data ..")
        prediction_data = read_csv(self.prediction_data_path)
        self.prediction_data = prediction_data

    def load_disruptive_prediction_data(self):
        log.debug("Loading the disruptive prediction data ..")
        disruptive_prediction_data = self.prediction_data[
            self.prediction_data["Prediction"] == 0
            ]
        self.disruptive_prediction_data = disruptive_prediction_data

    def get_disruptive_interactors(self, protein, mutation):
        disruptive_predicted_interactors = self.disruptive_prediction_data[
            (self.disruptive_prediction_data["UniProt_ID"] == protein) &
            (self.disruptive_prediction_data["Mutation"] == mutation) &
            (self.disruptive_prediction_data["Prediction"] == 0)  # Just to be double sure.
            ]["Interactor_UniProt_ID"].to_list()

        return disruptive_predicted_interactors

    def get_non_disruptive_interactors(self, protein, mutation):
        non_disruptive_predicted_interactors = self.prediction_data[
            (self.prediction_data["UniProt_ID"] == protein) &
            (self.prediction_data["Mutation"] == mutation) &
            (self.prediction_data["Prediction"] == 1)  # Class 1
            ]["Interactor_UniProt_ID"].to_list()

        return non_disruptive_predicted_interactors

    def get_all_interactors(self, protein, mutation):
        interactors = self.prediction_data[
            (self.prediction_data["UniProt_ID"] == protein) &
            (self.prediction_data["Mutation"] == mutation)
            ]["Interactor_UniProt_ID"].to_list()

        return interactors

    # def find_disruptive_interactions_single_patient(self, patient):
    #  # # # INCOMPLETE FUNCTION -- NOT NEEDED.
    #     log.info(f"Finding disruptive interactions for patient: {patient} ..")
    #     patient_snv_data = self.patient_to_snv_data[patient]
    #     disruptive_interactions = []
    #     for index, row in patient_snv_data.iterrows():
    #         protein = row["SWISSPROT"]
    #         mutation = row["HGVSp_Short"]
    #         disruptive_predicted_interactions = self.get_disruptive_predicted_interactions(
    #             protein, mutation
    #         )
    #         for interactor in disruptive_predicted_interactions:
    #             print(f"{protein}, {mutation}, {interactor}")

    def find_all_interactions(self):
        log.info("Finding all interactions (disruptive and non-disruptive) for each patient ..")
        patient_to_interactions = {}
        for patient in tqdm(self.patients):
            if self.verbose:
                log.debug(f"\tPATIENT: {patient}")
            interactions = []
            patient_snv_data = self.patient_to_snv_data[patient]

            for index, row in patient_snv_data.iterrows():
                protein = row["SWISSPROT"]
                mutation = row["HGVSp_Short"]
                # for current protein.mutation
                all_interactors = self.get_all_interactors(
                    protein, mutation
                )

                # Add triplets
                for interactor in all_interactors:
                    # log.debug(f"\t\interaction: {protein, mutation, interactor}")

                    interactions.append(
                        (protein, mutation, interactor)
                    )

            patient_to_interactions[patient] = interactions

        self.patient_to_interactions = patient_to_interactions

    def find_disruptive_interactions(self):
        log.info("Finding disruptive interactions for each patient ..")
        patient_to_disruptive_interactions = {}
        for patient in tqdm(self.patients):
            if self.verbose:
                log.debug(f"\tPATIENT: {patient}")
            disruptive_interactions = []
            patient_snv_data = self.patient_to_snv_data[patient]

            for index, row in patient_snv_data.iterrows():
                protein = row["SWISSPROT"]
                mutation = row["HGVSp_Short"]
                # for current protein.mutation
                disruptive_predicted_interactions = self.get_disruptive_interactors(
                    protein, mutation
                )

                # Add disruptive predicted triplets
                for interactor in disruptive_predicted_interactions:
                    # log.debug(f"\t\tDisruptive interaction: {protein, mutation, interactor}")

                    disruptive_interactions.append(
                        (protein, mutation, interactor)
                    )

            patient_to_disruptive_interactions[patient] = disruptive_interactions

        self.patient_to_disruptive_interactions = patient_to_disruptive_interactions

    def print_disruptive_interactions_per_patient(self):
        """
        Prints the disruptive interactions for each patient.
        """
        for patient in self.patients:
            if self.identifier == "uniprot":
                print(f"{patient} -> {self.patient_to_disruptive_interactions[patient]}")

            elif self.identifier == "hugo":
                protein, mutation, interactor = self.patient_to_disruptive_interactions[patient]
                disruptive_interaction_hugo = (
                    f"{protein}:{self.gene_retriever.fetch(protein)}",
                    mutation,
                    f"{interactor}:{self.gene_retriever.fetch(interactor)}",
                )

                print(f"{patient} -> {disruptive_interaction_hugo}")

    def load_uniprot_to_gene_id(self):
        log.info("Loading UniProt ID to Gene ID ..")
        uniprot_to_gene_id = {}

        prediction_data = self.prediction_data.copy()
        self_proteins = list(prediction_data["UniProt_ID"].unique())
        interactor_proteins = list(prediction_data["Interactor_UniProt_ID"].unique())
        proteins = sorted(set(self_proteins + interactor_proteins))

        for protein in tqdm(proteins, desc="Retrieving Gene IDs from UniProt API .. "):
            gene = self.gene_retriever.fetch(protein=protein)
            uniprot_to_gene_id[protein] = gene

        self.uniprot_to_gene_id = uniprot_to_gene_id

        log.info("`uniprot_to_gene_id` loaded. ")

    # def load_uniprot_to_gene_id(self):
    #     log.info("Loading UniProt ID to Gene ID ..")
    #     uniprot_to_gene_id = {}
    #
    #     prediction_data = self.prediction_data.copy()
    #     self_proteins = list(prediction_data["UniProt_ID"].unique())
    #     interactor_proteins = list(prediction_data["Interactor_UniProt_ID"].unique())
    #     proteins = sorted(set(self_proteins + interactor_proteins))
    #
    #     for protein in tqdm(proteins, desc="Retrieving Gene IDs from UniProt API .. "):
    #         gene = self.get_gene_id_from_uniprot(protein)
    #         uniprot_to_gene_id[protein] = gene
    #
    #     self.uniprot_to_gene_id = uniprot_to_gene_id
    #
    #     log.info("`uniprot_to_gene_id` loaded. ")

    # def get_gene_id_from_uniprot(self, uniprot_id):
    #     # log.debug("Retrieving sequence {} ...".format(uniprot_id))
    #     address = "http://www.uniprot.org/uniprot/{}.fasta".format(uniprot_id)
    #     n_attempt = 3
    #     attempt = 0
    #     while attempt < n_attempt:
    #         r = requests.get(address)
    #         if r.status_code == 200:
    #             gene = self.get_gene_from_fasta(r.text)
    #             return gene
    #
    #         attempt += 1
    #         log.warning(f"attempt: {attempt}")
    #         log.warning(f"status_code: {r.status_code}")
    #         time.sleep(1)
    #
    #     log.critical(f"COULD NOT RETRIEVE GENE: {attempt}")
    #     return "N/A"

    # @staticmethod
    # def get_gene_from_fasta(fasta_text):
    #     info_line = fasta_text.split('\n')[0]
    #
    #     # pattern = re.compile(r"GN=(.+)")
    #     pattern = re.compile(r"GN=(\S+)(\s)")
    #
    #     # Does not exists in UniProt server.
    #     if re.search(pattern, info_line) is None:
    #         return "N/A"
    #
    #     gene = re.search(pattern, info_line).group(1)
    #
    #     return gene

    def construct_analysis_table(self):
        log.info("Constructing the analysis table ..")
        patient_to_table_entry_set = defaultdict(list)  # Dictionary of list containing unique dictionaries.
        patient_to_seen_protein_mutation_pairs = defaultdict(list)

        for patient in tqdm(self.patients):
            patient_interactions = self.patient_to_interactions[patient]

            for protein, mutation, _ in patient_interactions:
                """
                # unnecessary duplication: but DOES NOT work fine.
                TCGA-A8-A093	P28062	R216W	['P40306']
                TCGA-A8-A093	Q15842	E237K	['Q14654', 'P63252']  
                TCGA-A8-A093	Q15842	E237K	['Q14654', 'P63252']  <- duplicated entries
                """
                # Skipping unnecessary duplication.
                if (protein, mutation) in patient_to_seen_protein_mutation_pairs[patient]:
                    continue

                interactors = self.get_all_interactors(protein, mutation)
                disruptive_interactors = self.get_disruptive_interactors(protein, mutation)
                non_disruptive_interactors = self.get_non_disruptive_interactors(protein, mutation)

                print(
                    f"{patient}\t{protein}\t{mutation}\tAll interactors: {interactors}"
                )

                patient_to_seen_protein_mutation_pairs[patient].append((protein, mutation))

                # Adding table info to table_values_set.
                patient_to_table_entry_set[patient].append(
                    {
                        "PATIENT": patient,

                        "PROTEIN_GENE": f"{protein}:{self.gene_retriever.fetch(protein)}",

                        "MUTATION": mutation,

                        # All Interactors
                        "INTERACTORS": ','.join(
                            map(str, map(lambda x: f"{x}:{self.gene_retriever.fetch(x)}", interactors))
                        ),
                        "NUM_INTERACTORS": len(interactors),

                        # Disruptive Interactors
                        "DISRUPTIVE_INTERACTORS": ','.join(
                            map(str, map(lambda x: f"{x}:{self.gene_retriever.fetch(x)}", disruptive_interactors))
                        ),
                        "NUM_DISRUPTIVE_INTERACTORS": len(disruptive_interactors),

                        # Non-Disruptive Interactors
                        "NON_DISRUPTIVE_INTERACTORS": ','.join(
                            map(str, map(lambda x: f"{x}:{self.gene_retriever.fetch(x)}", non_disruptive_interactors))
                        ),
                        "NUM_NON_DISRUPTIVE_INTERACTORS": len(non_disruptive_interactors),
                    }
                )

            print('-' * 100)

        analysis_table_entries = []
        for patient, table_entry_set in patient_to_table_entry_set.items():
            for table_entry_dict in table_entry_set:
                # Add table entry dictionary.
                analysis_table_entries.append(table_entry_dict)

        analysis_table = pd.DataFrame(
            analysis_table_entries,
            columns=[
                "PATIENT",
                "PROTEIN_GENE",
                "MUTATION",
                "INTERACTORS",
                "NUM_INTERACTORS",
                "DISRUPTIVE_INTERACTORS",
                "NUM_DISRUPTIVE_INTERACTORS",
                "NON_DISRUPTIVE_INTERACTORS",
                "NUM_NON_DISRUPTIVE_INTERACTORS",
            ]
        )

        self.analysis_table = analysis_table
        log.info("Analysis table constructed.")

    def extract(self, folder=None):
        output_file_date = datetime.today().strftime('%Y-%m-%d')
        filename = f"{self.tcga}_patient_interactions_analysis_table_{output_file_date}.xlsx"

        if folder is not None:
            filename = op.join(folder, filename)

        # Ensure the file is not exists before creating to prevent overwriting.
        if op.isfile(filename):
            log.warning(f'File {filename} is already exist.')

        else:
            # Export
            self.analysis_table.to_excel(filename, index=False)
            log.info(f'{filename} is exported.')
