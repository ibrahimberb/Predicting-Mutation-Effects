# Disruptive Interactions for Each Patients.

# Let's do it for BRCA

# read the prediction data

# read the SNV data and get `patient_to_mutations` dictionary.

# For each patient, obtain the disruptive interactions located in prediction.
# We look for protein.mutation only here.


import os

# os.chdir('../../')
import pandas as pd
from pandas import DataFrame
from pandas import read_csv
from pathlib import Path

from .loaders import load_snv_datasets

from ..mylogger import get_handler
import logging

from tqdm.auto import tqdm

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class DisruptiveInteractionsPerPatient:
    def __init__(
            self,
            tcga: str,
            prediction_data_path: Path,
            tcga_snv_path: Path
    ):
        self.tcga = tcga.upper()
        self.prediction_data_path = prediction_data_path
        self.tcga_snv_path = tcga_snv_path

        self.snv_data_simplified = None
        self.patients = None
        self.patient_to_snv_data = None
        self.disruptive_prediction_data = None
        self.prediction_data = None
        self.patient_to_disruptive_interactions = None

        self.load_materials()
        self.find_disruptive_interactions()

    def load_materials(self):
        log.info("Loading materials ..")
        self.load_snv_data_simplified()
        self.load_patient_ids()
        self.load_patient_to_snv_data()
        self.load_prediction_data()
        self.load_disruptive_prediction_data()
        log.info("Materials loaded.")

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

    def get_disruptive_predicted_interactions(self, protein, mutation):
        disruptive_predicted_interactions = self.disruptive_prediction_data[
            (self.disruptive_prediction_data["UniProt_ID"] == protein) &
            (self.disruptive_prediction_data["Mutation"] == mutation) &
            (self.disruptive_prediction_data["Prediction"] == 0)  # Just to be double sure.
        ]["Interactor_UniProt_ID"].to_list()

        return disruptive_predicted_interactions

    def find_disruptive_interactions(self):
        log.info("Finding disruptive interactions for each patient ..")
        patient_to_disruptive_interactions = {}
        for patient in tqdm(self.patients):
            log.debug(f"\tPATIENT: {patient}")
            disruptive_interactions = []
            patient_snv_data = self.patient_to_snv_data[patient]

            for index, row in patient_snv_data.iterrows():
                protein = row["SWISSPROT"]
                mutation = row["HGVSp_Short"]
                # for current protein.mutation
                disruptive_predicted_interactions = self.get_disruptive_predicted_interactions(
                    protein, mutation
                )

                # Add disruptive predicted triplets
                for interactor in disruptive_predicted_interactions:
                    log.debug(f"\t\tDisruptive interaction: {protein, mutation, interactor}")
                    disruptive_interactions.append(
                        (protein, mutation, interactor)
                    )

            patient_to_disruptive_interactions[patient] = disruptive_interactions

        self.patient_to_disruptive_interactions = patient_to_disruptive_interactions
        # todo: <-- remained here.