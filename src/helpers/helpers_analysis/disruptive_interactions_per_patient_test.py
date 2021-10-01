import os
from pathlib import Path
from unittest import TestCase
import re

from .disruptive_interactions_per_patient import DisruptiveInteractionsPerPatient


class TestDisruptiveInteractionsPerPatient(TestCase):
    os.chdir("../")

    PREDICTION_ID = "acf35ed1/"
    PREDICTIONS_COMMON_PATH = "../data/predictions_datasets/brca_prediction_2021-09-28/" + PREDICTION_ID
    PREDICTION_BRCA_REDUCED_PATH = PREDICTIONS_COMMON_PATH + "predictions_soft_2021-09-28.csv"
    PREDICTION_BRCA_REDUCED_PATH = Path(PREDICTION_BRCA_REDUCED_PATH)

    SNV_COMMON_PATH = "C:/Users/ibrah/Desktop/TUSEB_Study/Data_Collection_and_Filtering/SNV/"
    SNV_PATH = os.path.join(SNV_COMMON_PATH, "SNV_BRCA_hg38.csv")
    SNV_PATH = Path(SNV_PATH)

    def setUp(self) -> None:
        self.disruptive_interactions_per_patient = DisruptiveInteractionsPerPatient(
            tcga="brca",
            prediction_data_path=self.PREDICTION_BRCA_REDUCED_PATH,
            tcga_snv_path=self.SNV_PATH
        )

    def test_load_snv_data_simplified(self):
        snv_data_simplified = self.disruptive_interactions_per_patient.snv_data_simplified
        columns_list = ["Hugo_Symbol", "SWISSPROT", "HGVSp_Short", "Tumor_Sample_Barcode"]
        self.assertEqual(snv_data_simplified.columns.to_list(), columns_list)

    def test_get_patient_snv_data(self):
        patients = self.disruptive_interactions_per_patient.patients
        for patient in patients:
            patient_snv_data = self.disruptive_interactions_per_patient.get_patient_snv_data(
                patient
            )

            self.assertEqual(patient_snv_data["Tumor_Sample_Barcode"].nunique(), 1)

            self.assertEqual(
                list(patient_snv_data["Tumor_Sample_Barcode"].unique())[0], patient
            )

    def test_patient_ids_are_valid(self):
        patient_ids = self.disruptive_interactions_per_patient.patients
        match = re.compile(r"^TCGA-(\w\w)-(\w\w\w\w)$")
        for patient in patient_ids:
            self.assertIsNotNone(
                match.match(patient), "Test value is not none."
            )

    def test_load_patient_to_snv_data(self):
        patient_to_snv_data = self.disruptive_interactions_per_patient.patient_to_snv_data
        patients = self.disruptive_interactions_per_patient.patients

        self.assertEqual(
            type(patient_to_snv_data), dict
        )

        self.assertEqual(
            sorted(patient_to_snv_data.keys()), sorted(patients)
        )

    def test_load_prediction_data(self):
        prediction_data = self.disruptive_interactions_per_patient.prediction_data

        self.assertEqual(prediction_data.empty, False)

        self.assertEqual(
            set(prediction_data["Prediction"].value_counts().index),
            {1, 0}
        )

    def test_load_disruptive_prediction_data(self):
        disruptive_prediction_data = self.disruptive_interactions_per_patient.disruptive_prediction_data

        self.assertEqual(
            set(disruptive_prediction_data["Prediction"].value_counts().index),
            {0}
        )

    def test_get_disruptive_predicted_interactions(self):

        if self.disruptive_interactions_per_patient.tcga == "BRCA":

            self.assertNotEqual(
                self.disruptive_interactions_per_patient.get_disruptive_predicted_interactions(
                    "Q01196", "G95R"
                ), []
            )

            self.assertNotEqual(
                self.disruptive_interactions_per_patient.get_disruptive_predicted_interactions(
                    "P00747", "D665H"
                ), []
            )

        else:
            raise ValueError("TCGA falls outside test cases.")

        # A fail case.
        self.assertEqual(
            self.disruptive_interactions_per_patient.get_disruptive_predicted_interactions(
                "P123123", "MUT123"
            ), []
        )

    def test_patient_to_disruptive_interactions(self):
        patient_to_disruptive_interactions = self.disruptive_interactions_per_patient.patient_to_disruptive_interactions

        if self.disruptive_interactions_per_patient.tcga == "BRCA":
            self.assertEqual(
                len(patient_to_disruptive_interactions), 985
            )

        else:
            raise ValueError("TCGA falls outside test cases.")

    def test_find_disruptive_interactions(self):

        patient_to_disruptive_interactions = self.disruptive_interactions_per_patient.patient_to_disruptive_interactions
        prediction_data = self.disruptive_interactions_per_patient.prediction_data
        disruptive_prediction_data = self.disruptive_interactions_per_patient.disruptive_prediction_data

        for patient, disruptive_triplets in patient_to_disruptive_interactions.items():
            patient_snv_data = self.disruptive_interactions_per_patient.get_patient_snv_data(patient)
            for disruptive_triplet in disruptive_triplets:
                protein, mutation, interactor = disruptive_triplet

                patient_search_data = patient_snv_data[
                    (patient_snv_data["SWISSPROT"] == protein) &
                    (patient_snv_data["HGVSp_Short"] == mutation)
                ]

                self.assertFalse(patient_search_data.empty, msg="data is not empty")

                prediction_search_data = prediction_data[
                    (prediction_data["UniProt_ID"] == protein) &
                    (prediction_data["Mutation"] == mutation) &
                    (prediction_data["Interactor_UniProt_ID"] == interactor)
                ]

                # Assert there is only one row.
                self.assertEqual(len(prediction_search_data), 1)

                # Ensure there is only one value in `Prediction` column and get it.
                self.assertEqual(prediction_search_data["Prediction"].nunique(), 1)
                [val] = prediction_search_data["Prediction"].values

                # Assert the prediction is "disruptive" (class 0)
                self.assertEqual(val, 0)

                disruptive_prediction_search_data = disruptive_prediction_data[
                    (disruptive_prediction_data["UniProt_ID"] == protein) &
                    (disruptive_prediction_data["Mutation"] == mutation) &
                    (disruptive_prediction_data["Interactor_UniProt_ID"] == interactor)
                ]

                # Assert there is only one row.
                self.assertEqual(len(disruptive_prediction_search_data), 1)

                # Ensure there is only one value in `Prediction` column and get it.
                self.assertEqual(disruptive_prediction_search_data["Prediction"].nunique(), 1)
                [val] = disruptive_prediction_search_data["Prediction"].values

                # Assert the prediction is "disruptive" (class 0)
                self.assertEqual(val, 0)

                # Assert the search datasets are identical.
                self.assertTrue(
                    prediction_search_data.equals(disruptive_prediction_search_data),
                    "datasets are not identical."
                )
