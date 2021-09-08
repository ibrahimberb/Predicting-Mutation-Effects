from pathlib import Path
from typing import Tuple, List
from IPython.display import display

import pandas as pd
import numpy as np
from sklearn import metrics

from helpers.helpers_analysis.loaders import (
    load_snv_datasets,
    load_prediction_dataset,
    load_elaspic_datasets,
    load_reference_dataset,
)

from helpers.mylogger import get_handler
import logging

from helpers.helpers_analysis.get_elaspic_proteins import get_elaspic_proteins
from helpers.helpers_analysis.get_protein_to_gene_dict import get_protein_to_gene_dict
from helpers.helpers_analysis.get_protein_to_num_elaspic_entries_dict import \
    get_protein_to_num_elaspic_interface_entires_dict
from helpers.helpers_analysis.get_protein_to_num_unique_interactors import get_protein_to_num_unique_interactors
from helpers.helpers_analysis.add_baseline import add_baseline
from helpers.helpers_analysis.add_elaspic_coverage import add_elaspic_coverage
from helpers.helpers_analysis.add_genes import add_genes
from helpers.helpers_analysis.add_num_disruptive_entries import add_num_disruptive_entries
from helpers.helpers_analysis.add_num_elaspic_interface_entries import add_num_elaspic_interface_entries
from helpers.helpers_analysis.add_num_incr_noeff_entries import add_num_incr_noeff_entries
from helpers.helpers_analysis.add_num_unique_interactors import add_num_unique_interactors
from helpers.helpers_analysis.add_our_method import add_our_method
from helpers.helpers_analysis.add_patient_core_count import add_patient_core_count
from helpers.helpers_analysis.add_patient_interface_count import add_patient_interface_count
from helpers.helpers_analysis.counts_baseline_vs_our_method import counts_baseline_vs_our_method
from helpers.helpers_analysis.add_cancermine_status import add_cancermine_status

from helpers.helpers_analysis.plot_roc_curve import roc_curve_analysis

from helpers.helpers_analysis.loaders import ReferenceDataset
from helpers.helpers_analysis.get_fpr_tpr_ths import get_fpr_tpr_ths

import matplotlib.pyplot as plt
import seaborn as sns

from helpers.helpers_analysis.add_cgc_status import add_cgc_status

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)

TCGA_CODE = str
REFERENCE_DATA_NAME = str
CohortSpecificReferenceDataPath, ReferenceDataPath = Path, Path

SNV_COMMON_PATH = "C:/Users/ibrah/Desktop/TUSEB_Study/Data_Collection_and_Filtering/SNV/"
SNV_BRCA_PATH = Path(SNV_COMMON_PATH + "SNV_BRCA_hg38.csv")
SNV_COAD_PATH = Path(SNV_COMMON_PATH + "SNV_COAD_hg38.csv")
SNV_OV_PATH = Path(SNV_COMMON_PATH + "SNV_OV_hg38.csv")

## PREDICTIONS_DATASETS
PREDICTIONS_COMMON_PATH = "predictions_datasets/ML_6_shap_SPSM/"
# PREDICTION_BRCA_PATH = PREDICTIONS_COMMON_PATH + "___.csv"
### PREDICTION_REDUCED_DATASETS
PREDICTION_BRCA_REDUCED_PATH = Path(PREDICTIONS_COMMON_PATH + "brca_predicted_reduced_data_shap_top-11_2021-07-27.csv")

## ELASPIC_RESULTS_DATASETS
ELASPIC_RESULTS_COMMON_PATH = "elaspic_results_datasets/"
BRCA_CORE_PATH = Path(ELASPIC_RESULTS_COMMON_PATH + "BRCA_CORE_A2.txt")
BRCA_INTERFACE_PATH = Path(ELASPIC_RESULTS_COMMON_PATH + "BRCA_INTERFACE_A2.txt")

# CANCER MINE GENES
CANCER_MINE_ALL_PATH = Path("cancer_mine_genes/all_genes_2021-07-18.txt")
CANCER_MINE_BREAST_PATH = Path("cancer_mine_genes/breast_genes_2021-07-18.txt")


class PredatorAnalysis:
    def __init__(
            self,
            tcga: TCGA_CODE,
            snv_path: Path,
            prediction_data_path: Path,
            elaspic_core_path: Path,
            elaspic_interface_path: Path,
            reference_data_name: REFERENCE_DATA_NAME,
            reference_data_spec_cohort_path: Path,
            reference_data_path: Path
    ):
        self.tcga = tcga
        self.snv_path = snv_path
        self.prediction_data_path = prediction_data_path
        self.elaspic_core_path = elaspic_core_path
        self.elaspic_interface_path = elaspic_interface_path
        self.reference_data_name = reference_data_name
        self.reference_data_spec_cohort_path = reference_data_spec_cohort_path
        self.reference_data_path = reference_data_path
        self.data_materials = {}

        # setattr(self, self.tcga, {})

        self.load_datasets()

    def load_datasets(self):
        load_snv_datasets(self.tcga, self.snv_path, self.data_materials)
        load_prediction_dataset(self.tcga, self.prediction_data_path, self.data_materials)
        load_elaspic_datasets(
            self.tcga, self.elaspic_core_path, self.elaspic_interface_path, self.data_materials
        )
        load_reference_dataset(
            self.tcga,
            self.reference_data_name,
            self.reference_data_spec_cohort_path,
            self.reference_data_path,
            self.data_materials
        )

    def prepare_analysis(self):
        # Proteins
        elaspic_proteins = get_elaspic_proteins(
            self.data_materials[f"{self.tcga}_elaspic_core_data"],
            self.data_materials[f"{self.tcga}_elaspic_interface_data"]
        )
        self.data_materials[f"{self.tcga}_elaspic_proteins"] = elaspic_proteins
        log.debug(f"{self.tcga}_elaspic_proteins loaded.")
        num_proteins = len(self.data_materials[f"{self.tcga}_elaspic_proteins"])
        log.debug(f"Number of proteins in ELASPIC {self.tcga}: {num_proteins}")

        # Genes
        protein_to_gene_dict = get_protein_to_gene_dict(
            self.data_materials[f"{self.tcga}_elaspic_proteins"],
            self.data_materials[f"{self.tcga}_snv_data_simplified"]
        )
        self.data_materials[f"{self.tcga}_protein_to_gene_dict"] = protein_to_gene_dict
        log.debug(f"{self.tcga}_protein_to_gene_dict loaded.")

        # ELASPIC Number of Interface Entries
        protein_to_num_elaspic_interface_entries = get_protein_to_num_elaspic_interface_entires_dict(
            self.data_materials[f"{self.tcga}_elaspic_proteins"],
            self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"]
        )
        self.data_materials[f"{self.tcga}_protein_to_num_elaspic_interface_entries"] = protein_to_num_elaspic_interface_entries
        log.debug(f"{self.tcga}_protein_to_num_elaspic_interface_entries loaded.")

        # ELASPIC Number of Unique Interactors
        protein_to_num_unique_interactors = get_protein_to_num_unique_interactors(
            self.data_materials[f"{self.tcga}_elaspic_proteins"],
            self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"]
        )
        self.data_materials[f"{self.tcga}_protein_to_num_unique_interactors"] = protein_to_num_unique_interactors
        log.debug(f"{self.tcga}_protein_to_num_unique_interactors loaded.")

        # Patients
        patients = list(
            self.data_materials[f"{self.tcga}_snv_data_simplified"]['Tumor_Sample_Barcode'].unique()
        )
        log.debug(f'Number of patients in {self.tcga}: {len(patients)}.')
        self.data_materials[f"{self.tcga}_patients"] = patients

    def construct_analysis_table(self):
        # 1. Adding `PROTEIN` Column
        log.debug(f"Adding `PROTEIN` column ..")
        preliminary_data = pd.DataFrame(
            self.data_materials[f"{self.tcga}_elaspic_proteins"], columns=['PROTEIN']
        )

        # 2. Adding `GENE` Column
        log.debug(f"Adding `GENE` column ..")
        add_genes(preliminary_data, self.data_materials[f"{self.tcga}_protein_to_gene_dict"])

        # 3. Adding `NUM_ELASPIC_INTERFACE_ENTRIES` Column
        log.debug(f"Adding `NUM_ELASPIC_INTERFACE_ENTRIES` column ..")
        add_num_elaspic_interface_entries(
            preliminary_data, self.data_materials[f"{self.tcga}_protein_to_num_elaspic_interface_entries"]
        )

        # 4. Adding `NUM_DISRUPTIVE_ENTRIES` Column
        log.debug(f"Adding `NUM_DISRUPTIVE_ENTRIES` column ..")
        add_num_disruptive_entries(
            preliminary_data, self.data_materials[f"{self.tcga}_prediction_data"]
        )

        # 5. Adding `NUM_INCR_NOEFF_ENTRIES` Column
        log.debug(f"Adding `NUM_INCR_NOEFF_ENTRIES` column ..")
        add_num_incr_noeff_entries(
            preliminary_data, self.data_materials[f"{self.tcga}_prediction_data"]
        )

        # 6. Adding `NUM_UNIQUE_INTERACTORS` Column
        log.debug(f"Adding `NUM_UNIQUE_INTERACTORS` column ..")
        add_num_unique_interactors(
            preliminary_data, self.data_materials[f"{self.tcga}_protein_to_num_unique_interactors"]
        )

        # 7. Adding `PATIENT_CORE_COUNT` Column
        log.debug(f"Adding `PATIENT_CORE_COUNT` column ..")
        add_patient_core_count(
            preliminary_data,
            self.data_materials[f"{self.tcga}_snv_data_simplified"],
            elaspic_interface_data=self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"],
            elaspic_core_data=self.data_materials[f"{self.tcga}_elaspic_core_data"]
        )

        # 8. Adding `PATIENT_INTERFACE_COUNT` Column
        log.debug(f"Adding `PATIENT_INTERFACE_COUNT` column ..")
        add_patient_interface_count(
            preliminary_data,
            self.data_materials[f"{self.tcga}_snv_data_simplified"],
            elaspic_interface_data=self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"],
            elaspic_core_data=self.data_materials[f"{self.tcga}_elaspic_core_data"]
        )

        # 9. Adding `BASELINE` and `OUR_METHOD` Columns
        log.debug(f"Adding `BASELINE` and `OUR_METHOD` columns ..")
        proteins_to_counts_baseline_dict, proteins_to_counts_our_method_dict = counts_baseline_vs_our_method(
            proteins=self.data_materials[f"{self.tcga}_elaspic_proteins"],
            patients=self.data_materials[f"{self.tcga}_patients"],
            snv_data=self.data_materials[f"{self.tcga}_snv_data_simplified"],
            elaspic_core_data=self.data_materials[f"{self.tcga}_elaspic_core_data"],
            elaspic_interface_data=self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"],
            prediction_data=self.data_materials[f"{self.tcga}_prediction_data"],
            add_core_flag_1_case_dict=None
        )
        add_baseline(preliminary_data, proteins_to_counts_baseline_dict)
        add_our_method(preliminary_data, proteins_to_counts_our_method_dict)

        # 10. Adding `OUR_METHOD / BASELINE` Column
        log.debug(f"Adding `OUR_METHOD / BASELINE` column ..")
        preliminary_data["OUR_METHOD/BASELINE"] = preliminary_data["OUR_METHOD"] / preliminary_data["BASELINE"]

        # 11. ELASPIC_COVERAGE
        log.debug(f"Adding `ELASPIC_COVERAGE` column ..")
        add_elaspic_coverage(
            preliminary_data,
            self.data_materials[f"{self.tcga}_elaspic_core_and_interface_data"],
            self.data_materials[f"{self.tcga}_snv_data_simplified"]
        )

        # 12. Adding Reference Dataset Columns: General and Cohort Specific
        # TODO Code refactor: add_reference_data_status
        log.debug(f"Adding Reference Dataset Columns: General and Cohort Specific columns ..")
        if self.reference_data_name == ReferenceDataset.CANCERMINE:
            add_cancermine_status(
                preliminary_data,
                cancermine_genes=self.data_materials["cancermine_all_genes"],
                cancermine_cohort_genes=self.data_materials[f"cancermine_{self.tcga}_genes"],
                tcga_type=self.tcga.upper()
            )

        elif self.reference_data_name == ReferenceDataset.CGC:
            add_cgc_status(
                preliminary_data,
                cgc_genes=self.data_materials["cgc_all_genes"],
                cgc_cohort_genes=self.data_materials[f"cgc_{self.tcga}_genes"],
                tcga_type=self.tcga.upper()
            )

        else:
            raise ValueError("Invalid Reference Dataset.")

        # Finally, assign the attribute.
        self.data_materials[f"{self.tcga}_preliminary_data"] = preliminary_data

        log.debug(f"{self.tcga}_preliminary_data is constructed.")

    def run_roc_curve_analysis(
            self,
            preliminary_data_name: str,
            state_variables=List[str]
    ):
        preliminary_data = self.data_materials[f"{preliminary_data_name}"].copy()

        ref_gene_column = state_variables[0]
        ref_gene_column_cohort = state_variables[1]

        roc_curve_analysis(
            reference_data_name=self.reference_data_name,
            preliminary_data=preliminary_data,
            ref_gene_column=ref_gene_column,
            cohort_specific=None
        )

        roc_curve_analysis(
            reference_data_name=self.reference_data_name,
            preliminary_data=preliminary_data,
            ref_gene_column=ref_gene_column_cohort,
            cohort_specific=self.tcga
        )

        # Baseline Non-Zero
        preliminary_data_baseline_nonzero = preliminary_data[preliminary_data['BASELINE'] != 0].copy()

        roc_curve_analysis(
            reference_data_name=self.reference_data_name,
            preliminary_data=preliminary_data_baseline_nonzero,
            ref_gene_column=ref_gene_column,
            cohort_specific=None
        )

        roc_curve_analysis(
            reference_data_name=self.reference_data_name,
            preliminary_data=preliminary_data_baseline_nonzero,
            ref_gene_column=ref_gene_column_cohort,
            cohort_specific=self.tcga
        )

