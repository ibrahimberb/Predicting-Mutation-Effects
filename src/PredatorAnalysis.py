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

from helpers.helpers_analysis.loaders import ReferenceDataset
from helpers.helpers_analysis.get_fpr_tpr_ths import get_fpr_tpr_ths

import matplotlib.pyplot as plt
import seaborn as sns

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
            reference_data_paths: Tuple[CohortSpecificReferenceDataPath, ReferenceDataPath]
    ):
        self.tcga = tcga
        self.snv_path = snv_path
        self.prediction_data_path = prediction_data_path
        self.elaspic_core_path = elaspic_core_path
        self.elaspic_interface_path = elaspic_interface_path
        self.reference_data_name = reference_data_name
        self.reference_data_spec_cohort_path = reference_data_paths[0]
        self.reference_data_path = reference_data_paths[1]
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
        log.debug(f"Adding Reference Dataset Columns: General and Cohort Specific columns ..")
        if self.reference_data_name == ReferenceDataset.CANCERMINE:
            add_cancermine_status(
                preliminary_data,
                cancermine_genes=self.data_materials["cancermine_all_genes"],
                cancermine_cohort_genes=self.data_materials[f"cancermine_{self.tcga}_genes"],
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
            test_variables=List[str]
    ):
        preliminary_data = self.data_materials[f"{preliminary_data_name}"].copy()

        ref_gene_column = test_variables[0]
        ref_gene_column_cohort = test_variables[1]

        baseline_counts_vs_cgc_status_data = preliminary_data[["BASELINE", ref_gene_column]].copy()
        our_method_counts_vs_cgc_status_data = preliminary_data[["OUR_METHOD", ref_gene_column]].copy()
        elaspic_cov_vs_cgc_status_data = preliminary_data[["ELASPIC_COVERAGE", ref_gene_column]].copy()

        baseline_counts_vs_cgc_status_tcga_data = preliminary_data[["BASELINE", ref_gene_column_cohort]].copy()
        our_method_counts_vs_cgc_status_tcga_data = preliminary_data[["OUR_METHOD", ref_gene_column_cohort]].copy()
        elaspic_cov_vs_cgc_status_tcga_data = preliminary_data[["ELASPIC_COVERAGE", ref_gene_column_cohort]].copy()

        fpr_baseline, tpr_baseline, ths_baseline = get_fpr_tpr_ths(
            np.array(baseline_counts_vs_cgc_status_data["CancerMine_STATUS"]),  # FIXME: ref_gene_column
            np.array(baseline_counts_vs_cgc_status_data["BASELINE"]))

        fpr_our_method, tpr_our_method, ths_our_method = get_fpr_tpr_ths(
            np.array(our_method_counts_vs_cgc_status_data["CancerMine_STATUS"]),
            np.array(our_method_counts_vs_cgc_status_data["OUR_METHOD"]))

        fpr_elaspic_cov, tpr_elaspic_cov, ths_elaspic_cov = get_fpr_tpr_ths(
            np.array(elaspic_cov_vs_cgc_status_data["CancerMine_STATUS"]),
            np.array(elaspic_cov_vs_cgc_status_data["ELASPIC_COVERAGE"]))

        fpr_baseline_brca, tpr_baseline_brca, ths_baseline_brca = get_fpr_tpr_ths(
            np.array(baseline_counts_vs_cgc_status_tcga_data["CancerMine_STATUS (BRCA)"]),
            np.array(baseline_counts_vs_cgc_status_tcga_data["BASELINE"]))

        fpr_our_method_brca, tpr_our_method_brca, ths_our_method_brca = get_fpr_tpr_ths(
            np.array(our_method_counts_vs_cgc_status_tcga_data["CancerMine_STATUS (BRCA)"]),
            np.array(our_method_counts_vs_cgc_status_tcga_data["OUR_METHOD"]))

        fpr_elaspic_cov_brca, tpr_elaspic_cov_brca, ths_elaspic_cov_brca = get_fpr_tpr_ths(
            np.array(elaspic_cov_vs_cgc_status_tcga_data["CancerMine_STATUS (BRCA)"]),
            np.array(elaspic_cov_vs_cgc_status_tcga_data["ELASPIC_COVERAGE"]))

        roc_auc_baseline = metrics.auc(fpr_baseline, tpr_baseline)
        roc_auc_our_method = metrics.auc(fpr_our_method, tpr_our_method)
        roc_auc_elaspic_cov = metrics.auc(fpr_elaspic_cov, tpr_elaspic_cov)

        roc_auc_baseline_tcga = metrics.auc(fpr_baseline_brca, tpr_baseline_brca)
        roc_auc_our_method_tcga = metrics.auc(fpr_our_method_brca, tpr_our_method_brca)
        roc_auc_elaspic_cov_tcga = metrics.auc(fpr_elaspic_cov_brca, tpr_elaspic_cov_brca)

        # LAST HERE ....
        # fixme: refactor
        # Plotting CancerMine_STATUS
        plt.figure(figsize=(7, 5))
        sns.set(style='white', font_scale=1.50)

        plt.title('Receiver Operating Characteristic (ROC)\n$CancerMine\ STATUS$ vs Various Columns')
        plt.plot(fpr_baseline, tpr_baseline, label='baseline (%0.3f)' % roc_auc_baseline)
        plt.plot(fpr_our_method, tpr_our_method, label='our_method (%0.3f)' % roc_auc_our_method)
        plt.plot(fpr_elaspic_cov, tpr_elaspic_cov, label='elaspic_cov (%0.3f)' % roc_auc_elaspic_cov)

        plt.legend(loc='lower right')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.ylabel('True Positive Rate\n(Sensitivity)')
        plt.xlabel('False Positive Rate\n(1-Specificity)')
        plt.show()

        # Plotting CancerMine_STATUS_brca
        plt.figure(figsize=(7, 5))
        sns.set(style='white', font_scale=1.50)

        plt.title('Receiver Operating Characteristic (ROC)\n$CancerMine\ BRCA\ STATUS$ vs Various Columns')
        plt.plot(fpr_baseline_brca, tpr_baseline_brca, label='baseline BRCA  (%0.3f)' % roc_auc_baseline_tcga)
        plt.plot(fpr_our_method_brca, tpr_our_method_brca, label='our_method BRCA  (%0.3f)' % roc_auc_our_method_tcga)
        plt.plot(fpr_elaspic_cov_brca, tpr_elaspic_cov_brca,
                 label='elaspic_cov BRCA  (%0.3f)' % roc_auc_elaspic_cov_tcga)

        plt.legend(loc='lower right')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.ylabel('True Positive Rate\n(Sensitivity)')
        plt.xlabel('False Positive Rate\n(1-Specificity)')
        plt.show()

        # # ROC with `BASELINE` Non-zero Columns # # # # # # # # # # # # # # #
        preliminary_data_baseline_nonzero = preliminary_data[preliminary_data['BASELINE'] != 0].copy()
        print(preliminary_data_baseline_nonzero.shape)

        baseline_counts_vs_cgc_status_bnz_data = preliminary_data_baseline_nonzero[["BASELINE", ref_gene_column]].copy()
        our_method_counts_vs_cgc_status_bnz_data = preliminary_data_baseline_nonzero[["OUR_METHOD", ref_gene_column]].copy()
        elaspic_cov_vs_cgc_status_bnz_data = preliminary_data_baseline_nonzero[["ELASPIC_COVERAGE", ref_gene_column]].copy()

        baseline_counts_vs_cgc_status_brca_bnz_data = preliminary_data_baseline_nonzero[
            ["BASELINE", ref_gene_column_cohort]].copy()
        our_method_counts_vs_cgc_status_brca_bnz_data = preliminary_data_baseline_nonzero[
            ["OUR_METHOD", ref_gene_column_cohort]].copy()
        elaspic_cov_vs_cgc_status_brca_bnz_data = preliminary_data_baseline_nonzero[
            ["ELASPIC_COVERAGE", ref_gene_column_cohort]].copy()

        fpr_baseline_bnz, tpr_baseline_bnz, ths_baseline_bnz = get_fpr_tpr_ths(
            np.array(baseline_counts_vs_cgc_status_bnz_data["CancerMine_STATUS"]),
            np.array(baseline_counts_vs_cgc_status_bnz_data["BASELINE"]))

        fpr_our_method_bnz, tpr_our_method_bnz, ths_our_method_bnz = get_fpr_tpr_ths(
            np.array(our_method_counts_vs_cgc_status_bnz_data["CancerMine_STATUS"]),
            np.array(our_method_counts_vs_cgc_status_bnz_data["OUR_METHOD"]))

        fpr_elaspic_cov_bnz, tpr_elaspic_cov_bnz, ths_elaspic_cov_bnz = get_fpr_tpr_ths(
            np.array(elaspic_cov_vs_cgc_status_bnz_data["CancerMine_STATUS"]),
            np.array(elaspic_cov_vs_cgc_status_bnz_data["ELASPIC_COVERAGE"]))

        fpr_baseline_bnz_brca, tpr_baseline_bnz_brca, ths_baseline_bnz_brca = get_fpr_tpr_ths(
            np.array(baseline_counts_vs_cgc_status_brca_bnz_data["CancerMine_STATUS (BRCA)"]),
            np.array(baseline_counts_vs_cgc_status_brca_bnz_data["BASELINE"]))

        fpr_our_method_bnz_brca, tpr_our_method_bnz_brca, ths_our_method_bnz_brca = get_fpr_tpr_ths(
            np.array(our_method_counts_vs_cgc_status_brca_bnz_data["CancerMine_STATUS (BRCA)"]),
            np.array(our_method_counts_vs_cgc_status_brca_bnz_data["OUR_METHOD"]))

        fpr_elaspic_cov_bnz_brca, tpr_elaspic_cov_bnz_brca, ths_elaspic_cov_bnz_brca = get_fpr_tpr_ths(
            np.array(elaspic_cov_vs_cgc_status_brca_bnz_data["CancerMine_STATUS (BRCA)"]),
            np.array(elaspic_cov_vs_cgc_status_brca_bnz_data["ELASPIC_COVERAGE"]))

        roc_auc_baseline_bnz = metrics.auc(fpr_baseline_bnz, tpr_baseline_bnz)
        roc_auc_our_method_bnz = metrics.auc(fpr_our_method_bnz, tpr_our_method_bnz)
        roc_auc_elaspic_cov_bnz = metrics.auc(fpr_elaspic_cov_bnz, tpr_elaspic_cov_bnz)

        roc_auc_baseline_bnz_brca = metrics.auc(fpr_baseline_bnz_brca, tpr_baseline_bnz_brca)
        roc_auc_our_method_bnz_brca = metrics.auc(fpr_our_method_bnz_brca, tpr_our_method_bnz_brca)
        roc_auc_elaspic_cov_bnz_brca = metrics.auc(fpr_elaspic_cov_bnz_brca, tpr_elaspic_cov_bnz_brca)

        # Plotting BNZ CancerMine
        plt.figure(figsize=(7, 5))
        sns.set(style='white', font_scale=1.50)

        plt.title('Receiver Operating Characteristic (ROC)\n$CancerMine\ STATUS$ vs Various Columns')
        plt.plot(fpr_baseline_bnz, tpr_baseline_bnz, label='baseline bnz (%0.3f)' % roc_auc_baseline_bnz)
        plt.plot(fpr_our_method_bnz, tpr_our_method_bnz, label='our_method bnz (%0.3f)' % roc_auc_our_method_bnz)
        plt.plot(fpr_elaspic_cov_bnz, tpr_elaspic_cov_bnz, label='elaspic_cov bnz (%0.3f)' % roc_auc_elaspic_cov_bnz)

        plt.legend(loc='lower right')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.ylabel('True Positive Rate\n(Sensitivity)')
        plt.xlabel('False Positive Rate\n(1-Specificity)')
        plt.show()

        # Plotting BNZ CancerMine_STATUS_brca
        plt.figure(figsize=(7, 5))
        sns.set(style='white', font_scale=1.25)

        plt.title('Receiver Operating Characteristic (ROC)\n$CancerMine\ BRCA\ STATUS$ vs Various Columns')
        plt.plot(fpr_baseline_bnz_brca, tpr_baseline_bnz_brca, label='baseline (%0.3f)' % roc_auc_baseline_bnz_brca)
        plt.plot(fpr_our_method_bnz_brca, tpr_our_method_bnz_brca,
                 label='disruptive interactions (%0.3f)' % roc_auc_our_method_bnz_brca)
        plt.plot(fpr_elaspic_cov_bnz_brca, tpr_elaspic_cov_bnz_brca,
                 label='coverage (%0.3f)' % roc_auc_elaspic_cov_bnz_brca)

        plt.legend(loc='lower right')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.ylabel('$True Positive Rate\n(Sensitivity)$')
        plt.xlabel('False Positive Rate\n(1-Specificity)')
        plt.show()

        # displayers
        display(pd.DataFrame({
            "roc_auc_baseline": [roc_auc_baseline],
            "roc_auc_our_method": [roc_auc_our_method],
            "roc_auc_elaspic_cov": [roc_auc_elaspic_cov],
        }, index=['AUC']).T)

        display(pd.DataFrame({
            "roc_auc_baseline_brca": [roc_auc_baseline_tcga],
            "roc_auc_our_method_brca": [roc_auc_our_method_tcga],
            "roc_auc_elaspic_cov_brca": [roc_auc_elaspic_cov_tcga],
        }, index=['AUC']).T)

        display(pd.DataFrame({
            "roc_auc_baseline_bnz": [roc_auc_baseline_bnz],
            "roc_auc_our_method_bnz": [roc_auc_our_method_bnz],
            "roc_auc_elaspic_cov_bnz": [roc_auc_elaspic_cov_bnz],
        }, index=['AUC']).T)

        display(pd.DataFrame({
            "roc_auc_baseline_bnz_brca": [roc_auc_baseline_bnz_brca],
            "roc_auc_our_method_bnz_brca": [roc_auc_our_method_bnz_brca],
            "roc_auc_elaspic_cov_bnz_brca": [roc_auc_elaspic_cov_bnz_brca],
        }, index=['AUC']).T)



# predator_analysis = PredatorAnalysis(
#     tcga_code="brca",
#     snv_path=SNV_BRCA_PATH,
#     prediction_data_path=PREDICTION_BRCA_REDUCED_PATH,
#     elaspic_core_path=BRCA_CORE_PATH,
#     elaspic_interface_path=BRCA_INTERFACE_PATH,
#     reference_data="cancermine",
#     reference_data_path=CANCER_MINE_BREAST_PATH,
# )
