from datetime import datetime
from pathlib import Path
import os.path as op

from pandas import DataFrame

from IPython.display import display

from src.dev.HighConfidenceDisruptivesCosmicCheck.CosmicAPI.utils.misc import get_residue_position
from src.helpers.helpers_analysis.gene_id_retrieval import GeneIDFetcher

UNIPROT_GENE_MAPPING_PATH = "../../helpers/helpers_analysis/gene_retrieval/UNIPROT_GENE_MAPPING.csv"


class HighConfidenceDisruptiveMutationsHelper:
    def __init__(self, tcga: str, data: DataFrame, confidence: float):
        """
        :param tcga: The TCGA cohort name
        :param data: Prediction data
        :param confidence: Disruptive confidence should be between 0.50 and 1.
        """
        if confidence < 0.5:
            raise ValueError("If confidence is less then 0.50, it is not predicted as disruptive.")

        self.tcga = tcga.upper()
        self.data = data
        self.confidence = confidence
        self.gene_id_fetcher = GeneIDFetcher(UNIPROT_GENE_MAPPING_PATH)

        self.high_confidence_data = None
        self.prepare_high_confidence_disruptive_mutations()

    def get_high_confidence_disruptive_mutations(self):
        return self.high_confidence_data

    def prepare_high_confidence_disruptive_mutations(self):
        """
        Confidence is between 0 and 1
        """

        # Since Median probability is the probability of being class 1, that is increasing or no effect.
        self.data["Disruptive_probability"] = 1 - self.data["Median_Probability"]

        high_confidence_data = self.data[self.data["Disruptive_probability"] >= self.confidence]

        assert len(high_confidence_data["Prediction"].unique()) == 1
        assert high_confidence_data["Prediction"].unique()[0] == 0

        high_confidence_data.insert(
            loc=0,
            column="TCGA",
            value=self.tcga
        )

        high_confidence_data.insert(
            loc=1,
            column="GENE",
            value=high_confidence_data["UniProt_ID"].apply(lambda x: self.gene_id_fetcher.fetch(x))
        )

        high_confidence_data.insert(
            loc=4,
            column="INTERACTOR_GENE",
            value=high_confidence_data["Interactor_UniProt_ID"].apply(lambda x: self.gene_id_fetcher.fetch(x))
        )

        self.high_confidence_data = high_confidence_data

    def extract_high_confidence_disruptive_mutations(self, view=False):
        folder_path = op.join(
            "HighConfidenceDisruptiveData", f"confidence_{self.confidence:.2f}"
        )
        file_date = datetime.today().strftime('%Y-%m-%d')
        Path(f"{folder_path}").mkdir(parents=True, exist_ok=True)
        file_name = f"{self.tcga}_confidence_{self.confidence:.2f}_{file_date}.csv"
        file_path = op.join(folder_path, file_name)

        if op.isfile(file_path):
            raise FileExistsError(f"You already have the file {file_path}")

        high_confidence_data = self.get_high_confidence_disruptive_mutations()

        if view:
            display(high_confidence_data)

        high_confidence_data.to_csv(file_path, index=False)

        print(f"{self.tcga} data is extracted to {file_path} successfully.")


class CosmicResultsAttaching:
    def __init__(self, cosmic_results_data):
        self.cosmic_results_data = cosmic_results_data

    @staticmethod
    def find_in_cosmic_results(
            cosmic_results_data: DataFrame,
            gene: str,
            mut: str
    ) -> dict:

        query = cosmic_results_data[
            (cosmic_results_data["GENE"] == gene) &
            (cosmic_results_data["RESIDUE_POSITION"] == int(get_residue_position(mut)))
        ]

        if query.empty:
            query_result = {
                "CGC_status": "NOT_FOUND",
                "most_significant_codon_tier": "NOT_FOUND",
            }

        else:
            [CGC_status] = query["CGC_STATUS"]
            [most_significant_codon_tier] = query["MOST_SIGNIFICANT_CODON_TIER"]

            query_result = {
                "CGC_status": CGC_status,
                "most_significant_codon_tier": most_significant_codon_tier,
            }

        return query_result

    def attach_results(
            self,
            tcga_prediction_data: DataFrame,
    ) -> DataFrame:

        tcga_prediction_data_cosmic_results = tcga_prediction_data.copy()
        cosmic_results = tcga_prediction_data_cosmic_results.apply(
            lambda row: self.find_in_cosmic_results(
                cosmic_results_data=self.cosmic_results_data,
                gene=row["GENE"],
                mut=row["Mutation"]
            ), axis=1
        )

        tcga_prediction_data_cosmic_results["CGC_status"] = cosmic_results.apply(
            lambda row: row["CGC_status"]
        )

        tcga_prediction_data_cosmic_results["MOST_SIGNIFICANT_CODON_TIER"] = cosmic_results.apply(
            lambda row: row["most_significant_codon_tier"]
        )

        return tcga_prediction_data_cosmic_results
