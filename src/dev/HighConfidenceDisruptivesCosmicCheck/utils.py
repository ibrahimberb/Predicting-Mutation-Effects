from datetime import datetime
from pathlib import Path
import os.path as op

from pandas import DataFrame

from IPython.display import display

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
