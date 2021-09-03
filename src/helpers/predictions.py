from typing import List

import pandas as pd
from pandas import DataFrame
from .mylogger import get_handler
import logging
import matplotlib.pyplot as plt

from tqdm.notebook import tqdm

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class Predictions(dict):
    def __init__(self, n_experiment):
        super().__init__()
        self.n_experiment = n_experiment
        self.value_counts = {}
        self.predictions_distributions_per_exp = {}

    def init_value_counts(self, tcga):
        log.debug("Initializing value counts ..")
        value_counts = []
        for exp in range(self.n_experiment):
            pred = self[tcga][exp]
            value_counts.append((len(pred[pred == 0]), len(pred[pred == 1])))
        self.value_counts[tcga] = value_counts

    def prepare_predictions_distribution_data(self, tcga):
        predictions_distributions = pd.DataFrame(
            self.value_counts[tcga],
            index=[f"EXP_{exp}" for exp in range(1, self.n_experiment + 1)],
        )
        predictions_distributions.columns = ["Disrupting", "NoEffect+Increasing"]
        predictions_distributions["Total_entry"] = (
                predictions_distributions["Disrupting"]
                + predictions_distributions["NoEffect+Increasing"]
        )

        self.predictions_distributions_per_exp[tcga] = predictions_distributions

    def plot_distributions(self, tcga):
        predictions_distributions_data = self.predictions_distributions_per_exp[tcga].copy()
        predictions_distributions_data = (
            predictions_distributions_data
            .rename_axis("EXPERIMENT")
            .reset_index()
        )
        predictions_distributions_data.plot(
            x="EXPERIMENT",
            y=["Disrupting", "NoEffect+Increasing"],
            kind="bar",
            figsize=(15, 5),
        )
        plt.grid(zorder=0, axis="y")
        plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
        plt.title("Distribution of Class Predictions per Experiment")
        plt.xlabel("Experiments")
        plt.ylabel("Counts")
        plt.show()

    def plot_predictions_distributions(self, tcga):
        self.init_value_counts(tcga)
        self.prepare_predictions_distribution_data(tcga)
        self.plot_distributions(tcga)

    def add_predictions(self, tcga, tcga_predictions):
        log.debug(f"Adding key `{tcga}` to self.predictions")
        self[tcga] = tcga_predictions

    def plot_distribution_valid_vs_invalid(self, tcga):
        valid_vs_invalid_counts_data = pd.DataFrame({
            "Num_Valid": [len(data) for data in self[f"{tcga}_predicted_valid_datasets"]],
            "Num_Invalid": [len(data) for data in self[f"{tcga}_predicted_invalid_datasets"]]
        })

        valid_vs_invalid_counts_data.index = [f'EXP_{i}' for i in range(1, self.n_experiment + 1)]
        valid_vs_invalid_counts_data.index.name = 'EXPERIMENT'
        valid_vs_invalid_counts_data.reset_index().plot(
            x="EXPERIMENT", y=["Num_Valid", "Num_Invalid"], kind="bar",
            figsize=(15, 5), color=['#9bb7bf', '#ed3e4f']
        )
        plt.grid(zorder=0, axis='y')
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.title(f"{tcga.upper()} Prediction Distribution of Number of Valid Entries vs Invalid (Dropped) Entries")
        plt.xlabel("Experiments")
        plt.ylabel("Number of Entries")
        plt.show()

    def plot_num_finalized_predictions(self, tcga):
        log.debug("Plotting number of finalized predictions per experiment.\n"
                  "Note that following plot shows the unique (protein, mutation, interactor) "
                  "triplets which had valid prediction.")
        finalized_prediction_counts_data = pd.DataFrame({
            "Num_Entries": [len(data) for data in self[f"{tcga}_finalized_prediction_dataframes"]],
        })

        finalized_prediction_counts_data.index = [f'EXP_{i}' for i in range(1, self.n_experiment + 1)]
        finalized_prediction_counts_data.index.name = 'EXPERIMENT'
        finalized_prediction_counts_data.reset_index().plot(
            x="EXPERIMENT", y=["Num_Entries"], kind="bar",
            figsize=(15, 5), color=['#24BFA5']
        )
        plt.grid(zorder=0, axis='y')
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.title(f"{tcga.upper()} Number of Unique Entries per Experiment in Finalized Prediction Dataframes")
        plt.xlabel("Experiments")
        plt.ylabel("Number of Entries")
        plt.show()

    def merge_predictions_cancer_datasets(self, tcga, tcga_data):
        log.debug(f"Merging predictions with {tcga} cancer dataset ..")
        tcga_predicted_datasets = []
        for exp in range(self.n_experiment):
            tcga_predicted_data = tcga_data.copy(deep=True)
            # Insert prediction array values into data as first (0th) column.
            tcga_predicted_data.insert(0, 'Predictions', self[f"{tcga}"][exp])
            tcga_predicted_datasets.append(tcga_predicted_data)

        self[f"{tcga}_predicted_datasets"] = tcga_predicted_datasets

    def post_process_predictions(self, tcga):
        log.debug(f"Post processing predictions for cohort {tcga} ..")
        tcga_predicted_valid_datasets = []
        tcga_predicted_invalid_datasets = []
        for exp in tqdm(range(self.n_experiment)):
            isomer_converted_data = self.convert_primary_isomer(
                column_name="Interactor_UniProt_ID",
                data=self[f"{tcga}_predicted_datasets"][exp]
            )
            predicted_valid_data, predicted_invalid_data = self.drop_invalid_predicted_entries(
                isomer_converted_data
            )

            tcga_predicted_valid_datasets.append(predicted_valid_data)
            tcga_predicted_invalid_datasets.append(predicted_invalid_data)

        self[f"{tcga}_predicted_valid_datasets"] = tcga_predicted_valid_datasets
        self[f"{tcga}_predicted_invalid_datasets"] = tcga_predicted_invalid_datasets
        log.debug(f"Preparing finalized prediction datasets for {tcga} ..")
        self.prepare_finalized_prediction_datasets(tcga, tcga_predicted_valid_datasets)
        log.debug(f"Post processing completed for {tcga}.")

    @staticmethod
    def get_predictive_columns_removed_data(data: DataFrame) -> DataFrame:
        """
        Remove predictive columns (i.e. feature columns) from given dataset, leaving them
        with predictions and triplets.

        Parameters
        ----------
            data : <DataFrame>
                The input dataframe whose predictive columns to be removed.

        Returns
        -------
            features_removed_data : <DataFrame>
                The dataframe containing prediction values along with triplet information.
                I.e. ["Predictions", "UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
        """

        features_removed_data = data[["Predictions", "UniProt_ID", "Mutation", "Interactor_UniProt_ID"]].copy()
        return features_removed_data

    @staticmethod
    def get_triplet_columns(data: DataFrame) -> DataFrame:
        """
        Remove all columns except (protein, mutation, interactor) triplets.

        Parameters
        ----------
            data : <DataFrame>
                The input dataframe whose non-triplet columns to be removed.

        Returns
        -------
            triplet_data : <DataFrame>
                The dataframe containing triplet information with following columns:
                ["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
        """

        triplet_data = data[["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]].copy()
        return triplet_data

    def get_predictive_columns_removed_datasets(
            self,
            prediction_datasets: List[DataFrame]
    ) -> List[DataFrame]:
        """
        Remove predictive columns (i.e. feature columns) from each datasets, leaving them
        with predictions and triplets.

        Parameters
        ----------
            prediction_datasets : <List[DataFrame]>
                A list of dataframe containing cancer data with predictions.

        Returns
        -------
            features_removed_datasets : <List[DataFrame]>
                List of dataframes, each containing following columns:
                ["Predictions", "UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
        """

        features_removed_datasets = []
        for exp in range(self.n_experiment):
            features_removed_data = self.get_predictive_columns_removed_data(prediction_datasets[exp])
            features_removed_datasets.append(features_removed_data)

        return features_removed_datasets

    def prepare_finalized_prediction_datasets(self, tcga, prediction_datasets):
        tcga_finalized_prediction_dataframes = self.get_predictive_columns_removed_datasets(
            prediction_datasets
        )
        self.drop_duplicated_entries_datasets(tcga_finalized_prediction_dataframes)
        self[f"{tcga}_finalized_prediction_dataframes"] = tcga_finalized_prediction_dataframes

    @staticmethod
    def drop_invalid_predicted_entries(data: DataFrame):
        """
        Prediction data contains entries which for the same (PROTEIN, MUTATION, INTERACTOR), the predicted
        class is both 0 and 1. Find such instances, and drop them.

        Parameters
        ----------
            data : <DataFrame>
                The dataframe whose invalid predicted entries will be dropped.

        Returns
        -------
            data : <DataFrame>
                Processed version of input dataframe.

            removed_entries_data : <DataFrame>
                A dataframe which contains removed entires.
        """

        entries = []
        entries_ix = []

        # For each (PROTEIN, MUTATION, INTERACTOR), capture the predicted class numbers.
        # If they are not all the same, then that (PROTEIN, MUTATION, INTERACTOR) row will be dropped.
        for index, row in data.iterrows():
            # Predicted class number(s) for current (PROTEIN, MUTATION, INTERACTOR) triplet.
            # Ideally, should be all the same. If not, then it will contain two class names.
            seach_data_predictions = data[(data["UniProt_ID"] == row["UniProt_ID"]) &
                                          (data["Mutation"] == row["Mutation"]) &
                                          (data["Interactor_UniProt_ID"] == row["Interactor_UniProt_ID"])][
                "Predictions"].unique()

            # If seach_data_predictions contains class-0 and class-1 together, then it is an invalid predicted entries.
            if len(seach_data_predictions) > 1:
                entries.append((row["Predictions"], row["UniProt_ID"], row["Mutation"], row["Interactor_UniProt_ID"]))
                entries_ix.append(index)

        removed_entries_data = pd.DataFrame(entries,
                                            columns=["PREDICTIONS", "PROTEIN", "MUTATION", "INTERACTOR"])
        log.debug('Removed entries: \n{}'.format(removed_entries_data))

        # Drop invalid predicted entries based on their index.
        data_dropped = data.drop(entries_ix, axis='index')

        # Reset index of the dataframe to avoid any possible errors.
        data_dropped.reset_index(drop=True, inplace=True)

        return data_dropped, removed_entries_data

    def drop_duplicated_entries_datasets(self, datasets: List[DataFrame]) -> None:
        """
        Dropes the duplicated entires in each dataframes.

        Parameters
        ----------
            datasets : <List[DataFrame]>

        Returns
        -------
            None, modifies the dataframes inplace.
        """
        for exp in range(self.n_experiment):
            datasets[exp].drop_duplicates(keep="first", inplace=True)

    @staticmethod
    def convert_primary_isomer(column_name: str, data: DataFrame) -> DataFrame:
        """
        Converts proteins into primary form representation (dash-free from) in given column name of the given dataframe.
        E.g.
            P16473-2 → P16473

        Parameters
        ----------
            column_name : <string>
                Name of the column where protein is stored.

            data : <DataFrame>
                The dataframe whose proteins will be processed in `column_name` column.

        Returns
        -------
            data : <DataFrame>
                Processed version of input dataframe.

        """

        # Protein names will be converted dashed-free version, if they contain.
        data[column_name] = data[column_name].apply(lambda x: x.split('-')[0])

        return data

    @staticmethod
    def get_prediction_entry(protein, mutation, interactor, data):
        predictions_array = data[
            (data['UniProt_ID'] == protein) &
            (data['Mutation'] == mutation) &
            (data['Interactor_UniProt_ID'] == interactor)
        ]['Predictions'].values

        # The entry such that it is predicted both `1` and `0` for that triplet, hence dropped.
        if predictions_array.size == 0:
            return 'NO_VOTE'

        # Prediction array contains one element, and return that element.
        elif len(predictions_array) == 1:
            [prediction] = predictions_array  # extracting single value from length-1 list.
            return prediction

        else:
            raise ValueError('There should be one entry, thus one prediction value, '
                             'but contains {} elements.'.format(len(predictions_array)))

    def add_votes(
            self,
            data: DataFrame,
            final_prediction_datasets: List[DataFrame]
    ) -> DataFrame:
        """
        Add the votes from final prediction datasets to given data.
        For each entry in data to become ensambled prediction data, the corresponding predictions will be looked up
        and placed at that entry. If corresponding prediction does not exists (indicating that prediction removed from
        that prediction data because it was an invalid prediction), "NO_VOTE" label will be assigned.

        Data will have the following form:
            +---------+----------+------------+--------------+--------------+-----+--------------+
            | Protein | Mutation | Interactor | Prediction 1 | Prediction 2 | ... | Prediction N |
            +---------+----------+------------+--------------+--------------+-----+--------------+
            |         |          |            | 1            | 0            | ... | 1            |
            +---------+----------+------------+--------------+--------------+-----+--------------+
            |         |          |            | 0            | 0            | ... | NO_VOTE      |
            +---------+----------+------------+--------------+--------------+-----+--------------+
            |         |          |            | 1            | 1            | ... | 1            |
            +---------+----------+------------+--------------+--------------+-----+--------------+

        Parameters
        ----------
            data : <DataFrame>
                The data to be ensambled prediction data, which final predictions will be added to.

            final_prediction_datasets : <List[DataFrame]>
                A list of datasets containing predictions.

        Returns
        -------
            data : <DataFrame>
                Corresponding predictions added version of input data.
        """

        final_votes = []
        for index, row in tqdm(data.iterrows(), total=len(data)):
            protein, mutation, interactor = row['UniProt_ID'], row['Mutation'], row['Interactor_UniProt_ID']

            votes = []
            for final_prediction_data in final_prediction_datasets:
                prediction = self.get_prediction_entry(protein, mutation, interactor, final_prediction_data)
                votes.append(prediction)

            final_votes.append((votes.count(0), votes.count(1), votes.count('NO_VOTE')))

        data['Num_preds_0'] = [item for item, _, _ in final_votes]
        data['Num_preds_1'] = [item for _, item, _ in final_votes]
        data['Num_preds_NO_VOTE'] = [item for _, _, item in final_votes]

        return data

    def prepare_ensambled_prediction_data(self, tcga, tcga_data):
        log.debug(f"Preparing ensambled prediction data for {tcga} ..")
        tcga_ensambled_prediction_data = self.get_triplet_columns(tcga_data)
        tcga_ensambled_prediction_data = self.convert_primary_isomer(
            "Interactor_UniProt_ID", tcga_ensambled_prediction_data
        )
        tcga_ensambled_prediction_data.drop_duplicates(keep="first", inplace=True)
        tcga_ensambled_prediction_data.reset_index(drop=True, inplace=True)
        # brca_ensamble_prediction_data.apply(lambda voting: voting()) ## todo vectorication using 3 column as param.
        tcga_ensambled_prediction_data = self.add_votes(
            tcga_ensambled_prediction_data, self[f"{tcga}_finalized_prediction_dataframes"]
        )
        self[f"{tcga}_ensambled_prediction_data"] = tcga_ensambled_prediction_data
        log.debug(f"Ensambled prediction data for {tcga} is prepared.")
