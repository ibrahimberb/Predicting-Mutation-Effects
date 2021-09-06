# from helpers.common import print_annotation
from typing import List, Tuple
from pathlib import Path

from helpers.machine_learning_utils import get_default_classifier, evaluate_valid

from helpers.log_script import ColorHandler
import logging
from helpers.loaders import load_train_data, load_tcga_data
from helpers.preprocessing import preprocess_train_data, preprocess_tcga_data
from tqdm.notebook import tqdm

from helpers.prepare_data_for_ML import prepare_data_machine_learning

from helpers.data_sampling import prepare_data_spsm

from helpers.data_materials import DataMaterialsML

# from helpers.feature_selection import AggregatedFeatureSelector

log = logging.Logger("debug_runner", level=logging.INFO)
log.addHandler(ColorHandler())

# PATHS
PROJECT_COMMON_FILE_DIR = Path("../../data/")
MUTATIONS_PATH = Path("training_data_M1.txt")
INITIAL_COLUMNS_PATH = Path("../../data/initial_columns_59.csv")
BRCA_PATH = Path("../../data/BRCA_INTERFACE_A2.txt")

TCGA_CODE = str


class Predator:
    def __init__(
        self,
        project_common_file_dir: Path,
        mutations_path: Path,
        tcga_code_path_pairs: List[Tuple[TCGA_CODE, Path]],
        initial_columns_path: Path,
        n_experiment: int,
        feature_selection_args=None,
    ):

        self.n_experiment = n_experiment
        self.random_seeds = list(
            range(1, self.n_experiment + 1)
        )  # todo: this will be a random array, too.

        self.train_data = load_train_data(project_common_file_dir, mutations_path)
        self.train_data_processed = preprocess_train_data(
            self.train_data, initial_columns_path
        )

        for tcga_code_path_pair in tcga_code_path_pairs:
            tcga_name, tcga_path = tcga_code_path_pair
            setattr(
                self,
                tcga_name.lower(),
                load_tcga_data(project_common_file_dir, tcga_path),
            )
            setattr(
                self,
                "target_" + tcga_name.lower() + "_data",
                preprocess_tcga_data(
                    getattr(self, tcga_name.lower()), initial_columns_path
                ),
            )

        self.sampled_train_data_list = None
        self.data_materials_ML = DataMaterialsML()

        self.selected_features = None
        # self.aggregated_feature_selector = AggregatedFeatureSelector(features_list=[[]], aggregation_method=None)

        # todo
        self.Xs_train_benchmark_feature_names_dataframes_list = None

        # fixme: code refactor.
        # feature_selection_args = {n_features, }

        # -------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------

        # variant = spsm
        # self.variant = variant

    def sample_spsm(self):
        log.debug("sampling ..")
        sampled_train_data_list = []
        for i in tqdm(range(self.n_experiment)):
            sampled_train_data = prepare_data_spsm(
                self.train_data_processed, random_seed=self.random_seeds[i]
            )
            sampled_train_data_list.append(sampled_train_data)

        self.sampled_train_data_list = sampled_train_data_list

    def prepare_data_matarials_ML(self):
        self.data_materials_ML = DataMaterialsML()
        log.debug("preparing data materials for ML ..")
        for i in tqdm(range(self.n_experiment)):
            data_materials = prepare_data_machine_learning(
                self.sampled_train_data_list[i], random_seed=self.random_seeds[i]
            )
            self.data_materials_ML.append_data_materials(data_materials)

    def run_evaluate_valid(self, model_type="default", confusion_matrix=False) -> Tuple[List[float], List[float]]:
        if model_type == "default":
            classifiers = [get_default_classifier()] * self.n_experiment
        else:
            # TODO
            classifiers = [get_default_classifier()] * self.n_experiment
            # classifiers = self.classifiers

        log.debug("Evaluating on validation data ..")
        accuracy_scores = []
        balanced_accuracy_scores = []
        for i in tqdm(range(self.n_experiment)):
            print("-------- EXPERIMENT: {:>2} --------".format(i + 1))
            acc_score, balan_acc_score = evaluate_valid(
                classifiers[i],
                self.data_materials_ML.Xs_train[i],
                self.data_materials_ML.ys_train[i],
                self.data_materials_ML.Xs_valid[i],
                self.data_materials_ML.ys_valid[i],
                show_confusion_matrix=confusion_matrix,
            )
            accuracy_scores.append(acc_score)
            balanced_accuracy_scores.append(balan_acc_score)
            print("================================")

        return accuracy_scores, balanced_accuracy_scores

    def prepare_feature_selected_data_materials(self, shap_top_n: int):
        log.debug("preparing data materials for ML ..")
        for i in tqdm(range(self.n_experiment)):
            self.data_materials_ML.append_feature_selected_data_materials(shap_top_n)
            data_materials = prepare_data_machine_learning(
                self.sampled_train_data_list[i], random_seed=self.random_seeds[i]
            )


    def evaluation_metrics(self):
        # second thought, it's not needed.
        raise NotImplementedError


#
# predator = Predator(project_common_file_dir=PROJECT_COMMON_FILE_DIR,
#                     mutations_path=MUTATIONS_PATH,
#                     tcga_code_path_pairs=[('brca', BRCA_PATH)],
#                     initial_columns_path=INITIAL_COLUMNS_PATH,
#                     n_experiment=50)
