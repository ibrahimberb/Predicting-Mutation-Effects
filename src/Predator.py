# from helpers.common import print_annotation
from typing import List, Tuple
from pathlib import Path

from helpers.machine_learning_utils import get_default_classifier, evaluate_valid

from tqdm.notebook import tqdm

from helpers.data_sampling import prepare_data_spsm

# from helpers.data_materials import DataMaterialsML
from helpers.data_materials import DataMaterials

from helpers.feature_selection import ShapFeatureSelector

# from helpers.data_materials import initialize_data_materials_ML
# from helpers.data_materials import append_data_materials

# from helpers.feature_selection import AggregatedFeatureSelector

# from .helpers.feature_selection import

from helpers.mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


# PATHS
PROJECT_COMMON_FILE_DIR = Path("../data/")
MUTATIONS_PATH = Path("training_data_M1.txt")
INITIAL_COLUMNS_PATH = Path("../data/initial_columns_59.csv")
BRCA_PATH = Path("../data/BRCA_INTERFACE_A2.txt")

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

        log.debug('Initializing Predator ..')
        self.n_experiment = n_experiment
        self.random_seeds = list(
            range(1, self.n_experiment + 1)
        )  # todo: this will be a random array, too.

        self.data_materials = DataMaterials(n_experiment, self.random_seeds)

        self.data_materials.initialize_train_datasets(project_common_file_dir, initial_columns_path, mutations_path)
        self.data_materials.initialize_target_datasets(project_common_file_dir, initial_columns_path, tcga_code_path_pairs)

        self.shap_feature_selector = None

        ##########################
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
                self.data_materials["train_data_processed"], random_seed=self.random_seeds[i]
            )
            sampled_train_data_list.append(sampled_train_data)

        self.data_materials["sampled_train_data_list"] = sampled_train_data_list

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
                self.data_materials["Xs_train"][i],
                self.data_materials["ys_train"][i],
                self.data_materials["Xs_valid"][i],
                self.data_materials["ys_valid"][i],
                show_confusion_matrix=confusion_matrix,
            )
            accuracy_scores.append(acc_score)
            balanced_accuracy_scores.append(balan_acc_score)
            print("================================")

        return accuracy_scores, balanced_accuracy_scores

    def init_shap_feature_selector(self, shap_top_ns):
        self.shap_feature_selector = ShapFeatureSelector(self.n_experiment, shap_top_ns)
        self.shap_feature_selector.load_shap_values(self.data_materials)
        self.shap_feature_selector.get_selected_features(self.data_materials)

        #
        # for shap_top_n in shap_top_ns:
        #     self.data_materials.initialize_feature_selected_data_materials(shap_top_n)
        #     self.data_materials.append_feature_selected_data_materials(shap_top_n, )

    def aggregate_selected_features(self, method):
        self.shap_feature_selector.aggregate_selected_features(method)
        for shap_top_n, aggregated_feature in self.shap_feature_selector.n_features_to_aggregated_features.items():
            self.data_materials.initialize_feature_selected_data_materials(shap_top_n)
            self.data_materials.append_feature_selected_data_materials(shap_top_n, aggregated_feature)



    # def prepare_feature_selected_data_materials(self, shap_top_n: int):
    #     log.debug("preparing data materials for ML ..")
    #     for i in tqdm(range(self.n_experiment)):
    #         self.data_materials_ML.append_feature_selected_data_materials(shap_top_n)
    #         data_materials = prepare_data_machine_learning(
    #             self.sampled_train_data_list[i], random_seed=self.random_seeds[i]
    #         )


    def evaluation_metrics(self):
        # second thought, it's not needed.
        raise NotImplementedError


#
# predator = Predator(project_common_file_dir=PROJECT_COMMON_FILE_DIR,
#                     mutations_path=MUTATIONS_PATH,
#                     tcga_code_path_pairs=[('brca', BRCA_PATH)],
#                     initial_columns_path=INITIAL_COLUMNS_PATH,
#                     n_experiment=50)
