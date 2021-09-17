# from helpers.common import print_annotation
from typing import List, Tuple
from pathlib import Path

from tqdm.notebook import tqdm

from helpers.helpers_predator.data_sampling import prepare_data_spsm

from helpers.helpers_predator.paths import Paths

# from helpers.data_materials import DataMaterialsML
from helpers.helpers_predator.data_materials import DataMaterials

from helpers.helpers_predator.feature_selection import ShapFeatureSelector

# from helpers.data_materials import initialize_data_materials_ML
# from helpers.data_materials import append_data_materials

# from helpers.feature_selection import AggregatedFeatureSelector

# from .helpers.feature_selection import

from helpers.helpers_predator.evaluation import EvaluationMetrics, EvaluationValid

from helpers.helpers_predator.fine_tuning import FineTuner

from helpers.helpers_predator.predictions import predictions_object

from helpers.helpers_predator.models import (
    DefaultModels,
    TunedModels,
    QualifiedModels,
    FinalizedModels,
    EnsembleVotingClassifier)


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
    ):

        log.debug('Initializing Predator ..')
        self.n_experiment = n_experiment
        self.random_seeds = list(
            range(1, self.n_experiment + 1)
        )  # todo: this will be a random array, too.

        self.paths = Paths(
            project_common_file_dir,
            mutations_path,
            tcga_code_path_pairs,
            initial_columns_path,
        )

        self.tcga_cohorts = []
        for tcga, _ in tcga_code_path_pairs:
            self.tcga_cohorts.append(tcga)

        self.data_materials = DataMaterials(n_experiment, self.random_seeds)

        self.data_materials.initialize_train_datasets(
            project_common_file_dir, initial_columns_path, mutations_path
        )
        self.data_materials.initialize_target_datasets(
            project_common_file_dir, initial_columns_path, tcga_code_path_pairs
        )

        self.default_models = DefaultModels(self.n_experiment)
        self.tuned_models = None
        self.finalized_models = None
        self.qualified_models = None

        self.eval_valid = EvaluationValid(
            self.n_experiment,
            self.data_materials,
            self.default_models,
            self.tuned_models,
            self.qualified_models,
        )

        self.determined_feature_set = None
        self.determined_features = None

        self.shap_feature_selector = None
        self.eval_metrics = None
        self.fine_tuner = None

        self.ensemble_voting_classifier = None
        self.predictions = None

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

    def run_evaluate_valid(
        self,
        models_type="default",
        show_confusion_matrix=False,
    ):

        if models_type == 'tuned':
            self.eval_valid.tuned_models = self.tuned_models

        self.eval_valid.evaluate(
            models_type=models_type,
            show_confusion_matrix=show_confusion_matrix,
            determined_features=self.determined_features
        )

    def init_shap_feature_selector(self, shap_top_ns):
        self.shap_feature_selector = ShapFeatureSelector(self.n_experiment, shap_top_ns)
        self.shap_feature_selector.load_shap_values(self.data_materials)  # fixme -> select_features
        self.shap_feature_selector.get_selected_features(self.data_materials)

    def aggregate_selected_features(self, method):
        self.shap_feature_selector.shap_aggregate_selected_features(method)
        for shap_top_n, aggregated_feature in self.shap_feature_selector.n_features_to_aggregated_features.items():
            self.data_materials.initialize_feature_selected_data_materials(shap_top_n)
            self.data_materials.append_feature_selected_data_materials(shap_top_n, aggregated_feature)

    def initialize_evaluation_metrics(self):
        self.eval_metrics = EvaluationMetrics(self.n_experiment, self.data_materials, self.shap_feature_selector)

    def set_determined_feature_set(self, feature_set: str):
        _, top_n = feature_set.split('_')
        determined_features = (
            self
            .shap_feature_selector
            .n_features_to_aggregated_features[int(top_n)]
        )
        log.debug(f"Setting determined feature set to `{feature_set}`.")
        self.determined_feature_set = feature_set
        log.debug(f"Setting determined features to \n{determined_features}.")
        self.determined_features = determined_features

    def init_fine_tuner(self, n_iter, n_repeats_cv, n_jobs, verbose):
        # Fine tuning will be applied to training set, not all training data.
        Xs_determined = f"Xs_train_{self.determined_feature_set}"

        self.fine_tuner = FineTuner(
            n_iter=n_iter, n_repeats_cv=n_repeats_cv, n_jobs=n_jobs, verbose=verbose,
            Xs_determined=Xs_determined, n_experiment=self.n_experiment,
            random_seeds=self.random_seeds, data_materials=self.data_materials
        )

    def run_search(self, search_type):
        self.fine_tuner.run_search(search_type)
        self.tuned_models = TunedModels(self.fine_tuner.best_estimators)

    def fit_finalized_models(self):
        log.debug("Fitting finalized models with all training data ..")
        self.finalized_models = FinalizedModels(self.tuned_models)
        self.finalized_models.fit_all(self.data_materials, self.determined_feature_set)

    def initialize_target_data_materials(self):
        self.data_materials.initialize_target_data_materials(
            self.determined_features, self.paths.tcga_code_path_pairs
        )

    def predict(self, voting="hard"):
        log.debug("Predicting on cancer datasets ..")
        self.ensemble_voting_classifier = EnsembleVotingClassifier(
            self.finalized_models, voting=voting
        )

        self.predictions = predictions_object(voting, self.n_experiment)

        for tcga in self.tcga_cohorts:
            log.debug(f"Predicting on {tcga} cohort ..")
            tcga_predictions = self.ensemble_voting_classifier.predict(
                self.data_materials[f"Xs_{tcga}"]
            )
            self.predictions.add_predictions(tcga, tcga_predictions)

    def predictions_post_process(self):
        for tcga in self.tcga_cohorts:
            self.predictions.post_process_predictions(
                tcga,
                self.data_materials[f"{tcga}"]
            )

    def prepare_ensemble_prediction_data(self):
        for tcga in self.tcga_cohorts:
            self.predictions.prepare_ensemble_prediction_data(
                tcga, self.data_materials[f"{tcga}"]
            )

#
# predator = Predator(project_common_file_dir=PROJECT_COMMON_FILE_DIR,
#                     mutations_path=MUTATIONS_PATH,
#                     tcga_code_path_pairs=[('brca', BRCA_PATH)],
#                     initial_columns_path=INITIAL_COLUMNS_PATH,
#                     n_experiment=50)
