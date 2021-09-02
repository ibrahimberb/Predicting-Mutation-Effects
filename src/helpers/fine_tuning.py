from typing import List

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import RepeatedStratifiedKFold

import pandas as pd

from .machine_learning_utils import get_default_classifier

from tqdm.notebook import tqdm

from .mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)

# PARAM_GRID_SIMPLE = {
#     'bootstrap': [True, False],
#     'max_depth': [2, 5, 10],
#     'max_features': ['auto', 'sqrt'],
#     'n_estimators': [50, 100, 200, 400],
#     'min_samples_split': [2, 5, 10]
# }

PARAM_GRID_SIMPLE = {
    'max_depth': [2, 5, 10],
    'n_estimators': [10, 25, 50, 75, 100, 200, 400],
    'min_samples_split': [2, 5],
    'max_features': ['auto', 'sqrt', None],
    'class_weight': ['balanced', None]
}


class FineTuner:
    def __init__(
            self,
            n_iter: int,
            n_repeats_cv: int,
            n_jobs: int,
            verbose: int,
            Xs_determined,
            n_experiment: int,
            random_seeds: List[int],
            data_materials,
    ):
        self.n_iter = n_iter
        self.n_repeats_cv = n_repeats_cv
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.n_experiment = n_experiment
        self.random_seeds = random_seeds
        self.data_materials = data_materials
        self.Xs_determined = Xs_determined

        self.randomized_search_objects = None
        self.classifiers_attributes_data = None
        self.best_estimators = None

    def run_randomized_search(self) -> None:
        log.debug("Running randomized search for each experiment ..")
        randomized_search_objects = []
        for exp in tqdm(range(self.n_experiment)):
            random_seed = self.random_seeds[exp]
            randomized_search = self.get_randomized_search(random_seed)
            X_fine_tuning = self.data_materials[self.Xs_determined][exp]
            y_fine_tuning = self.data_materials["ys_train"][exp]
            randomized_search.fit(X_fine_tuning, y_fine_tuning)
            randomized_search_objects.append(randomized_search)
        self.randomized_search_objects = randomized_search_objects
        self.save_randomized_search_info()
        self.save_best_estimators()

    def run_search(self, search_type="randomized") -> None:
        log.debug(f"Running {search_type} search for each experiment ..")
        log.debug(f"PARAM_GRID: {PARAM_GRID_SIMPLE}")
        search_objects = []
        if search_type == "randomized":
            searcher = self.get_randomized_search
        else:
            searcher = self.get_grid_search
        # searcher = self.get_randomized_search if search_type == "randomized" else self.get_grid_search
        for exp in tqdm(range(self.n_experiment)):
            random_seed = self.random_seeds[exp]
            X_fine_tuning = self.data_materials[self.Xs_determined][exp]
            y_fine_tuning = self.data_materials["ys_train"][exp]
            randomized_search = searcher(random_seed)
            randomized_search.fit(X_fine_tuning, y_fine_tuning)
            search_objects.append(randomized_search)
        self.randomized_search_objects = search_objects
        self.save_randomized_search_info()
        self.save_best_estimators()

    def get_randomized_search(self, random_seed) -> RandomizedSearchCV:
        param_grid = PARAM_GRID_SIMPLE

        clf = get_default_classifier(random_state=random_seed)

        randomized_search = RandomizedSearchCV(clf, param_grid, n_iter=self.n_iter,
                                               ## TODO, change n_iter to 10, (or some other num.)
                                               random_state=random_seed,
                                               cv=RepeatedStratifiedKFold(n_splits=10, n_repeats=self.n_repeats_cv,
                                                                          random_state=random_seed),
                                               scoring='balanced_accuracy',
                                               return_train_score=True, n_jobs=self.n_jobs, verbose=self.verbose)

        return randomized_search

    def get_grid_search(self, random_seed) -> GridSearchCV:
        param_grid = PARAM_GRID_SIMPLE

        clf = get_default_classifier(random_state=random_seed)

        grid_search = GridSearchCV(clf, param_grid,
                                   cv=RepeatedStratifiedKFold(n_splits=10, n_repeats=self.n_repeats_cv,
                                                              random_state=random_seed),
                                   scoring='balanced_accuracy',
                                   return_train_score=True, n_jobs=self.n_jobs, verbose=self.verbose)

        return grid_search

    def save_randomized_search_info(self):
        experiment_repeat_to_randomized_search_info = {}

        for exp in range(self.n_experiment):
            experiment_repeat_to_randomized_search_info[F'EXP_{exp + 1}'] = [
                self.randomized_search_objects[exp].best_params_,
                self.randomized_search_objects[exp].best_estimator_,
                self.randomized_search_objects[exp].best_score_]

        classifiers_attributes_data = pd.DataFrame(experiment_repeat_to_randomized_search_info,
                                                   index=['best_params_', 'best_estimator_', 'best_score_']).T

        self.classifiers_attributes_data = classifiers_attributes_data

    def save_best_estimators(self):
        best_estimators = []
        for search_obj in self.randomized_search_objects:
            best_estimators.append(search_obj.best_estimator_)
        self.best_estimators = best_estimators
