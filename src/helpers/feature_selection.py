from typing import List
from .data_materials import DataMaterials
from tqdm.notebook import tqdm
from .machine_learning_utils import get_default_classifier

from collections import defaultdict

from pandas import DataFrame

import numpy as np

import shap


from .mylogger import get_handler
import logging

handler = get_handler()
handler_simple = get_handler('simple')

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)

log_simple = logging.getLogger('Feature_selection')
log_simple.handlers[:] = []
log_simple.addHandler(handler_simple)
log_simple.setLevel(logging.DEBUG)

# ## i was here.
# @dataclass  # fixme: namedtuple really necessary?
# class AggregationMethods(namedtuple):
#     occurance = 'occurrance'
#     concordance = 'concordance'


class ShapFeatureSelector:
    def __init__(self, n_experiment, shap_top_ns: List[int]):
        log.debug("Initializing ShapFeatureSelector ..")
        self.shap_top_ns = shap_top_ns
        self.n_experiment = n_experiment
        shap.initjs()

        self.shap_values_train_list = None

        self.n_features_to_selected_features_list = defaultdict(list)  # Dict[int, List[List[str]]]
        self.n_features_to_combined_features = defaultdict(list)  # Dict[int, List[str]]

    def load_shap_values(self, data_materials: DataMaterials):
        """Fitted on trainsets, not on all train data."""
        log.debug("Loading ShapFeatureSelector ..")
        shap_values_train_list = []
        for i in tqdm(range(self.n_experiment)):
            model = get_default_classifier()
            model.fit(data_materials["Xs_train"][i], data_materials["ys_train"][i])
            explainer = shap.TreeExplainer(model)
            shap_values_train_list.append(explainer.shap_values(data_materials["Xs_train"][i],
                                                                approximate=False, check_additivity=False))

        self.shap_values_train_list = shap_values_train_list

    @staticmethod
    def get_selected_features_single_data(shap_values, features_data: DataFrame, top_n) -> List[str]:
        column_list = features_data.columns
        feature_ratio = (np.abs(shap_values).sum(0) / np.abs(shap_values).sum()) * 100
        column_list = column_list[np.argsort(feature_ratio)[::-1]]
        column_list = column_list[:top_n]

        return list(column_list)

    def get_selected_features(self, data_materials: DataMaterials):
        log_simple.debug(" === SELECTED FEATURES === ")
        for top_n in self.shap_top_ns:
            log_simple.debug(f" --- SHAP TOP {top_n} ---")
            for exp in range(self.n_experiment):
                shap_values = self.shap_values_train_list[exp][1]
                features_data = data_materials["Xs_train"][exp]
                selected_features = self.get_selected_features_single_data(shap_values, features_data, top_n=top_n)
                log_simple.debug(f"[Experiment {exp+1}] ")
                print(f"{selected_features}\n")
                self.n_features_to_selected_features_list[top_n].append(selected_features)

    def combine_selected_features(self):
        pass


    # def init_shap_feature_selector(self, data_materials: DataMaterials):
    #     self.load_shap_values(data_materials)
    #     self.print_selected_features(data_materials)

    # def append_top_n_features(self, data_materials: DataMaterials):
    #     log.debug("Appending feature datasets with selected columns to data_materials ..")
    #     for exp in range(self.n_experiment):
    #         for top_n in self.shap_top_ns:
    #             shap_values = self.shap_values_train_list[exp][1]
    #             features_data = data_materials["Xs_train"][exp]
    #             top_n_features = self.get_selected_features_single_data(shap_values, features_data, top_n=top_n)
    #             data_materials[f"Xs_train_shap_{self.shap_top_ns}"].append(
    #                 features_data[top_n_features]
    #             )

    def select_features(self, data_materials: DataMaterials):
        self.load_shap_values(data_materials)
        # self.append_top_n_features(data_materials)



#
#
# class AggregatedFeatureSelector:                   # fixme
#     def __init__(self, features_list: List[List], aggregation_method=None):
#         self.aggregation_method = AggregationMethods.occurance
#
#
#
#     def aggregator(self):
#         if self.aggregation_method == AggregationMethods.occurance:
#
