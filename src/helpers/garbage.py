# def type_coercion_train_data(train_data_binned):
#     # Some columns have been interpreted as object type, eventhough they are actually numeric.
#     log.debug('{}'.format(set(train_data_binned.dtypes)))
#
#     # These non-numeric interpereted columns will be coerced.  NaN  values will be converted to  0 .
#     features = [column for column in train_data_binned.columns if
#                 column not in ['UniProt_ID', 'Mutation', 'Interactor_UniProt_ID']]
#
#     # Get column names where its type is *not* int or float, i.e. whose type is object.
#     coerce_numeric_cols = set(
#         [cname for cname in features if train_data_binned[cname].dtype not in ['int64', 'float64']]
#     )
#
#     # Remove target variable from the list
#     coerce_numeric_cols = coerce_numeric_cols - {"Mutation_Effect_Label", "UniProt_ID", "Mutation",
#                                                  "Interactor_UniProt_ID"}
#
#     for cname in coerce_numeric_cols:
#         train_data_binned[cname] = pd.to_numeric(train_data_binned[cname], errors='coerce')
#
#     train_data_binned.fillna(0, inplace=True)
#
#     # Now all columns are interpreted as numeric type, except "UniProt_ID", "Mutation", "Interactor_UniProt_ID".
#     log.debug('{}'.format(set(train_data_binned[features].dtypes)))
#     assert set(train_data_binned.dtypes[features].values) == {np.dtype('int64'), np.dtype('float64')}
#
#     data_processed = train_data_binned.copy()
#     return data_processed
#

#
# setattr(self, f"Xs_shap_{n_top}", [])
# setattr(self, f"ys_shap_{n_top}", [])
# setattr(self, f"Xs_train_shap_{n_top}", [])
# setattr(self, f"ys_train_shap_{n_top}", [])
# setattr(self, f"Xs_valid_shap_{n_top}", [])
# setattr(self, f"ys_valid_shap_{n_top}", [])

###############

from typing import List


# from ..Predator import Predator

#
# def initialize_data_materials_ML(predator):
#     predator.datasets["prepared_dataframes"] = []
#     predator.datasets["label_proportions_dataframes"] = []
#     predator.datasets["Xs"] = []
#     predator.datasets["ys"] = []
#     predator.datasets["Xs_train"] = []
#     predator.datasets["ys_train"] = []
#     predator.datasets["Xs_valid"] = []
#     predator.datasets["ys_valid"] = []
#     predator.datasets["Xs_train_random"] = []
#     predator.datasets["ys_train_random"] = []
#     predator.datasets["Xs_valid_random"] = []
#     predator.datasets["ys_valid_random"] = []
#
#
# def append_data_materials(predator, data_materials: dict):
#     predator.datasets["prepared_dataframes"].append(data_materials['data_prepared'])
#     predator.datasets["label_proportions_dataframes"].append(data_materials['label_proportions_data'])
#     predator.datasets["Xs"].append(data_materials['X'])
#     predator.datasets["ys"].append(data_materials['y'])
#     predator.datasets["Xs_train"].append(data_materials['X_train'])
#     predator.datasets["ys_train"].append(data_materials['y_train'])
#     predator.datasets["Xs_valid"].append(data_materials['X_valid'])
#     predator.datasets["ys_valid"].append(data_materials['y_valid'])
#     predator.datasets["Xs_train_random"].append(data_materials['X_train_random'])
#     predator.datasets["ys_train_random"].append(data_materials['y_train_random'])
#     predator.datasets["Xs_valid_random"].append(data_materials['X_valid_random'])
#     predator.datasets["ys_valid_random"].append(data_materials['y_valid_random'])
#
#
# def initialize_feature_selected_data_materials(predator, n_top: int):
#     predator.datasets[f"Xs_shap_{n_top}"] = []
#     predator.datasets[f"ys_shap_{n_top}"] = []
#     predator.datasets[f"Xs_train_shap_{n_top}"] = []
#     predator.datasets[f"ys_train_shap_{n_top}"] = []
#     predator.datasets[f"Xs_valid_shap_{n_top}"] = []
#     predator.datasets[f"ys_valid_shap_{n_top}"] = []
#
#
# class DataMaterials:
#     n_experiment = 0
#
#     def __init__(self):
#         self.prepared_dataframes = []
#         self.label_proportions_dataframes = []
#         self.Xs = []
#         self.ys = []
#         self.Xs_train = []
#         self.ys_train = []
#         self.Xs_valid = []
#         self.ys_valid = []
#         self.Xs_train_random = []
#         self.ys_train_random = []
#         self.Xs_valid_random = []
#         self.ys_valid_random = []
#         # self.
#
#     def append_data_materials(self, data_materials: dict):
#         self.prepared_dataframes.append(data_materials['data_prepared'])
#         self.label_proportions_dataframes.append(data_materials['label_proportions_data'])
#         self.Xs.append(data_materials['X'])
#         self.ys.append(data_materials['y'])
#         self.Xs_train.append(data_materials['X_train'])
#         self.ys_train.append(data_materials['y_train'])
#         self.Xs_valid.append(data_materials['X_valid'])
#         self.ys_valid.append(data_materials['y_valid'])
#         self.Xs_train_random.append(data_materials['X_train_random'])
#         self.ys_train_random.append(data_materials['y_train_random'])
#         self.Xs_valid_random.append(data_materials['X_valid_random'])
#         self.ys_valid_random.append(data_materials['y_valid_random'])
#         self.n_experiment += 1
#
#     def initialize_feature_selected_data_materials(self, n_top: int):
#         setattr(self, f"Xs_shap_{n_top}", [])
#         setattr(self, f"ys_shap_{n_top}", [])
#         setattr(self, f"Xs_train_shap_{n_top}", [])
#         setattr(self, f"ys_train_shap_{n_top}", [])
#         setattr(self, f"Xs_valid_shap_{n_top}", [])
#         setattr(self, f"ys_valid_shap_{n_top}", [])
#
#     def append_feature_selected_data_materials(self, n_top, selected_features):
#         for exp in range(self.n_experiment):
#             getattr(self, f"Xs_shap_{n_top}").append(
#                 self.Xs[exp][selected_features]
#             )
#
#             # setattr(self, f"Xs_shap_top_{n_top}", self.Xs[selected_features])
#             # setattr(self, f"ys_shap_top_{n_top}", self.ys[selected_features])
#             # setattr(self, f"Xs_train_shap_top_{n_top}", self.Xs_train[selected_features])
#             # setattr(self, f"ys_train_shap_top_{n_top}", self.ys_train[selected_features])
#             # setattr(self, f"Xs_valid_shap_top_{n_top}", self.Xs_valid[selected_features])
#             # setattr(self, f"ys_valid_shap_top_{n_top}", self.ys_valid[selected_features])
#
#
