# Aug 19th, 20201
# Evaluation Metrics

from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import cross_val_score
import pandas as pd
from tqdm.notebook import tqdm
from .machine_learning_utils import get_default_classifier

EVAL_METRICS = [
    "f1",
    "balanced_accuracy",
    "accuracy",
    "f1_macro",
    "f1_micro",
    "precision",
    "recall",
    "roc_auc",
    "precision_macro",
    "precision_micro"

]


# class EvaluationMetrics:
#     def __init__(self, n_experiment):
#         self.n_experiment = n_experiment
#         self.Xs_train_benchmark_feature_names_dataframes_list = None
#
#     def initialize_benchmark_dataframes(self, data_materials):
#         Xs_train_benchmark_feature_names_dataframes_list = []
#         for exp in range(self.n_experiment):
#             benchmark_feature_names_dataframes = [
#                 (f"Xs_train_{exp}_provean", data_materials["Xs_train"][exp][['Provean_score']]),
#                 (f"Xs_train_{exp}_ddG", data_materials["Xs_train"][exp][['Final_ddG']]),
#                 (f"Xs_train_{exp}", data_materials["Xs_train"][exp]),
#                 (f"Xs_train_{exp}_shap_HSF_10", data_materials["Xs_train"][exp][highly_selected_10_features]),
#             ]
#             Xs_train_benchmark_feature_names_dataframes_list.append(benchmark_feature_names_dataframes)
#         self.Xs_train_benchmark_feature_names_dataframes_list = Xs_train_benchmark_feature_names_dataframes_list


def cross_val_score_feature_comparison(X, y, scoring, n_repeats, n_jobs):
    # In calculation of scores, cross-validation is repeated n times, which yields a total of 10*n folds.
    # E.g. if n=10, it means cross-validation is repeated 10 times with a total of 100 folds.
    clf = get_default_classifier(random_state=42)
    return (round(cross_val_score(clf, X, y,
                                  cv=RepeatedStratifiedKFold(n_splits=10, n_repeats=n_repeats),
                                  scoring=scoring, n_jobs=n_jobs).mean(), 4))


def evaluate_metric(X_benchmark_feature_names_dataframes: dict, y, metric, n_repeats, n_jobs, verbose):
    scores_comparison = []
    for X_item_name, X_item in X_benchmark_feature_names_dataframes.items():
        scores = cross_val_score_feature_comparison(X_item, y, metric, n_repeats=n_repeats, n_jobs=n_jobs)
        scores_comparison.append(scores)
        if verbose:
            print("{: <28}: {}".format(X_item_name, scores))

    return scores_comparison


# def evaluate_metric(X_benchmark_feature_names_dataframes, y, clf, metric, n_repeats, n_jobs, verbose):
#     scores_comparison = []
#     for X_item_name, X_item in X_benchmark_feature_names_dataframes:
#         scores = cross_val_score_feature_comparison(X_item, y, metric, clf, n_repeats=n_repeats, n_jobs=n_jobs)
#         scores_comparison.append(scores)
#         if verbose:
#             print("{: <28}: {}".format(X_item_name, scores))
#
#     return scores_comparison


def evaluate_metrics(X_benchmark_feature_names_dataframes, y, n_repeats, n_jobs, verbose, eval_metrics=None):
    if eval_metrics is None:
        eval_metrics = EVAL_METRICS

    scoring_metrics = {}

    for metric in eval_metrics:
        if verbose:
            print(F"\nEVALUATION METRIC: {metric.upper()}")
            print("------------------------------------")
        scores_comparison = evaluate_metric(X_benchmark_feature_names_dataframes, y, metric=metric,
                                            n_repeats=n_repeats,
                                            n_jobs=n_jobs, verbose=verbose)
        scoring_metrics[metric] = scores_comparison
        if verbose:
            print("====================================")

    scoring_metrics_table = pd.DataFrame(scoring_metrics, index=[feature_names for
                                                                 feature_names in X_benchmark_feature_names_dataframes])

    return scoring_metrics_table

# def evaluate_metrics(X_benchmark_feature_names_dataframes, y, clf, n_repeats, n_jobs, verbose, eval_metrics=None):
#     if eval_metrics is None:
#         eval_metrics = EVAL_METRICS
#
#     scoring_metrics = {}
#
#     for metric in eval_metrics:
#         if verbose:
#             print(F"\nEVALUATION METRIC: {metric.upper()}")
#             print("------------------------------------")
#         scores_comparison = evaluate_metric(X_benchmark_feature_names_dataframes, y, clf=clf, metric=metric,
#                                             n_repeats=n_repeats,
#                                             n_jobs=n_jobs, verbose=verbose)
#         scoring_metrics[metric] = scores_comparison
#         if verbose:
#             print("====================================")
#
#     scoring_metrics_table = pd.DataFrame(scoring_metrics, index=[e[0] for e in X_benchmark_feature_names_dataframes])
#
#     return scoring_metrics_table
