# Aug 19th, 20201
# Evaluation Metrics

from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import cross_val_score
import pandas as pd
from tqdm.notebook import tqdm

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


def cross_val_score_feature_comparison(X, y, scoring, clf, n_repeats, n_jobs):
    # In calculation of scores, cross-validation is repeated n times, which yields a total of 10*n folds.
    # E.g. if n=10, it means cross-validation is repeated 10 times with a total of 100 folds.
    return (round(cross_val_score(clf, X, y,
                                  cv=RepeatedStratifiedKFold(n_splits=10, n_repeats=n_repeats),
                                  scoring=scoring, n_jobs=n_jobs).mean(), 4))


def evaluate_metric(X_benchmark_feature_names_dataframes, y, clf, metric, n_repeats, n_jobs, verbose):
    scores_comparison = []
    for X_item_name, X_item in X_benchmark_feature_names_dataframes:
        scores = cross_val_score_feature_comparison(X_item, y, metric, clf, n_repeats=n_repeats, n_jobs=n_jobs)
        scores_comparison.append(scores)
        if verbose:
            print("{: <28}: {}".format(X_item_name, scores))

    return scores_comparison


def evaluate_metrics(X_benchmark_feature_names_dataframes, y, clf, n_repeats, n_jobs, verbose, eval_metrics=None):
    if eval_metrics is None:
        eval_metrics = EVAL_METRICS

    scoring_metrics = {}

    for metric in eval_metrics:
        if verbose:
            print(F"\nEVALUATION METRIC: {metric.upper()}")
            print("------------------------------------")
        scores_comparison = evaluate_metric(X_benchmark_feature_names_dataframes, y, clf=clf, metric=metric,
                                            n_repeats=n_repeats,
                                            n_jobs=n_jobs, verbose=verbose)
        scoring_metrics[metric] = scores_comparison
        if verbose:
            print("====================================")

    scoring_metrics_table = pd.DataFrame(scoring_metrics, index=[e[0] for e in X_benchmark_feature_names_dataframes])

    return scoring_metrics_table
