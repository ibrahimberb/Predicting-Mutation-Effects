from typing import List

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from .mylogger import get_handler
import logging

handler_simple = get_handler('simple')
log_simple = logging.getLogger('Visualizers')
log_simple.handlers[:] = []
log_simple.addHandler(handler_simple)
log_simple.setLevel(logging.DEBUG)


def get_sampled_datasets_label_counts(sampled_train_data_list):
    sampled_train_data_to_label_counts = {}
    for i, sampled_train_data in enumerate(sampled_train_data_list):
        label_counts = sampled_train_data["Mutation_Effect_Label"].value_counts()
        sampled_train_data_to_label_counts["SAMPLED_TRAIN_DATA_" + str(i + 1)] = [
            label_counts.loc[0],
            label_counts.loc[1],
        ]
    sampled_datasets_label_counts = pd.DataFrame(sampled_train_data_to_label_counts).T
    sampled_datasets_label_counts.columns = ["Disrupting", "Increasing+NoEff"]
    sampled_datasets_label_counts.index.name = "SAMPLED_TRAIN_DATA"
    sampled_datasets_label_counts.reset_index(inplace=True)
    return sampled_datasets_label_counts


def visualize_sampled_train_datasets_label_counts(sampled_train_data_list, kind):

    sampled_datasets_label_counts = get_sampled_datasets_label_counts(
        sampled_train_data_list
    )

    if kind == "strip":

        experiment_statistics_data_melted = pd.melt(
            sampled_datasets_label_counts,
            id_vars=["SAMPLED_TRAIN_DATA"],
            value_vars=["Disrupting", "Increasing+NoEff"],
            var_name="MUTATION_EFFECT",
            value_name="LABEL_COUNT",
        )

        plt.figure(figsize=(3, 4))
        sns.stripplot(
            x="MUTATION_EFFECT",
            y="LABEL_COUNT",
            data=experiment_statistics_data_melted,
            palette="ch:s=-.2,r=.6",
            jitter=True,
        )

    elif kind == "bar":
        sampled_datasets_label_counts.plot(
            figsize=(25, 4), kind="bar", color=["#E3D9C1", "#27213F"],
            title='Label Counts per Experiment', xlabel='Experiment', ylabel='Counts',
            rot=0
        )

    else:
        log_simple.error(f"Parameter `kind` must be either `strip` or `bar`, not `{kind}`")


# def visualize_sampled_train_datasets_label_counts_stripplot(sampled_train_data_list):
#     sampled_datasets_label_counts = get_sampled_datasets_label_counts(
#         sampled_train_data_list
#     )
#
#     experiment_statistics_data_melted = pd.melt(
#         sampled_datasets_label_counts,
#         id_vars=["SAMPLED_TRAIN_DATA"],
#         value_vars=["Disrupting", "Increasing+NoEff"],
#         var_name="MUTATION_EFFECT",
#         value_name="LABEL_COUNT",
#     )
#
#     plt.figure(figsize=(3, 4))
#     sns.stripplot(
#         x="MUTATION_EFFECT",
#         y="LABEL_COUNT",
#         data=experiment_statistics_data_melted,
#         palette="ch:s=-.2,r=.6",
#         jitter=True,
#     )


# def visualize_sampled_train_datasets_label_counts_barplot(sampled_train_data_list):
#     sampled_datasets_label_counts = get_sampled_datasets_label_counts(
#         sampled_train_data_list
#     )
#     sampled_datasets_label_counts.plot(
#         figsize=(25, 4), kind="bar", color=["#E3D9C1", "#27213F"],
#         title='Label Counts per Experiment', xlabel='Experiment', ylabel='Counts',
#         rot=0
#     )


def visualize_accuracy_metrics(
    acc_scores: List[float],
    balan_acc_scores: List[float],
    kind='strip'
):
    data = pd.DataFrame({'ACCURACY': acc_scores,
                         'BALANCED_ACCURACY': balan_acc_scores})
    data_melted = pd.melt(data, var_name='METRIC', value_name='SCORE')
    # plt.figure(figsize=(4, 5))
    sns.catplot(x='METRIC', y='SCORE', data=data_melted,
                palette=sns.color_palette(['#265191', '#9F2945']),
                kind=kind)
    plt.show()
