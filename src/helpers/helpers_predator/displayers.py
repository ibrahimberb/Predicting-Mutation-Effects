import pandas as pd
from IPython.display import display
import seaborn as sns

from ..mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


# TODO: unittest
def display_label_counts(data_param):
    """
    Display a dataframe that contains label categories and their counts.
    """
    label_counts = pd.DataFrame(data_param["Mutation_Effect_Label"].value_counts())
    label_counts.reset_index(inplace=True)
    label_counts.columns = ["Mutation_Effect_Label", "Counts"]
    label_counts.rename(index={0: 'Disrupting', 1: 'Increasing + No Effect'}, inplace=True)
    display(label_counts)


# TODO: unittest
def display_labels(data_param):
    """
    Display a dataframe that contains label categories.
    """
    label_counts = pd.DataFrame(data_param["Mutation_Effect_Label"].value_counts().index)
    label_counts.columns = ["Mutation_Effect_Label"]
    display(label_counts)


def visualize_label_counts(data_param, label_name_param="Mutation_Effect_Label"):
    # noinspection SpellCheckingInspection
    sns.set(style="white", font_scale=1.15)  # white, dark, whitegrid, darkgrid, ticks
    val_counts = data_param[label_name_param].value_counts().sort_index()
    val_counts = val_counts.rename({0: 'Disrupting', 1: 'Increasing + No Effect'})
    log.debug(f"Label counts:\n{val_counts}")
    ax = sns.barplot(x=val_counts.index,
                     y=val_counts,
                     palette="ch:s=-.2,r=.6")
    ax.set_title('Disrupting vs Increasing & No Effect')  # ch:s=-.2,r=.6, ocean
    ax.set_ylabel('Value counts')
