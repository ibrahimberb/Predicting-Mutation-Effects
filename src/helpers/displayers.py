import pandas as pd
from IPython.display import display
import seaborn as sns


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


# TODO: unittest
def visualize_label_counts(data_param, label_name_param="Mutation_Effect_Label"):
    sns.set(style="white", font_scale=1.15)  # white, dark, whitegrid, darkgrid, ticks
    ax = sns.barplot(x=data_param[label_name_param].value_counts().index,
                     y=data_param[label_name_param].value_counts(),
                     palette="ch:s=-.2,r=.6")
    ax.set_title('Disrupting vs Increasing & No Effect')  # ch:s=-.2,r=.6, ocean
    ax.set_ylabel('Value counts')
    ax.set_xticklabels(['Distrupting', 'Increasing + No Effect'])
