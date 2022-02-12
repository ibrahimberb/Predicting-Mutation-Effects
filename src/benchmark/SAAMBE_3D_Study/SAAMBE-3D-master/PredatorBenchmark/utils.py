import urllib.request
import os
from datetime import datetime
import os.path as op


def mutation_effect_label_binner(train_data):
    """
    Mutation Effect label binning is only applied to train_data.
    Apply Label binning.
        - Disruptive → 0
        - No effect + Increasing → 1
        - Decreasing → dropped
        - Causing → dropped
    :param train_data:
    :return:
    """

    labels_to_bins = {
        "mutation disrupting(MI:0573)": 0,
        "mutation decreasing(MI:0119)": "IGNORED",
        "mutation disrupting strength(MI:1128)": 0,
        "mutation decreasing strength(MI:1133)": "IGNORED",
        "mutation with no effect(MI:2226)": 1,
        "disrupting": 0,
        "mutation increasing(MI:0382)": 1,
        "mutation increasing strength(MI:1132)": 1,
        "mutation decreasing rate(MI:1130)": "IGNORED",
        "mutation disrupting rate(MI:1129)": 0,
        "mutation causing(MI:2227)": "IGNORED",
        "mutation increasing rate(MI:1131)": 1,
    }

    replace_map = {"Mutation_Effect_Label": labels_to_bins}

    # Size of dataframe before binning.
    print("Size of dataframe before binning: {}".format(train_data.shape))

    # Modifications will be done on train_data_binned.
    train_data_binned = train_data.copy()

    # Replace the labels as described above.
    train_data_binned.replace(replace_map, inplace=True)

    # Drop the entries with "IGNORED": 'mutation cusing' in this case.
    train_data_binned = train_data_binned[
        train_data_binned["Mutation_Effect_Label"] != "IGNORED"
    ]

    # Reset index of the dataframe to avoid any possible errors
    train_data_binned.reset_index(drop=True, inplace=True)

    # Size of dataframe after binning.
    print("Size of dataframe after binning: {}".format(train_data_binned.shape))

    # First 5 rows of binned data.
    print("Train Data Binned:")
    print("Dataframe head: {}".format(train_data_binned.head()))

    # Confirming replacement of values are properly done. Mutation_Effect_Label only contains of 0 or 1.
    assert set(train_data_binned["Mutation_Effect_Label"].value_counts().index) == {0, 1}

    return train_data_binned


# download the pdb file from Protein Data Bank
def download_pdb(pdb_id):
    # print the status of the download
    print("Downloading " + pdb_id + "...")
    url = "http://www.rcsb.org/pdb/files/" + pdb_id + ".pdb"
    file_name = pdb_id + ".pdb"
    
    # create a folder named pdb_files if it doesn't exist
    if not os.path.exists("pdb_files"):
        os.makedirs("pdb_files")

    file_path = os.path.join("pdb_files", file_name)
        
    if not os.path.exists(file_path):
        urllib.request.urlretrieve(url, file_path)
        
    else:
        print("File already exists")
        return


def save_prediction_data(predator_benchmark_dir, prediction_file_name, prediction_data):
    file_date = datetime.now().strftime("%Y-%m-%d")
    prediction_file_name = "{}_{}.csv".format(prediction_file_name, file_date)
    prediction_data.to_csv(op.join(predator_benchmark_dir, prediction_file_name), index=False)
    print("Prediction data `{}`is exported.".format(op.join(predator_benchmark_dir, prediction_file_name)))


def save_errors(errors, file_path):
    with open(file_path, "w") as f:
        for error in errors:
            f.write("%s\n" % error)
            
    print("Errors written to file: " + file_path)
