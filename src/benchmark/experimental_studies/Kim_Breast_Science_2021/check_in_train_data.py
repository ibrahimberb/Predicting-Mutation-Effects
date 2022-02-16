# Check if protein.mutation pairs are already present in training dataset.

import pandas as pd


TRAIN_DATA_PATH = r"../../../"


from src.helpers.helpers_analysis.convert_primary_isomer import convert_primary_isomer

train_data_740 = pd.read_csv("../../../processed_data_740.csv")
triplets_data = train_data_740[["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]].copy()
triplets_data = convert_primary_isomer("Interactor_UniProt_ID", triplets_data)
triplets_data = convert_primary_isomer("UniProt_ID", triplets_data)

triplets = set(
    zip(
        triplets_data["UniProt_ID"], triplets_data["Mutation"], triplets_data["Interactor_UniProt_ID"]
    )
)

print(triplets_data.head())
print(len(triplets))
print(triplets)

unique_proteins = sorted(triplets_data["UniProt_ID"].unique())
print(f"{len(unique_proteins)=}")
print(unique_proteins)

# def load_pairs(file_path):
#     with open(file_path) as fin:
#         pairs = [line.strip() for line in fin.readlines()]
#         pairs = [tuple(pair.split('.')) for pair in pairs]
#
#     pairs = set(pairs)
#
#     return pairs
#
#
# pairs_1 = load_pairs("ELASPIC_Input/S3_2021-11-19.txt")
# pairs_2 = load_pairs("ELASPIC_Input/S4_2021-11-19.txt")
#
# print(f"{len(pairs_1)=}")
# print(pairs_1)
#
# print(f"{len(pairs_2)=}")
# print(pairs_2)
#
#
# all_pairs = pairs_1.union(pairs_2)


# check if in training data
# why am I doing this?
