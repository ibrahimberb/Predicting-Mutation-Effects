import pandas as pd
from pathlib import Path
import os
import os.path as op
from datetime import datetime
from src.helpers.helpers_analysis.gene_id_retrieval import GeneIDFetcher


def get_sheet_names(excel_path):
    excel = pd.ExcelFile(excel_path)
    excel_sheet_names = excel.sheet_names  # see all sheet names
    return excel_sheet_names


def excel_to_csv_cell_line(file_name, excel_path, sheet_name, folder_name="Parsed", use_different_name=None):
    folder_path = op.join(os.curdir, folder_name)
    Path(f"{folder_path}").mkdir(parents=True, exist_ok=True)

    data = pd.read_excel(excel_path, sheet_name=sheet_name)

    if use_different_name is None:
        file_name = f"{file_name}_{sheet_name}.csv"
    else:
        file_name = f"{file_name}_{use_different_name}.csv"
    file_path = op.join(folder_path, file_name)
    if op.isfile(file_path):
        raise FileExistsError(f"You already have the file {file_path}")

    data.to_csv(file_path)
    print(f"{file_name} Exported successfully.")


class StandardizedColumnNames:
    SELF = "SELF_PROTEIN"
    INTERACTOR = "INTERACTOR_PROTEIN"


def standardize_file(
        data_path,
        self_protein_col_name,
        interactor_protein_col_name,
        parsed_path="Parsed",
        standardized_path="ParsedStandardized"
):
    data = pd.read_csv(fr"{parsed_path}\{data_path}", index_col=0)
    standardized_data = data[[self_protein_col_name, interactor_protein_col_name]].copy()
    standardized_data.columns = [StandardizedColumnNames.SELF, StandardizedColumnNames.INTERACTOR]
    output_file_name = fr"{standardized_path}\{data_path}"
    if op.isfile(output_file_name):
        raise FileExistsError("You already have the file!")

    standardized_data.to_csv(output_file_name, index=False)
    print(f"{data_path} is exported successfully.")


def find_in_single_standardized_data(standardized_data, protein, interactor):
    query = standardized_data[
        (standardized_data[StandardizedColumnNames.SELF] == protein) &
        (standardized_data[StandardizedColumnNames.INTERACTOR] == interactor)
        ]

    if query.empty:
        return False

    return True


def search_in_parsed_standardized(protein, interactor, standardized_path="ParsedStandardized"):
    files = os.listdir(standardized_path)
    found_file_names = []
    for file in files:
        data_path = op.join(standardized_path, file)
        std_data = pd.read_csv(data_path)
        if find_in_single_standardized_data(std_data, protein, interactor):
            found_file_names.append(file)

    return found_file_names


class FindInScienceArticles:
    def __init__(self, tcga, standardized_path):
        self.tcga = tcga.upper()
        self.entries = []
        self.standardized_path = standardized_path
        self.data = None
        self.gene_id_fetcher = GeneIDFetcher()

    def look_up(self, protein, interactor):
        found_files = search_in_parsed_standardized(protein, interactor, self.standardized_path)
        found_flag = bool(len(found_files))
        self.entries.append(
            (
                self.tcga,
                protein,
                self.gene_id_fetcher.fetch(protein),
                interactor,
                self.gene_id_fetcher.fetch(interactor),
                found_flag,
                found_files
            )
        )

    def construct_table(self):
        data = pd.DataFrame(
            self.entries, columns=[
                "TCGA", "PROTEIN", "GENE", "INTERACTOR_PROTEIN", "INTERACTOR_GENE", "FOUND_FLAG", "FOUND_FILES"
            ]
        )
        self.data = data
        print("Table constructed.")

    def export_summary_table(self, file_path, file_name="FindInScienceArticles"):
        if self.data is None:
            self.construct_table()

        file_date = datetime.today().strftime('%Y-%m-%d')
        file_name = f"{file_name}_{self.tcga}_{file_date}.csv"
        file_name = op.join(file_path, file_name)
        if op.isfile(file_name):
            raise FileExistsError

        self.data.to_csv(file_name, index=False)
        print("Table extracted.")
