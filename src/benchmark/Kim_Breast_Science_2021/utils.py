from src.benchmark.benchmark_utils import (
    SupplementaryExcelParser,
    get_mutation_corrected
)
from src.helpers.mylogger import get_handler
import re

import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def is_valid_aa_change_in_experiment_id(item) -> bool:
    """
    Given items in the <gene>_<mutant> form, check whether mutant is in
    <letter><position><letter> form.
    """

    if "_" not in item:
        return False

    valid_form = re.match(r"(.+)_([A-Za-z][0-9]+[A-Za-z]_?[0-9]*)$", item)
    if not valid_form:
        return False

    return True


class Kim2021Table(SupplementaryExcelParser):
    def __init__(
            self,
            supp_excel_path,
            sheet_range=None,
            protein_col="Bait ID",
            mutation_col="Experiment.ID"
    ):
        super().__init__(supp_excel_path)

        self.pairs = set()
        self.protein_col = protein_col
        self.mutation_col = mutation_col

        self.sheet_names = self.get_sheet_names(supp_excel_path)
        log.info(f"Sheets found ({len(self.sheet_names)}): {self.sheet_names}")

        if sheet_range is None:
            # Sheet 0 is README
            self.sheet_names = self.sheet_names[1:]
        else:
            log.info(f"Sheets range {sheet_range[0]} - {sheet_range[1]}")
            self.sheet_names = self.sheet_names[sheet_range[0]: sheet_range[1]]

        log.info(f"Sheets selected ({len(self.sheet_names)}): {self.sheet_names}")

        for sheet in self.sheet_names:
            sheet_data = self.load_excel(supp_excel_path, sheet_name=sheet)
            self.pairs.update(self.get_inputs(sheet_data))

        print(f"Number of pairs: {len(self.pairs)}")
        print(self.pairs)

    def get_inputs(self, data):
        # Drop na values in Experiment.ID column
        data = data.dropna(subset=["Experiment.ID"])

        # Get entries where missense mutation is provided.
        mutation_available_data = data[
            (data["Experiment.ID"].apply(lambda x: is_valid_aa_change_in_experiment_id(x)))
        ]
        protein_series = mutation_available_data[self.protein_col]
        mutation_series = mutation_available_data[self.mutation_col].apply(lambda x: x.split('_')[1])
        mutation_series = mutation_series.apply(lambda x: get_mutation_corrected(x))
        pairs = self.get_input_pairs_via_series(
            protein_series, mutation_series
        )

        return pairs

# class Kim2021TableS3(Kim2021TableTemplate):
#     def __init__(self, supp_excel_path, sheet_range=None):
#         log.info(f"Initializing S3 ..")
#         super().__init__(supp_excel_path, sheet_range)
#
#
# class Kim2021TableS4(Kim2021TableTemplate):
#     def __init__(self, supp_excel_path, sheet_range=None):
#         log.info(f"Initializing S4 ..")
#         super().__init__(supp_excel_path, sheet_range)


# s3 = Kim2021TableS3(r"science.abf3066_tables_s2_to_s12\science.abf3066_Table_S3.xlsx")
# print("S3 PAIRS:")
# print(s3.pairs)

# s3.extract_pairs("test.txt")

# L = {('P42336', 'h1047r'), ('P42336', 'e545k'), ('Q86UE4', 'wt'), ('P04637', 'wt'), ('P38398', 'iso5'),
#      ('Q9BX63', 'a745t'), ('P38398', 'c61g'), ('P42336', 'm1043v'), ('P04637', 'r175h'), ('P04637', 'r273h'),
#      ('Q86YC2', 'e837k'), ('P31749', 'e17k'), ('P38398', '5382insC'), ('O96017', 'k373e'), ('P12830', 'e243k'),
#      ('P42336', 'wt'), ('P04637', 'r248w'), ('Q9BX63', 'wt'), ('P38398', 's1655f'), ('Q86UE4', 'a78s'),
#      ('P31751', 'wt'), ('O96017', 'wt'), ('Q9Y243', 'e17k'), ('P01112', 'g12d'), ('Q86YC2', 'wt'), ('P60484', 'r130q'),
#      ('P31749', 'wt'), ('P38398', 'r71g'), ('P01112', 'wt'), ('P38398', 'wt'), ('P38398', 'm1775r'), ('P38398', 'i26a'),
#      ('P12830', 'wt'), ('O96017', '1100delC'), ('P60484', 'wt'), ('Q9Y243', 'wt')}
#
# print("asserting L == s3 (with the same items)")
# # assert L == s3.pairs
# assert L == s3.pairs.union(L)

# for L_item in L:
#     protein = L_item[0]
#     mutation = L_item[1]
#     val = f"{protein}_{mutation}"
#
#     print(f"VAL: {val}")
#     print(is_valid_aa_change_in_experiment_id(val))
#     if is_valid_aa_change_in_experiment_id(val):
#         print(f"CORRECTED: {get_mutation_corrected(mutation)}")
#
#     print()
