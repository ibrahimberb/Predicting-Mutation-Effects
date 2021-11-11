from .loaders import load_snv_datasets
from ..mylogger import get_handler
import logging
from tqdm.auto import tqdm

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class MutualExclusivity:
    def __init__(
            self,
            tcga,
            tcga_snv_path,
    ):
        self.tcga = tcga.upper()
        self.tcga_snv_path = tcga_snv_path
        self.snv_data = None
        self.patients = None
        self.patient_to_snv_data = None

        self.load_snv_data()
        self.load_patient_ids()
        self.load_patient_to_snv_data()

    def load_snv_data(self):
        log.debug("Loading SNV data simplified ..")
        data_materials = {}
        load_snv_datasets(self.tcga, self.tcga_snv_path, data_materials)
        self.snv_data = data_materials[f"{self.tcga}_snv_data_simplified"]

    def load_patient_ids(self):
        log.debug("Loading patient ids ..")
        patients = list(self.snv_data["Tumor_Sample_Barcode"].unique())
        self.patients = patients

    def get_patient_snv_data(self, patient_id):
        patient_snv_data = self.snv_data[
            self.snv_data["Tumor_Sample_Barcode"] == patient_id
            ]
        return patient_snv_data

    def load_patient_to_snv_data(self):
        log.debug("Loading patient to snv_data ..")
        patient_to_snv_data = {}
        for patient in tqdm(self.patients):
            patient_snv_data = self.get_patient_snv_data(patient)
            patient_to_snv_data[patient] = patient_snv_data

        self.patient_to_snv_data = patient_to_snv_data

    def get_patients_with(self, identifier, identifier_type):

        if identifier_type == "protein":
            patients = self.snv_data[
                self.snv_data["SWISSPROT"] == identifier
                ]["Tumor_Sample_Barcode"].unique()

            return list(patients)

        elif identifier_type == "gene":
            raise NotImplementedError

        else:
            raise ValueError(f"Invalid argument `{identifier_type}` for `identifier_type`.")

    # def get_patients_with(self, identifier, identifier_type):
    #
    #     if identifier_type == "gene":
    #         raise NotImplementedError
    #
    #     elif identifier_type == "protein":
    #         patients_with_input_protein = []
    #         for patient in self.patients:
    #             if is_identifier_occur_in_patient(identifier, patient):
    #                 patients_with_input_protein.append(patient)
    #
    #         log.debug(f"{} patients found for input {identifier_type} in {self.tcga} SNV files.")
    #         return patients_with_input_protein
    #
    #
    #     else:
    #         raise ValueError(f"Invalid argument `{identifier_type}` for `identifier_type`.")

    def calculate_mutual_exclusivity(self, identifier_1, identifier_2):
        # todo: using genes to SNV data may not be a good idea
        #  an alternative is using Protein, and retrieving gene id of that protein, and then querying in SNV data.
        s1 = set(self.get_patients_with(identifier_1, "protein"))
        s2 = set(self.get_patients_with(identifier_2, "protein"))

        print(f"S1: {s1} ({len(s1)} patients)")
        print(f"S2: {s2} ({len(s2)} patients)")

        """
        S1 ERBB'nun mutasyona uğradığı hasta seti olsun. 
        S2 SRPK1'nın mutasyona uğradığı hasta seti olsun. 
        
        |S1 union S2| / |S1| + |S2| değerini hesaplayalım. 
        """

        mutual_exclusivity_value = len(s1.union(s2)) / (len(s1) + len(s2))

        return mutual_exclusivity_value
