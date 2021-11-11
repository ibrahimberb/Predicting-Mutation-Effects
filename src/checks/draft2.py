from collections import defaultdict
from typing import List

from tqdm import tqdm
import pandas as pd
from IPython.display import display
import pickle

from src.helpers.helpers_analysis.get_patient_protein_to_mutations_dict import get_patient_protein_to_mutations_dict
from src.helpers.helpers_analysis.preprocessing import process_snv
from src.helpers.helpers_analysis.is_core import is_core
from src.helpers.helpers_analysis.is_in_elaspic import is_in_elaspic

# Paths
ELASPIC_RESULTS_COMMON_PATH = "../../../My-ELASPIC-Web-API/Elaspic_Results/Merged_Results/"  # elaspic_results_datasets
BRCA_CORE_PATH = ELASPIC_RESULTS_COMMON_PATH + "BRCA_Core_2021-09-28.txt"
BRCA_INTERFACE_PATH = ELASPIC_RESULTS_COMMON_PATH + "BRCA_Interface_2021-09-28.txt"
SNV_COMMON_PATH = "C:/Users/ibrah/Desktop/TUSEB_Study/Data_Collection_and_Filtering/SNV/"
SNV_BRCA_PATH = SNV_COMMON_PATH + "SNV_BRCA_hg38.csv"

# Read SNV file
brca_snv = pd.read_csv(SNV_BRCA_PATH, low_memory=False)
print(brca_snv.shape)
display(brca_snv.head(3))
print("======================================================================")

# Read core data
brca_core_data_original = pd.read_csv(BRCA_CORE_PATH, sep='\t', low_memory=False)
print(brca_core_data_original.shape)
display(brca_core_data_original.head(3))
print("======================================================================")

# Read interface data
brca_interface_data_original = pd.read_csv(BRCA_INTERFACE_PATH, sep='\t', low_memory=False)
print(brca_interface_data_original.shape)
display(brca_interface_data_original.head(3))
print("======================================================================")

# SNV Ready
brca_snv_processed = process_snv(brca_snv)
brca_snv_simplified = brca_snv_processed[["Hugo_Symbol", "SWISSPROT", "HGVSp_Short", "Tumor_Sample_Barcode"]].copy()
brca_snv_simplified["Tumor_Sample_Barcode"] = brca_snv_simplified["Tumor_Sample_Barcode"].apply(lambda x: x[:12])
print(brca_snv_simplified.shape)
display(brca_snv_simplified.head(3))
print("======================================================================")

# Actual BRCA Patients
brca_patients_actual = sorted(brca_snv_simplified["Tumor_Sample_Barcode"].unique())
brca_patients_actual = [patient[:12] for patient in brca_patients_actual]
print(f"Number of BRCA patients: {len(brca_patients_actual)}")
display(brca_patients_actual[:3])
print("======================================================================")


def get_patient_snv(patient, snv_data_param):
    patient_snv = snv_data_param[
        snv_data_param["Tumor_Sample_Barcode"] == patient
        ]

    return patient_snv


def is_protein_and_its_mutations_interface_only(
        protein: str,
        mutations: List[str],
        core_data=brca_core_data_original,
        interface_data=brca_interface_data_original
):
    protein_interface_only_flag = False
    for mutation in mutations:
        if is_in_elaspic(protein, mutation, core_data, interface_data):
            if is_core(protein, mutation, core_data):
                return False

            else:
                protein_interface_only_flag = True

    return protein_interface_only_flag


def actualization(patients, snv_data, core_data, interface_data):
    proteins_to_num_interface_only_patients = defaultdict(int)

    patient_protein_mutation_to_C_I_status = dict()

    for patient in tqdm(patients):
        # print(f"Patient: {patient}")
        patient_snv = get_patient_snv(patient, snv_data)
        patient_protein_to_mutations = get_patient_protein_to_mutations_dict(patient_snv)

        for protein, mutations in patient_protein_to_mutations.items():

            interface_only = is_protein_and_its_mutations_interface_only(
                protein, mutations, core_data, interface_data
            )

            if interface_only:
                proteins_to_num_interface_only_patients[protein] += 1
                patient_protein_mutation_to_C_I_status[(protein, patient)] = "I"
                if protein == "P04637":  # ["Q9Y616", "Q9C0D5", P04637]:
                    print(f" - patient adding -- {patient} - {protein}.{mutations}")

            else:
                patient_protein_mutation_to_C_I_status[(protein, patient)] = "C"

    with open("brca_patient_protein_mutation_to_C_I_status.pickle", 'wb') as file_output:
        pickle.dump(patient_protein_mutation_to_C_I_status, file_output)
        print("patient_protein_mutation_to_C_I_status is exported.")

    return proteins_to_num_interface_only_patients


counts = actualization(
    brca_patients_actual, brca_snv_simplified, brca_core_data_original, brca_interface_data_original
)

print(counts)
# print("P04637 (TP53): {}".format(counts["P04637"]))
# print("Q9C0D5 (TANC1): {}".format(counts["Q9C0D5"]))
# print("A0AVT1 (UBA6): {}".format(counts["A0AVT1"]))
#
# print("A0AVT1 (UBA6): {}".format(counts["A0AVT1"]))
# print("A0AVT1 (UBA6): {}".format(counts["A0AVT1"]))
# print("A0AVT1 (UBA6): {}".format(counts["A0AVT1"]))
# print("A0AVT1 (UBA6): {}".format(counts["A0AVT1"]))
# print("A0AVT1 (UBA6): {}".format(counts["A0AVT1"]))
#
# print("P04637 (TP53): {} \t {}".format(counts["P04637"], counts["P04637"] == 55))
# print("P04626 (ERBB2): {} \t {}".format(counts["P04626"], counts["P04626"] == 14))
# print("P62805 (H4C1): {} \t {}".format(counts["P62805"], counts["P62805"] == 11))
# print("P62807 (H2BC4): {} \t {}".format(counts["P62807"], counts["P62807"] == 11))
# print("P68431 (H3C1): {} \t {}".format(counts["P68431"], counts["P68431"] == 11))
# print("P45985 (MAP2K4): {} \t {}".format(counts["P45985"], counts["P45985"] == 8))
# print("Q13951 (CBFB): {} \t {}".format(counts["Q13951"], counts["Q13951"] == 8))
# print("O75582 (RPS6KA5): {} \t {}".format(counts["O75582"], counts["O75582"] == 6))
# print("P0CG48 (UBC): {} \t {}".format(counts["P0CG48"], counts["P0CG48"] == 6))
# print("P12830 (CDH1): {} \t {}".format(counts["P12830"], counts["P12830"] == 6))
# print("P21860 (ERBB3): {} \t {}".format(counts["P21860"], counts["P21860"] == 6))
# print("A0AVT1 (UBA6): {} \t {}".format(counts["A0AVT1"], counts["A0AVT1"] == 5))
# print("Q14790 (CASP8): {} \t {}".format(counts["Q14790"], counts["Q14790"] == 5))
# print("Q8TBB1 (LNX1): {} \t {}".format(counts["Q8TBB1"], counts["Q8TBB1"] == 5))
# print("O00522 (KRIT1): {} \t {}".format(counts["O00522"], counts["O00522"] == 4))
# print("O14936 (CASK): {} \t {}".format(counts["O14936"], counts["O14936"] == 4))
# print("O15111 (CHUK): {} \t {}".format(counts["O15111"], counts["O15111"] == 4))
# print("O95835 (LATS1): {} \t {}".format(counts["O95835"], counts["O95835"] == 4))
# print("O96017 (CHEK2): {} \t {}".format(counts["O96017"], counts["O96017"] == 4))
# print("P00533 (EGFR): {} \t {}".format(counts["P00533"], counts["P00533"] == 4))
# print("P01112 (HRAS): {} \t {}".format(counts["P01112"], counts["P01112"] == 4))
# print("P17612 (PRKACA): {} \t {}".format(counts["P17612"], counts["P17612"] == 4))
# print("P42684 (ABL2): {} \t {}".format(counts["P42684"], counts["P42684"] == 4))
# print("P54646 (PRKAA2): {} \t {}".format(counts["P54646"], counts["P54646"] == 4))
# print("P61586 (RHOA): {} \t {}".format(counts["P61586"], counts["P61586"] == 4))
# print("Q05397 (PTK2): {} \t {}".format(counts["Q05397"], counts["Q05397"] == 4))
# print("Q08881 (ITK): {} \t {}".format(counts["Q08881"], counts["Q08881"] == 4))
# print("Q13043 (STK4): {} \t {}".format(counts["Q13043"], counts["Q13043"] == 4))
# print("Q13620 (CUL4B): {} \t {}".format(counts["Q13620"], counts["Q13620"] == 4))
# print("Q16778 (H2BC21): {} \t {}".format(counts["Q16778"], counts["Q16778"] == 4))
# print("Q92736 (RYR2): {} \t {}".format(counts["Q92736"], counts["Q92736"] == 4))
# print("Q9NZM3 (ITSN2): {} \t {}".format(counts["Q9NZM3"], counts["Q9NZM3"] == 4))
# print("A9YTQ3 (AHRR): {} \t {}".format(counts["A9YTQ3"], counts["A9YTQ3"] == 3))
# print("O00311 (CDC7): {} \t {}".format(counts["O00311"], counts["O00311"] == 3))
# print("O00571 (DDX3X): {} \t {}".format(counts["O00571"], counts["O00571"] == 3))
# print("O14802 (POLR3A): {} \t {}".format(counts["O14802"], counts["O14802"] == 3))
# print("O14980 (XPO1): {} \t {}".format(counts["O14980"], counts["O14980"] == 3))
# print("O15144 (ARPC2): {} \t {}".format(counts["O15144"], counts["O15144"] == 3))
# print("O43390 (HNRNPR): {} \t {}".format(counts["O43390"], counts["O43390"] == 3))
# print("O43747 (AP1G1): {} \t {}".format(counts["O43747"], counts["O43747"] == 3))
# print("O60285 (NUAK1): {} \t {}".format(counts["O60285"], counts["O60285"] == 3))
# print("O75306 (NDUFS2): {} \t {}".format(counts["O75306"], counts["O75306"] == 3))
# print("P07948 (LYN): {} \t {}".format(counts["P07948"], counts["P07948"] == 3))
# print("P08069 (IGF1R): {} \t {}".format(counts["P08069"], counts["P08069"] == 3))
# print("P17655 (CAPN2): {} \t {}".format(counts["P17655"], counts["P17655"] == 3))
# print("P20810 (CAST): {} \t {}".format(counts["P20810"], counts["P20810"] == 3))
# print("P23443 (RPS6KB1): {} \t {}".format(counts["P23443"], counts["P23443"] == 3))
# print("P23469 (PTPRE): {} \t {}".format(counts["P23469"], counts["P23469"] == 3))
# print("P24821 (TNC): {} \t {}".format(counts["P24821"], counts["P24821"] == 3))
# print("P28827 (PTPRM): {} \t {}".format(counts["P28827"], counts["P28827"] == 3))
# print("P33981 (TTK): {} \t {}".format(counts["P33981"], counts["P33981"] == 3))
# print("P36897 (TGFBR1): {} \t {}".format(counts["P36897"], counts["P36897"] == 3))
# print("P38117 (ETFB): {} \t {}".format(counts["P38117"], counts["P38117"] == 3))
# print("P42680 (TEC): {} \t {}".format(counts["P42680"], counts["P42680"] == 3))
# print("P46663 (BDKRB1): {} \t {}".format(counts["P46663"], counts["P46663"] == 3))
# print("P49356 (FNTB): {} \t {}".format(counts["P49356"], counts["P49356"] == 3))
# print("P51170 (SCNN1G): {} \t {}".format(counts["P51170"], counts["P51170"] == 3))
# print("P51812 (RPS6KA3): {} \t {}".format(counts["P51812"], counts["P51812"] == 3))
# print("P51813 (BMX): {} \t {}".format(counts["P51813"], counts["P51813"] == 3))
# print("P51817 (PRKX): {} \t {}".format(counts["P51817"], counts["P51817"] == 3))
# print("P53350 (PLK1): {} \t {}".format(counts["P53350"], counts["P53350"] == 3))
# print("P54753 (EPHB3): {} \t {}".format(counts["P54753"], counts["P54753"] == 3))
# print("P54756 (EPHA5): {} \t {}".format(counts["P54756"], counts["P54756"] == 3))
# print("P62837 (UBE2D2): {} \t {}".format(counts["P62837"], counts["P62837"] == 3))
# print("Q01196 (RUNX1): {} \t {}".format(counts["Q01196"], counts["Q01196"] == 3))
# print("Q02763 (TEK): {} \t {}".format(counts["Q02763"], counts["Q02763"] == 3))
# print("Q04721 (NOTCH2): {} \t {}".format(counts["Q04721"], counts["Q04721"] == 3))
# print("Q06141 (REG3A): {} \t {}".format(counts["Q06141"], counts["Q06141"] == 3))
# print("Q08209 (PPP3CA): {} \t {}".format(counts["Q08209"], counts["Q08209"] == 3))
# print("Q08462 (ADCY2): {} \t {}".format(counts["Q08462"], counts["Q08462"] == 3))
# print("Q13131 (PRKAA1): {} \t {}".format(counts["Q13131"], counts["Q13131"] == 3))
# print("Q13177 (PAK2): {} \t {}".format(counts["Q13177"], counts["Q13177"] == 3))
# print("Q13237 (PRKG2): {} \t {}".format(counts["Q13237"], counts["Q13237"] == 3))
# print("Q13619 (CUL4A): {} \t {}".format(counts["Q13619"], counts["Q13619"] == 3))
# print("Q14155 (ARHGEF7): {} \t {}".format(counts["Q14155"], counts["Q14155"] == 3))
# print("Q14324 (MYBPC2): {} \t {}".format(counts["Q14324"], counts["Q14324"] == 3))
# print("Q16548 (BCL2A1): {} \t {}".format(counts["Q16548"], counts["Q16548"] == 3))
# print("Q6N069 (NAA16): {} \t {}".format(counts["Q6N069"], counts["Q6N069"] == 3))
# print("Q7LDG7 (RASGRP2): {} \t {}".format(counts["Q7LDG7"], counts["Q7LDG7"] == 3))
# print("Q969H0 (FBXW7): {} \t {}".format(counts["Q969H0"], counts["Q969H0"] == 3))
# print("Q96KN2 (CNDP1): {} \t {}".format(counts["Q96KN2"], counts["Q96KN2"] == 3))
# print("Q9BZL6 (PRKD2): {} \t {}".format(counts["Q9BZL6"], counts["Q9BZL6"] == 3))
# print("Q9H1R3 (MYLK2): {} \t {}".format(counts["Q9H1R3"], counts["Q9H1R3"] == 3))
# print("Q9H832 (UBE2Z): {} \t {}".format(counts["Q9H832"], counts["Q9H832"] == 3))
# print("Q9UBT2 (UBA2): {} \t {}".format(counts["Q9UBT2"], counts["Q9UBT2"] == 3))
# print("Q9UKX5 (ITGA11): {} \t {}".format(counts["Q9UKX5"], counts["Q9UKX5"] == 3))
# print("Q9UL54 (TAOK2): {} \t {}".format(counts["Q9UL54"], counts["Q9UL54"] == 3))
# print("Q9Y566 (SHANK1): {} \t {}".format(counts["Q9Y566"], counts["Q9Y566"] == 3))
# print("Q9Y616 (IRAK3): {} \t {}".format(counts["Q9Y616"], counts["Q9Y616"] == 3))
# print("Q9Y6R4 (MAP3K4): {} \t {}".format(counts["Q9Y6R4"], counts["Q9Y6R4"] == 3))
# print("A0FGR9 (ESYT3): {} \t {}".format(counts["A0FGR9"], counts["A0FGR9"] == 2))
# print("O00203 (AP3B1): {} \t {}".format(counts["O00203"], counts["O00203"] == 2))
# print("O00238 (BMPR1B): {} \t {}".format(counts["O00238"], counts["O00238"] == 2))
# print("O14775 (GNB5): {} \t {}".format(counts["O14775"], counts["O14775"] == 2))
# print("O14827 (RASGRF2): {} \t {}".format(counts["O14827"], counts["O14827"] == 2))
# print("O14933 (UBE2L6): {} \t {}".format(counts["O14933"], counts["O14933"] == 2))
# print("O15018 (PDZD2): {} \t {}".format(counts["O15018"], counts["O15018"] == 2))
# print("O15488 (GYG2): {} \t {}".format(counts["O15488"], counts["O15488"] == 2))
# print("O43283 (MAP3K13): {} \t {}".format(counts["O43283"], counts["O43283"] == 2))
# print("O43424 (GRID2): {} \t {}".format(counts["O43424"], counts["O43424"] == 2))
# print("O43684 (BUB3): {} \t {}".format(counts["O43684"], counts["O43684"] == 2))
# print("O43791 (SPOP): {} \t {}".format(counts["O43791"], counts["O43791"] == 2))
# print("O60346 (PHLPP1): {} \t {}".format(counts["O60346"], counts["O60346"] == 2))
# print("O60383 (GDF9): {} \t {}".format(counts["O60383"], counts["O60383"] == 2))
# print("O94973 (AP2A2): {} \t {}".format(counts["O94973"], counts["O94973"] == 2))
# print("O95259 (KCNH1): {} \t {}".format(counts["O95259"], counts["O95259"] == 2))
# print("O95760 (IL33): {} \t {}".format(counts["O95760"], counts["O95760"] == 2))
# print("O95782 (AP2A1): {} \t {}".format(counts["O95782"], counts["O95782"] == 2))
# print("P00519 (ABL1): {} \t {}".format(counts["P00519"], counts["P00519"] == 2))
# print("P00747 (PLG): {} \t {}".format(counts["P00747"], counts["P00747"] == 2))
# print("P01024 (C3): {} \t {}".format(counts["P01024"], counts["P01024"] == 2))
# print("P01031 (C5): {} \t {}".format(counts["P01031"], counts["P01031"] == 2))
# print("P01034 (CST3): {} \t {}".format(counts["P01034"], counts["P01034"] == 2))
# print("P01242 (GH2): {} \t {}".format(counts["P01242"], counts["P01242"] == 2))
# print("P01574 (IFNB1): {} \t {}".format(counts["P01574"], counts["P01574"] == 2))
# print("P02549 (SPTA1): {} \t {}".format(counts["P02549"], counts["P02549"] == 2))
# print("P03372 (ESR1): {} \t {}".format(counts["P03372"], counts["P03372"] == 2))
# print("P04049 (RAF1): {} \t {}".format(counts["P04049"], counts["P04049"] == 2))
# print("P04275 (VWF): {} \t {}".format(counts["P04275"], counts["P04275"] == 2))
# print("P04908 (H2AC4): {} \t {}".format(counts["P04908"], counts["P04908"] == 2))
# print("P05107 (ITGB2): {} \t {}".format(counts["P05107"], counts["P05107"] == 2))
# print("P05129 (PRKCG): {} \t {}".format(counts["P05129"], counts["P05129"] == 2))
# print("P05771 (PRKCB): {} \t {}".format(counts["P05771"], counts["P05771"] == 2))
# print("P06239 (LCK): {} \t {}".format(counts["P06239"], counts["P06239"] == 2))
# print("P06401 (PGR): {} \t {}".format(counts["P06401"], counts["P06401"] == 2))
# print("P06493 (CDK1): {} \t {}".format(counts["P06493"], counts["P06493"] == 2))
# print("P07332 (FES): {} \t {}".format(counts["P07332"], counts["P07332"] == 2))
# print("P08263 (GSTA1): {} \t {}".format(counts["P08263"], counts["P08263"] == 2))
# print("P08476 (INHBA): {} \t {}".format(counts["P08476"], counts["P08476"] == 2))
# print("P08514 (ITGA2B): {} \t {}".format(counts["P08514"], counts["P08514"] == 2))
# print("P08709 (F7): {} \t {}".format(counts["P08709"], counts["P08709"] == 2))
# print("P08727 (KRT19): {} \t {}".format(counts["P08727"], counts["P08727"] == 2))
# print("P09619 (PDGFRB): {} \t {}".format(counts["P09619"], counts["P09619"] == 2))
# print("P0C0S5 (H2AZ1): {} \t {}".format(counts["P0C0S5"], counts["P0C0S5"] == 2))
# print("P0C0S8 (H2AC11): {} \t {}".format(counts["P0C0S8"], counts["P0C0S8"] == 2))
# print("P10398 (ARAF): {} \t {}".format(counts["P10398"], counts["P10398"] == 2))
# print("P11308 (ERG): {} \t {}".format(counts["P11308"], counts["P11308"] == 2))
# print("P12259 (F5): {} \t {}".format(counts["P12259"], counts["P12259"] == 2))
# print("P12645 (BMP3): {} \t {}".format(counts["P12645"], counts["P12645"] == 2))
# print("P12757 (SKIL): {} \t {}".format(counts["P12757"], counts["P12757"] == 2))
# print("P13349 (MYF5): {} \t {}".format(counts["P13349"], counts["P13349"] == 2))
# print("P13637 (ATP1A3): {} \t {}".format(counts["P13637"], counts["P13637"] == 2))
# print("P15056 (BRAF): {} \t {}".format(counts["P15056"], counts["P15056"] == 2))
# print("P15170 (GSPT1): {} \t {}".format(counts["P15170"], counts["P15170"] == 2))
# print("P15498 (VAV1): {} \t {}".format(counts["P15498"], counts["P15498"] == 2))
# print("P16234 (PDGFRA): {} \t {}".format(counts["P16234"], counts["P16234"] == 2))
# print("P16885 (PLCG2): {} \t {}".format(counts["P16885"], counts["P16885"] == 2))
# print("P17081 (RHOQ): {} \t {}".format(counts["P17081"], counts["P17081"] == 2))
# print("P17787 (CHRNB2): {} \t {}".format(counts["P17787"], counts["P17787"] == 2))
# print("P18433 (PTPRA): {} \t {}".format(counts["P18433"], counts["P18433"] == 2))
# print("P19105 (MYL12A): {} \t {}".format(counts["P19105"], counts["P19105"] == 2))
# print("P19429 (TNNI3): {} \t {}".format(counts["P19429"], counts["P19429"] == 2))
# print("P19838 (NFKB1): {} \t {}".format(counts["P19838"], counts["P19838"] == 2))
# print("P22314 (UBA1): {} \t {}".format(counts["P22314"], counts["P22314"] == 2))
# print("P22626 (HNRNPA2B1): {} \t {}".format(counts["P22626"], counts["P22626"] == 2))
# print("P25391 (LAMA1): {} \t {}".format(counts["P25391"], counts["P25391"] == 2))
# print("P25786 (PSMA1): {} \t {}".format(counts["P25786"], counts["P25786"] == 2))
# print("P26038 (MSN): {} \t {}".format(counts["P26038"], counts["P26038"] == 2))
# print("P27986 (PIK3R1): {} \t {}".format(counts["P27986"], counts["P27986"] == 2))
# print("P29274 (ADORA2A): {} \t {}".format(counts["P29274"], counts["P29274"] == 2))
# print("P29323 (EPHB2): {} \t {}".format(counts["P29323"], counts["P29323"] == 2))
# print("P31749 (AKT1): {} \t {}".format(counts["P31749"], counts["P31749"] == 2))
# print("P32246 (CCR1): {} \t {}".format(counts["P32246"], counts["P32246"] == 2))
# print("P34932 (HSPA4): {} \t {}".format(counts["P34932"], counts["P34932"] == 2))
# print("P35226 (BMI1): {} \t {}".format(counts["P35226"], counts["P35226"] == 2))
# print("P35251 (RFC1): {} \t {}".format(counts["P35251"], counts["P35251"] == 2))
# print("P35626 (GRK3): {} \t {}".format(counts["P35626"], counts["P35626"] == 2))
# print("P36871 (PGM1): {} \t {}".format(counts["P36871"], counts["P36871"] == 2))
# print("P36896 (ACVR1B): {} \t {}".format(counts["P36896"], counts["P36896"] == 2))
# print("P37231 (PPARG): {} \t {}".format(counts["P37231"], counts["P37231"] == 2))
# print("P38919 (EIF4A3): {} \t {}".format(counts["P38919"], counts["P38919"] == 2))
# print("P40145 (ADCY8): {} \t {}".format(counts["P40145"], counts["P40145"] == 2))
# print("P40818 (USP8): {} \t {}".format(counts["P40818"], counts["P40818"] == 2))
# print("P41235 (HNF4A): {} \t {}".format(counts["P41235"], counts["P41235"] == 2))
# print("P41743 (PRKCI): {} \t {}".format(counts["P41743"], counts["P41743"] == 2))
# print("P42336 (PIK3CA): {} \t {}".format(counts["P42336"], counts["P42336"] == 2))
# print("P42345 (MTOR): {} \t {}".format(counts["P42345"], counts["P42345"] == 2))
# print("P42681 (TXK): {} \t {}".format(counts["P42681"], counts["P42681"] == 2))
# print("P43146 (DCC): {} \t {}".format(counts["P43146"], counts["P43146"] == 2))
# print("P43405 (SYK): {} \t {}".format(counts["P43405"], counts["P43405"] == 2))
# print("P43686 (PSMC4): {} \t {}".format(counts["P43686"], counts["P43686"] == 2))
# print("P45983 (MAPK8): {} \t {}".format(counts["P45983"], counts["P45983"] == 2))
# print("P46531 (NOTCH1): {} \t {}".format(counts["P46531"], counts["P46531"] == 2))
# print("P47755 (CAPZA2): {} \t {}".format(counts["P47755"], counts["P47755"] == 2))
# print("P48051 (KCNJ6): {} \t {}".format(counts["P48051"], counts["P48051"] == 2))
# print("P49760 (CLK2): {} \t {}".format(counts["P49760"], counts["P49760"] == 2))
# print("P49802 (RGS7): {} \t {}".format(counts["P49802"], counts["P49802"] == 2))
# print("P50613 (CDK7): {} \t {}".format(counts["P50613"], counts["P50613"] == 2))
# print("P51168 (SCNN1B): {} \t {}".format(counts["P51168"], counts["P51168"] == 2))
# print("P51795 (CLCN5): {} \t {}".format(counts["P51795"], counts["P51795"] == 2))
# print("P52333 (JAK3): {} \t {}".format(counts["P52333"], counts["P52333"] == 2))
# print("P53708 (ITGA8): {} \t {}".format(counts["P53708"], counts["P53708"] == 2))
# print("P54652 (HSPA2): {} \t {}".format(counts["P54652"], counts["P54652"] == 2))
# print("P54762 (EPHB1): {} \t {}".format(counts["P54762"], counts["P54762"] == 2))
# print("P55795 (HNRNPH2): {} \t {}".format(counts["P55795"], counts["P55795"] == 2))
# print("P56199 (ITGA1): {} \t {}".format(counts["P56199"], counts["P56199"] == 2))
# print("P60709 (ACTB): {} \t {}".format(counts["P60709"], counts["P60709"] == 2))
# print("P61077 (UBE2D3): {} \t {}".format(counts["P61077"], counts["P61077"] == 2))
# print("P61106 (RAB14): {} \t {}".format(counts["P61106"], counts["P61106"] == 2))
# print("P61587 (RND3): {} \t {}".format(counts["P61587"], counts["P61587"] == 2))
# print("P62158 (-): {} \t {}".format(counts["P62158"], counts["P62158"] == 2))
# print("P62714 (PPP2CB): {} \t {}".format(counts["P62714"], counts["P62714"] == 2))
# print("P62834 (RAP1A): {} \t {}".format(counts["P62834"], counts["P62834"] == 2))
# print("P62993 (GRB2): {} \t {}".format(counts["P62993"], counts["P62993"] == 2))
# print("P81133 (SIM1): {} \t {}".format(counts["P81133"], counts["P81133"] == 2))
# print("P84022 (SMAD3): {} \t {}".format(counts["P84022"], counts["P84022"] == 2))
# print("Q01813 (PFKP): {} \t {}".format(counts["Q01813"], counts["Q01813"] == 2))
# print("Q02153 (GUCY1B1): {} \t {}".format(counts["Q02153"], counts["Q02153"] == 2))
# print("Q02156 (PRKCE): {} \t {}".format(counts["Q02156"], counts["Q02156"] == 2))
# print("Q02388 (COL7A1): {} \t {}".format(counts["Q02388"], counts["Q02388"] == 2))
# print("Q02779 (MAP3K10): {} \t {}".format(counts["Q02779"], counts["Q02779"] == 2))
# print("Q04724 (TLE1): {} \t {}".format(counts["Q04724"], counts["Q04724"] == 2))
# print("Q06187 (BTK): {} \t {}".format(counts["Q06187"], counts["Q06187"] == 2))
# print("Q07812 (BAX): {} \t {}".format(counts["Q07812"], counts["Q07812"] == 2))
# print("Q07817 (BCL2L1): {} \t {}".format(counts["Q07817"], counts["Q07817"] == 2))
# print("Q07889 (SOS1): {} \t {}".format(counts["Q07889"], counts["Q07889"] == 2))
# print("Q12905 (ILF2): {} \t {}".format(counts["Q12905"], counts["Q12905"] == 2))
# print("Q12965 (MYO1E): {} \t {}".format(counts["Q12965"], counts["Q12965"] == 2))
# print("Q13002 (GRIK2): {} \t {}".format(counts["Q13002"], counts["Q13002"] == 2))
# print("Q13087 (PDIA2): {} \t {}".format(counts["Q13087"], counts["Q13087"] == 2))
# print("Q13233 (MAP3K1): {} \t {}".format(counts["Q13233"], counts["Q13233"] == 2))
# print("Q13347 (EIF3I): {} \t {}".format(counts["Q13347"], counts["Q13347"] == 2))
# print("Q13418 (ILK): {} \t {}".format(counts["Q13418"], counts["Q13418"] == 2))
# print("Q13490 (BIRC2): {} \t {}".format(counts["Q13490"], counts["Q13490"] == 2))
# print("Q13546 (RIPK1): {} \t {}".format(counts["Q13546"], counts["Q13546"] == 2))
# print("Q13547 (HDAC1): {} \t {}".format(counts["Q13547"], counts["Q13547"] == 2))
# print("Q13952 (NFYC): {} \t {}".format(counts["Q13952"], counts["Q13952"] == 2))
# print("Q14005 (IL16): {} \t {}".format(counts["Q14005"], counts["Q14005"] == 2))
# print("Q14161 (GIT2): {} \t {}".format(counts["Q14161"], counts["Q14161"] == 2))
# print("Q14164 (IKBKE): {} \t {}".format(counts["Q14164"], counts["Q14164"] == 2))
# print("Q14240 (EIF4A2): {} \t {}".format(counts["Q14240"], counts["Q14240"] == 2))
# print("Q14318 (FKBP8): {} \t {}".format(counts["Q14318"], counts["Q14318"] == 2))
# print("Q14678 (KANK1): {} \t {}".format(counts["Q14678"], counts["Q14678"] == 2))
# print("Q14896 (MYBPC3): {} \t {}".format(counts["Q14896"], counts["Q14896"] == 2))
# print("Q15366 (PCBP2): {} \t {}".format(counts["Q15366"], counts["Q15366"] == 2))
# print("Q15375 (EPHA7): {} \t {}".format(counts["Q15375"], counts["Q15375"] == 2))
# print("Q15386 (UBE3C): {} \t {}".format(counts["Q15386"], counts["Q15386"] == 2))
# print("Q15413 (RYR3): {} \t {}".format(counts["Q15413"], counts["Q15413"] == 2))
# print("Q15700 (DLG2): {} \t {}".format(counts["Q15700"], counts["Q15700"] == 2))
# print("Q15746 (MYLK): {} \t {}".format(counts["Q15746"], counts["Q15746"] == 2))
# print("Q15842 (KCNJ8): {} \t {}".format(counts["Q15842"], counts["Q15842"] == 2))
# print("Q16288 (NTRK3): {} \t {}".format(counts["Q16288"], counts["Q16288"] == 2))
# print("Q16665 (HIF1A): {} \t {}".format(counts["Q16665"], counts["Q16665"] == 2))
# print("Q16695 (H3)-4 \t {}: {}".format(counts["Q16695"], counts["Q16695"] == 2))
# print("Q53EL6 (PDCD4): {} \t {}".format(counts["Q53EL6"], counts["Q53EL6"] == 2))
# print("Q5HYK7 (SH3D19): {} \t {}".format(counts["Q5HYK7"], counts["Q5HYK7"] == 2))
# print("Q5JSH3 (WDR44): {} \t {}".format(counts["Q5JSH3"], counts["Q5JSH3"] == 2))
# print("Q5TCZ1 (SH3PXD2A): {} \t {}".format(counts["Q5TCZ1"], counts["Q5TCZ1"] == 2))
# print("Q6PIL6 (KCNIP4): {} \t {}".format(counts["Q6PIL6"], counts["Q6PIL6"] == 2))
# print("Q6ZU15 (SEPTIN14): {} \t {}".format(counts["Q6ZU15"], counts["Q6ZU15"] == 2))
# print("Q6ZVD8 (PHLPP2): {} \t {}".format(counts["Q6ZVD8"], counts["Q6ZVD8"] == 2))
# print("Q70CQ2 (USP34): {} \t {}".format(counts["Q70CQ2"], counts["Q70CQ2"] == 2))
# print("Q7L576 (CYFIP1): {} \t {}".format(counts["Q7L576"], counts["Q7L576"] == 2))
# print("Q7L7X3 (TAOK1): {} \t {}".format(counts["Q7L7X3"], counts["Q7L7X3"] == 2))
# print("Q7L804 (RAB11FIP2): {} \t {}".format(counts["Q7L804"], counts["Q7L804"] == 2))
# print("Q86UL8 (MAGI2): {} \t {}".format(counts["Q86UL8"], counts["Q86UL8"] == 2))
# print("Q86UR5 (RIMS1): {} \t {}".format(counts["Q86UR5"], counts["Q86UR5"] == 2))
# print("Q86YT6 (MIB1): {} \t {}".format(counts["Q86YT6"], counts["Q86YT6"] == 2))
# print("Q8IYD1 (GSPT2): {} \t {}".format(counts["Q8IYD1"], counts["Q8IYD1"] == 2))
# print("Q8N680 (ZBTB2): {} \t {}".format(counts["Q8N680"], counts["Q8N680"] == 2))
# print("Q8NFW9 (MYRIP): {} \t {}".format(counts["Q8NFW9"], counts["Q8NFW9"] == 2))
# print("Q8NFX7 (STXBP6): {} \t {}".format(counts["Q8NFX7"], counts["Q8NFX7"] == 2))
# print("Q8NI35 (PATJ): {} \t {}".format(counts["Q8NI35"], counts["Q8NI35"] == 2))
# print("Q8TBX8 (PIP4K2C): {} \t {}".format(counts["Q8TBX8"], counts["Q8TBX8"] == 2))
# print("Q8WXK3 (ASB13): {} \t {}".format(counts["Q8WXK3"], counts["Q8WXK3"] == 2))
# print("Q92526 (CCT6B): {} \t {}".format(counts["Q92526"], counts["Q92526"] == 2))
# print("Q92793 (CREBBP): {} \t {}".format(counts["Q92793"], counts["Q92793"] == 2))
# print("Q92841 (DDX17): {} \t {}".format(counts["Q92841"], counts["Q92841"] == 2))
# print("Q96DI7 (SNRNP40): {} \t {}".format(counts["Q96DI7"], counts["Q96DI7"] == 2))
# print("Q96GD4 (AURKB): {} \t {}".format(counts["Q96GD4"], counts["Q96GD4"] == 2))
# print("Q96J92 (WNK4): {} \t {}".format(counts["Q96J92"], counts["Q96J92"] == 2))
# print("Q96JB8 (MPP4): {} \t {}".format(counts["Q96JB8"], counts["Q96JB8"] == 2))
# print("Q96KB5 (PBK): {} \t {}".format(counts["Q96KB5"], counts["Q96KB5"] == 2))
# print("Q96RU8 (TRIB1): {} \t {}".format(counts["Q96RU8"], counts["Q96RU8"] == 2))
# print("Q99719 (SEPTIN5): {} \t {}".format(counts["Q99719"], counts["Q99719"] == 2))
# print("Q99743 (NPAS2): {} \t {}".format(counts["Q99743"], counts["Q99743"] == 2))
# print("Q99879 (H2BC14): {} \t {}".format(counts["Q99879"], counts["Q99879"] == 2))
# print("Q99962 (SH3GL2): {} \t {}".format(counts["Q99962"], counts["Q99962"] == 2))
# print("Q9BPU6 (DPYSL5): {} \t {}".format(counts["Q9BPU6"], counts["Q9BPU6"] == 2))
# print("Q9BY11 (PACSIN1): {} \t {}".format(counts["Q9BY11"], counts["Q9BY11"] == 2))
# print("Q9BYG5 (PARD6B): {} \t {}".format(counts["Q9BYG5"], counts["Q9BYG5"] == 2))
# print("Q9H7D7 (WDR26): {} \t {}".format(counts["Q9H7D7"], counts["Q9H7D7"] == 2))
# print("Q9NPH3 (IL1RAP): {} \t {}".format(counts["Q9NPH3"], counts["Q9NPH3"] == 2))
# print("Q9NPI9 (KCNJ16): {} \t {}".format(counts["Q9NPI9"], counts["Q9NPI9"] == 2))
# print("Q9NPJ1 (MKKS): {} \t {}".format(counts["Q9NPJ1"], counts["Q9NPJ1"] == 2))
# print("Q9NRM7 (LATS2): {} \t {}".format(counts["Q9NRM7"], counts["Q9NRM7"] == 2))
# print("Q9P0L2 (MARK1): {} \t {}".format(counts["Q9P0L2"], counts["Q9P0L2"] == 2))
# print("Q9P1U1 (ACTR3B): {} \t {}".format(counts["Q9P1U1"], counts["Q9P1U1"] == 2))
# print("Q9P212 (PLCE1): {} \t {}".format(counts["Q9P212"], counts["Q9P212"] == 2))
# print("Q9P2N7 (KLHL13): {} \t {}".format(counts["Q9P2N7"], counts["Q9P2N7"] == 2))
# print("Q9UBF8 (PI4KB): {} \t {}".format(counts["Q9UBF8"], counts["Q9UBF8"] == 2))
# print("Q9UBN7 (HDAC6): {} \t {}".format(counts["Q9UBN7"], counts["Q9UBN7"] == 2))
# print("Q9UDY2 (TJP2): {} \t {}".format(counts["Q9UDY2"], counts["Q9UDY2"] == 2))
# print("Q9UI47 (CTNNA3): {} \t {}".format(counts["Q9UI47"], counts["Q9UI47"] == 2))
# print("Q9UKR3 (KLK13): {} \t {}".format(counts["Q9UKR3"], counts["Q9UKR3"] == 2))
# print("Q9UKS6 (PACSIN3): {} \t {}".format(counts["Q9UKS6"], counts["Q9UKS6"] == 2))
# print("Q9UL25 (RAB21): {} \t {}".format(counts["Q9UL25"], counts["Q9UL25"] == 2))
# print("Q9UPN9 (TRIM33): {} \t {}".format(counts["Q9UPN9"], counts["Q9UPN9"] == 2))
# print("Q9UQB8 (BAIAP2): {} \t {}".format(counts["Q9UQB8"], counts["Q9UQB8"] == 2))
# print("Q9UQQ2 (SH2B3): {} \t {}".format(counts["Q9UQQ2"], counts["Q9UQQ2"] == 2))
# print("Q9Y297 (BTRC): {} \t {}".format(counts["Q9Y297"], counts["Q9Y297"] == 2))
# print("Q9Y2A7 (NCKAP1): {} \t {}".format(counts["Q9Y2A7"], counts["Q9Y2A7"] == 2))
# print("Q9Y2U5 (MAP3K2): {} \t {}".format(counts["Q9Y2U5"], counts["Q9Y2U5"] == 2))
# print("Q9Y463 (DYRK1B): {} \t {}".format(counts["Q9Y463"], counts["Q9Y463"] == 2))
# print("Q9Y5A6 (ZSCAN21): {} \t {}".format(counts["Q9Y5A6"], counts["Q9Y5A6"] == 2))
# print("Q9Y5K6 (CD2AP): {} \t {}".format(counts["Q9Y5K6"], counts["Q9Y5K6"] == 2))
# print("Q9Y5X4 (NR2E3): {} \t {}".format(counts["Q9Y5X4"], counts["Q9Y5X4"] == 2))
# print("Q9Y6N7 (ROBO1): {} \t {}".format(counts["Q9Y6N7"], counts["Q9Y6N7"] == 2))
# print("O00141 (SGK1): {} \t {}".format(counts["O00141"], counts["O00141"] == 1))
# print("O00217 (NDUFS8): {} \t {}".format(counts["O00217"], counts["O00217"] == 1))
# print("O00339 (MATN2): {} \t {}".format(counts["O00339"], counts["O00339"] == 1))
# print("O00444 (PLK4): {} \t {}".format(counts["O00444"], counts["O00444"] == 1))
# print("O00482 (NR5A2): {} \t {}".format(counts["O00482"], counts["O00482"] == 1))
# print("O00487 (PSMD14): {} \t {}".format(counts["O00487"], counts["O00487"] == 1))
# print("O00505 (KPNA3): {} \t {}".format(counts["O00505"], counts["O00505"] == 1))
# print("O00506 (STK25): {} \t {}".format(counts["O00506"], counts["O00506"] == 1))
# print("O00629 (KPNA4): {} \t {}".format(counts["O00629"], counts["O00629"] == 1))
# print("O00635 (TRIM38): {} \t {}".format(counts["O00635"], counts["O00635"] == 1))
# print("O14578 (CIT): {} \t {}".format(counts["O14578"], counts["O14578"] == 1))
# print("O14641 (DVL2): {} \t {}".format(counts["O14641"], counts["O14641"] == 1))
# print("O14662 (STX16): {} \t {}".format(counts["O14662"], counts["O14662"] == 1))
# print("O14682 (ENC1): {} \t {}".format(counts["O14682"], counts["O14682"] == 1))
# print("O14733 (MAP2K7): {} \t {}".format(counts["O14733"], counts["O14733"] == 1))
# print("O14745 (SLC9A3R1): {} \t {}".format(counts["O14745"], counts["O14745"] == 1))
# print("O14757 (CHEK1): {} \t {}".format(counts["O14757"], counts["O14757"] == 1))
# print("O14786 (NRP1): {} \t {}".format(counts["O14786"], counts["O14786"] == 1))
# print("O14818 (PSMA7): {} \t {}".format(counts["O14818"], counts["O14818"] == 1))
# print("O14830 (PPEF2): {} \t {}".format(counts["O14830"], counts["O14830"] == 1))
# print("O14893 (GEMIN2): {} \t {}".format(counts["O14893"], counts["O14893"] == 1))
# print("O14910 (LIN7A): {} \t {}".format(counts["O14910"], counts["O14910"] == 1))
# print("O14920 (IKBKB): {} \t {}".format(counts["O14920"], counts["O14920"] == 1))
# print("O14950 (MYL12B): {} \t {}".format(counts["O14950"], counts["O14950"] == 1))
# print("O14976 (GAK): {} \t {}".format(counts["O14976"], counts["O14976"] == 1))
# print("O15041 (SEMA3E): {} \t {}".format(counts["O15041"], counts["O15041"] == 1))
# print("O15055 (PER2): {} \t {}".format(counts["O15055"], counts["O15055"] == 1))
# print("O15116 (LSM1): {} \t {}".format(counts["O15116"], counts["O15116"] == 1))
# print("O15143 (ARPC1B): {} \t {}".format(counts["O15143"], counts["O15143"] == 1))
# print("O15146 (MUSK): {} \t {}".format(counts["O15146"], counts["O15146"] == 1))
# print("O15444 (CCL25): {} \t {}".format(counts["O15444"], counts["O15444"] == 1))
# print("O15492 (RGS16): {} \t {}".format(counts["O15492"], counts["O15492"] == 1))
