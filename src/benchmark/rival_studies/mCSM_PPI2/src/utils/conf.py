import pathlib

test_mode = False

ALLOWED_RAM_PERCENTAGE: int = 93

import os.path as op

# if test_mode:
#     RECORDS_FOLDER_PATH = r"test_files\Records_Test"
#     ELASPIC_RESULTS_FOLDER_PATH = r"test_files\ELASPIC_Results_TEST"
#     UPLOAD_FAILED_PATH = r"test_files\upload_failed_records_TEST"
#     UNEXPECTED_FAILED_PATH = r"test_files\unexpected_failed_records_TEST"
#     INPUT_FILES_PATH = r"test_files\input_files_test"
#     ALLMUTATIONS_FAILED_PATH = r"test_files\allmutations_failed_TEST"
#
# # ACTUAL RUN
# else:
#     RECORDS_FOLDER_PATH = "Record"
#     ELASPIC_RESULTS_FOLDER_PATH = "Elaspic_Results"
#     UPLOAD_FAILED_PATH = "Upload_fails"
#     UNEXPECTED_FAILED_PATH = "Unexpected_fails"
#     INPUT_FILES_PATH = "ELASPIC_Input"
#     ALLMUTATIONS_FAILED_PATH = r"Allmutations_fails"

COMPUTATION_TIME_ALLOWED = 1  # 3 # 10  # in seconds.

# Paths
CHUNKS_TO_RUN_FOLDER_PATH = "Chunks_to_run"
DRIVER_PATH = r"C:\webdrivers\geckodriver.exe"
TEMP_DOWNLOAD_FOLDER_PATH = op.join("..", pathlib.Path().resolve(), "Firefox_download")
MCSM_PPI_MANY_URL = "http://biosig.unimelb.edu.au/mcsm_ppi2/submit_prediction"
PDB_FILES_PATH = op.join(pathlib.Path().resolve(), "..", "data", "pdb_files")
RECORD_PATH = "Record"

HEADLESS = False

# r"C:\Users\ibrah\Desktop\Spaceship\ELASPIC_cancer_data_smaller_chunks\ELASPIC_Input\BRCA_10\SNV_BRCA_Chunk_22_0_test.txt"
# UPLOAD_FILE_PATH = r"C:\Users\ibrah\Documents\GitHub\My-ELASPIC-Web-API\input_files_test\mutations_input_test.txt"

# Firefox Log output disabled
PATH_TO_DEV_NULL = 'nul'
