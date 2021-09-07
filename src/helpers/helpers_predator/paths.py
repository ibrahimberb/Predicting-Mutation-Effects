from ..mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class Paths:
    def __init__(
            self,
            project_common_file_dir,
            mutations_path,
            tcga_code_path_pairs,
            initial_columns_path
    ):
        self.project_common_file_dir = project_common_file_dir
        self.mutations_path = mutations_path
        self.tcga_code_path_pairs = tcga_code_path_pairs
        self.initial_columns_path = initial_columns_path
