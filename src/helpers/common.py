import os


def print_annotation(s):
    print(f"\n{s}\n{'-' * len(s)}")


def get_file_path(project_common_file_dir, filename):
    return os.path.join(project_common_file_dir, filename)
