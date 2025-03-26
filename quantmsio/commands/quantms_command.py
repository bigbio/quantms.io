from pathlib import Path
import os
import subprocess
from typing import Optional

import click


def is_diann(directory: str) -> bool:
    dirs = [
        diann_dir
        for diann_dir in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, diann_dir))
    ]
    if "diannsummary" in dirs:
        return True
    else:
        return False


def check_dir(folder_path: str) -> None:
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def find_file(directory: str, kind: str) -> list[Path]:
    path = Path(directory)
    ae_files = list(path.rglob(f"*{kind}"))
    return ae_files


def run_task(command: list) -> bool:
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        print(result.returncode)
        return True
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.stderr)
        return False


def quantmsio_workflow(directory: str, output_root_folder: str) -> None:
    dirs = [
        os.path.join(directory, diann_dir)
        for diann_dir in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, diann_dir))
    ]
    print(dirs)
    for diann_dir in dirs:
        if is_diann(diann_dir):
            report_file = find_file(diann_dir, "diann_report.tsv")[0]
            mzmlstatistics = os.path.join(diann_dir, "mzmlstatistics")
            sdrf_file = find_file(diann_dir, ".sdrf.tsv")[0]
            filename = os.path.basename(diann_dir)
            output_folder = os.path.join(directory, output_root_folder, filename)
            check_dir(output_folder)
            command = [
                "quantmsioc",
                "convert-diann",
                "--report_path",
                report_file,
                "--qvalue_threshold",
                "0.05",
                "--mzml_info_folder",
                mzmlstatistics,
                "--sdrf_path",
                sdrf_file,
                "--output_folder",
                output_folder,
                "--duckdb_max_memory",
                "64GB",
                "--output_prefix_file",
                filename,
            ]
            run_task(command)
        else:
            mztab_file = find_file(diann_dir, ".mzTab")
            if len(mztab_file) > 0:
                mztab_file = mztab_file[0]
                msstatsin_file = find_file(diann_dir, "msstats_in.csv")[0]
                sdrf_file = find_file(diann_dir, ".sdrf.tsv")[0]
                filename = os.path.basename(diann_dir)
                output_folder = os.path.join(directory, output_root_folder, filename)
                check_dir(output_folder)
                command_feature = [
                    "quantmsioc",
                    "convert-feature",
                    "--sdrf_file",
                    sdrf_file,
                    "--msstats_file",
                    msstatsin_file,
                    "--mztab_file",
                    mztab_file,
                    "--file_num",
                    "30",
                    "--output_folder",
                    output_folder,
                    "--duckdb_max_memory",
                    "64GB",
                    "--output_prefix_file",
                    filename,
                ]
                run_task(command_feature)
                command_psm = [
                    "quantmsioc",
                    "convert-psm",
                    "--mztab_file",
                    mztab_file,
                    "--output_folder",
                    output_folder,
                    "--output_prefix_file",
                    filename,
                ]
                run_task(command_psm)
            else:
                continue


@click.command(
    "convert-quantms—project",
    short_help="Generate quantmsio files from multiple quantms reports." "format",
)
@click.option(
    "--base_folder",
    help="The directory for storing project files",
    required=True,
)
@click.option(
    "--save_folder",
    help="The directory for saving the results",
    required=False,
)
def convert_quantms_project(
    base_folder: str,
    save_folder: Optional[str],
) -> None:
    if not save_folder:
        save_folder = "quantms_data"
    quantmsio_workflow(base_folder, save_folder)