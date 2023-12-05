import os
import uuid
from pathlib import Path
import pandas as pd
from quantms_io.core.project import ProjectHandler
from quantms_io.core.sdrf import SDRFHandler
from quantms_io.utils.file_utils import delete_files_extension
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_ibaq_columns(path):
    with open(path) as f:
        line = f.readline()
        return line.split('\n')[0].split(',')


class AbsoluteExpressionHander:
    LABEL_MAP = {
        'ProteinName': 'protein',
        'SampleID': 'sample_accession',
        'Condition': 'condition',
        'Ibaq': 'ibaq',
        'IbaqLog': 'ribaq'
    }
    AE_HEADER = """#INFO=<ID=protein, Number=inf, Type=String, Description="Protein Accession">
#INFO=<ID=sample_accession, Number=1, Type=String, Description="Sample Accession in the SDRF">
#INFO=<ID=condition, Number=1, Type=String, Description="Value of the factor value">
#INFO=<ID=ibaq, Number=1, Type=Float, Description="Intensity based absolute quantification">
#INFO=<ID=ribaq, Number=1, Type=Float, Description="relative iBAQ">\n"""

    ABSOLUTE_EXPRESSION_EXTENSION = ".absolute.tsv"

    def __init__(self):
        self.ibaq_df = None
        self.ae_file_path = None
        self.project_manager = None
        self.sdrf_manager = None
        self.sdrf_file_path = None

    def load_project_file(self, project_file: str):
        """
        Load a project file that link the different files in the quamtms.io format
        https://github.com/bigbio/quantms.io/blob/main/docs/PROJECT.md
        :param project_file: project file path
        :return: none
        """
        self.project_file = project_file

        if not os.path.isfile(project_file):
            raise FileNotFoundError("Project file not found: " + project_file)

        self.project_manager = ProjectHandler()
        self.project_manager.load_project_info(project_file)

    def load_ibaq_file(self, path):
        usecols = ['ProteinName', 'SampleID', 'Condition', 'Ibaq', 'IbaqLog']
        ibaq_columns = get_ibaq_columns(path)
        for col in usecols:
            if col not in ibaq_columns:
                raise Exception(f"Not found {col} in ibaq file")
        ibaqs = pd.read_csv(path, usecols=usecols)
        ibaqs.rename(columns=AbsoluteExpressionHander.LABEL_MAP, inplace=True)
        self.ae_file_path = path
        self.ibaq_df = ibaqs

    def load_sdrf_file(self, sdrf_file: str):
        self.sdrf_file_path = sdrf_file

        if not os.path.isfile(sdrf_file):
            raise FileNotFoundError("SDRF file not found: " + sdrf_file)

        self.sdrf_manager = SDRFHandler(sdrf_file=sdrf_file)

    def convert_ibaq_to_quantms(
            self,
            output_folder: str = None,
            output_file_prefix: str = None,
            delete_existing: bool = False,
    ):
        output_lines = ''
        if self.project_manager:
            output_lines += (
                    "#project_accession: "
                    + self.project_manager.project.project_info["project_accession"]
                    + "\n"
            )
            output_lines += (
                    "#project_title: "
                    + self.project_manager.project.project_info["project_title"]
                    + "\n"
            )
            output_lines += (
                    "#project_description: "
                    + self.project_manager.project.project_info["project_description"]
                    + "\n"
            )
            output_lines += (
                    "#quantms_version: "
                    + self.project_manager.project.project_info["quantms_version"]
                    + "\n"
            )
        factor_value = self.get_factor_value()
        if factor_value is not None:
            output_lines += "#factor_value: " + factor_value + "\n"
        # Combine comments and DataFrame into a single list
        output_lines += AbsoluteExpressionHander.AE_HEADER + str(
            self.ibaq_df.to_csv(sep="\t", index=False, header=True)
        )
        output_lines = output_lines.replace('\r', '')
        # Create the output file name
        base_name = output_file_prefix
        if output_file_prefix is None:
            base_name = os.path.basename(self.ae_file_path).replace(".csv", "")

        # Create the output folder if it does not exist.
        if output_folder is not None and not os.path.exists(output_folder):
            Path(output_folder).mkdir(parents=True, exist_ok=True)

        # Delete existing SDRF file
        if delete_existing:
            delete_files_extension(
                output_folder,
                AbsoluteExpressionHander.ABSOLUTE_EXPRESSION_EXTENSION,
            )

        output_filename = f"{base_name}-{str(uuid.uuid4())}{AbsoluteExpressionHander.ABSOLUTE_EXPRESSION_EXTENSION}"
        if output_folder is None:
            output_filename_path = output_filename
        else:
            output_filename_path = f"{output_folder}/{output_filename}"

        # Save the combined lines to a TSV file
        with open(output_filename_path, "w", encoding='utf8') as f:
            f.write(output_lines)

        if self.project_manager:
            self.project_manager.add_quantms_file(
                file_category="absolute_file", file_name=output_filename
            )
        logger.info(
            f"Absolute expression file copied to {output_filename} and added to the project information"
        )

    def get_factor_value(self):
        """
        Get the factor value from the SDRF file
        """
        if self.sdrf_manager is None:
            return None
        return self.sdrf_manager.get_factor_value()

    def update_project_file(self, project_file: str = None):
        """
        Update the project file with the differential expression file
        :param project_file: project file path
        :return: none
        """
        if self.project_manager is None:
            raise Exception("Project file not loaded")

        if project_file is not None:
            project_file = self.project_file
        self.project_manager.save_updated_project_info(output_file_name=project_file)
