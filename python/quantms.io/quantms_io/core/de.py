import os
import uuid
from pathlib import Path

import pandas as pd

from quantms_io.core.project import ProjectHandler
from quantms_io.utils.file_utils import delete_files_extension


class DifferentialExpressionHandler:

    PROTEIN_ACCESSION_COLUMN = {"msstats_column": "protein", "quantms_column": "protein"}
    LABEL_COLUMN = {"msstats_column": "label", "quantms_column": "label"}
    log2FC_COLUMN = {"msstats_column": "log2fc", "quantms_column": "log2FC"}
    SE_COLUMN = {"msstats_column": "se", "quantms_column": "se"}
    DF_COLUMN = {"msstats_column": "df", "quantms_column": "df"}
    PVALUE_COLUMN = {"msstats_column": "pvalue", "quantms_column": "pvalue"}
    ADJUST_PVALUE_COLUMN = {"msstats_column": "adj.pvalue", "quantms_column": "adj.pvalue"}
    ISSUE_COLUMN = {"msstats_column": "issue", "quantms_column": "issue"}

    DE_HEADER = '''#INFO=<ID=protein, Number=inf, Type=String, Description="Protein Accession">
#INFO=<ID=label,Number=1, Type=String, Description="Label for the Conditions combination">
#INFO=<ID=log2fc, Number=1, Type=Double, Description="Log2 Fold Change">
#INFO=<ID=se, Number=1, Type=Double, Description="Standard error of the log2 fold change">
#INFO=<ID=df, Number=1, Type=Integer, Description="Degree of freedom of the Student test"> 
#INFO=<ID=pvalue, Number=1, Type=Double, Description="Raw p-values">
#INFO=<ID=adj.pvalue, Number=1, Type=Double, Description="P-values adjusted among all the proteins in the specific comparison using the approach by Benjamini and Hochberg">
#INFO=<ID=issue, Number=1, Type=String, Description="Issue column shows if there is any issue for inference in corresponding protein and comparison">\n'''
    
    DIFFERENTIAL_EXPRESSION_EXTENSION = ".differential.tsv"

    def __init__(self):
        """
        Differential Expression Handler class enable to create, update and retrieve information about the
        differential expression analysis performed in a project; documentation of the format can be found in
        https://github.com/bigbio/quantms.io/blob/main/docs/DE.md
        """
        self.project_file = None
        self.project_manager = None
        self.msstats_df = None
        self.de_file_path = None

    def load_msstats_file(self, msstats_file_path: str):
        """
        Load a MSstats differential file
        :param msstats_file_path: MSstats differential file path
        :return: none
        """
        self.de_file_path = msstats_file_path

        if not os.path.isfile(msstats_file_path):
            raise FileNotFoundError("MSstats differential file not found: " + msstats_file_path)

        self.msstats_df = pd.read_csv(msstats_file_path, sep="\t")
        # Rename columns to lower case
        self.msstats_df.columns = self.msstats_df.columns.str.lower()


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

    def convert_msstats_to_quantms(self, output_folder: str = None, output_file_prefix: str = None, delete_existing: bool = False):
        """
        Convert a MSstats differential file to quantms.io format
        :return: none
        """
        if self.msstats_df is None or self.project_manager is None:
            raise Exception("MSstats file or project file not loaded")

        quantms_df = self.msstats_df[[DifferentialExpressionHandler.PROTEIN_ACCESSION_COLUMN["msstats_column"],
                                        DifferentialExpressionHandler.LABEL_COLUMN["msstats_column"],
                                        DifferentialExpressionHandler.log2FC_COLUMN["msstats_column"],
                                        DifferentialExpressionHandler.SE_COLUMN["msstats_column"],
                                        DifferentialExpressionHandler.DF_COLUMN["msstats_column"],
                                        DifferentialExpressionHandler.PVALUE_COLUMN["msstats_column"],
                                        DifferentialExpressionHandler.ADJUST_PVALUE_COLUMN["msstats_column"],
                                        DifferentialExpressionHandler.ISSUE_COLUMN["msstats_column"]]].copy()

        # Combine comments and DataFrame into a single list
        output_lines = DifferentialExpressionHandler.DE_HEADER + str(quantms_df.to_csv(sep='\t', index=False, header=True))

        # Create the output file name
        base_name = output_file_prefix
        if output_file_prefix is None:
            base_name = os.path.basename(self.de_file_path).replace(".csv", "")
        
        # Create the output folder if it does not exist.
        if output_folder is not None and not os.path.exists(output_folder):
            Path(output_folder).mkdir(parents=True, exist_ok=True)

        ## Delete existing SDRF file
        if delete_existing:
            delete_files_extension(output_folder, DifferentialExpressionHandler.DIFFERENTIAL_EXPRESSION_EXTENSION)

        output_filename = f"{base_name}-{str(uuid.uuid4())}{DifferentialExpressionHandler.DIFFERENTIAL_EXPRESSION_EXTENSION}"
        if output_folder is None:
            output_filename_path = output_filename
        else:
            output_filename_path = f"{output_folder}/{base_name}-{str(uuid.uuid4())}{DifferentialExpressionHandler.DIFFERENTIAL_EXPRESSION_EXTENSION}"

        # Save the combined lines to a TSV file
        with open(output_filename_path, 'w') as f:
            f.write(output_lines)

        self.project_manager.add_quantms_file(file_category="differential_file", file_name = output_filename)
        print(f"Differential expression file copied to {output_filename} and added to the project information")

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
        self.project_manager.save_updated_project_info(output_file_name = project_file)
        
        

        










