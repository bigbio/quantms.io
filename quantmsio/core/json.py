import json
import logging
import re
import numpy as np
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.core.tools import load_de_or_ae
from quantmsio.core.tools import read_large_parquet


class Npencoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.int32):
            return int(obj)
        return json.JSONEncoder.default(self, obj)


def scores_to_json(scores_row) -> list:
    """
    List of string to list of dictionary
    :param scores_row: a row in scores parquet file
    :return: a list of dictionary
    """
    result = []
    for string in scores_row:
        if string is None:
            logging.warning("The score is None")
        else:
            score_string = string.split(":")
            if len(score_string) == 2:
                result.append({score_string[0]: score_string[1]})
            elif len(score_string) > 2:
                result.append({":".join(score_string[:-1]): score_string[-1]})
            else:
                logging.warning("The score is not in the right format score: value")
    return result


def feature_json(feature_row) -> dict:
    """
    this function converts a feature row in parquet to json dictionary.
    :param feature_row: a row in feature parquet file
    :return: a json dictionary
    """
    feature_dic = feature_row.to_dict()

    feature_dic["protein_accessions"] = (
        list(feature_dic["protein_accessions"]) if feature_dic["protein_accessions"].any() else []
    )
    feature_dic["protein_start_positions"] = (
        list(feature_dic["protein_start_positions"]) if feature_dic["protein_start_positions"].any() else []
    )
    feature_dic["protein_end_positions"] = (
        list(feature_dic["protein_end_positions"]) if feature_dic["protein_end_positions"].any() else []
    )
    feature_dic["modifications"] = (
        list(feature_dic["modifications"])
        if "modifications" in feature_dic
        and feature_dic["modifications"] is not None
        and feature_dic["modifications"].any()
        else []
    )
    feature_dic["gene_names"] = (
        list(feature_dic["gene_names"])
        if feature_dic["gene_names"] is not None and feature_dic["gene_names"].any()
        else []
    )
    feature_dic["gene_accessions"] = (
        list(feature_dic["gene_accessions"])
        if (feature_dic["gene_accessions"] is not None and feature_dic["gene_accessions"].any())
        else []
    )

    feature_dic["id_scores"] = scores_to_json(list(feature_dic["id_scores"]))

    return feature_dic


class JsonConverter:

    def __init__(self):
        super().__init__()

    def feature_to_json(self, parquet_feature_path: str, json_feature_path: str):
        """
        convert the feature file to json format.
        the method uses an iterative method to convert the feature file to json in batches.
        :param parquet_feature_path: the path of feature file
        :param json_feature_path: the path of json file
        :return: None
        """
        chunks = read_large_parquet(parquet_feature_path)
        json_file = open(json_feature_path, "w")
        for feature_df in chunks:
            for index, row in feature_df.iterrows():
                json_obj = feature_json(row)
                json.dump(json_obj, json_file, cls=Npencoder)
                json_file.write("\n")
            logging.info(f"Finished converting {len(feature_df)} rows")
        json_file.close()

    def psm_to_json(self, parquet_psm_path: str, json_psm_path: str):
        """
        Convert the psm file to json format.
        :param parquet_psm_path: PSM in parquet format
        :param json_psm_path: PSMs in json format
        :return: None
        """

        chunks = read_large_parquet(parquet_psm_path)
        for psm_df in chunks:
            psm_df.to_json(json_psm_path, orient="records", lines=True, compression="gzip")
        return json_psm_path

    def convert_tsv_to_json(self, file_path: str, json_path: str):
        """
        by providing the json format of AE and DE files for retrieval. return json
        """
        table, content = load_de_or_ae(file_path)
        output = {}
        pattern = r"[\\|\|//|/]"
        file_name = re.split(pattern, file_path)[-1]
        output["id"] = file_name
        output["metadata"] = content
        records = {}
        for col in table.columns:
            records[col] = table.loc[:, col].to_list()
        output["records"] = records
        b = json.dumps(output)
        f = open(json_path, "w")
        f.write(b)
        f.close()

    def sdrf_to_json(self, sdrf_path: str, json_path: str):
        sdrf_handler = SDRFHandler(sdrf_path)
        sdrf_handler._load_sdrf_info(sdrf_path)
        sdrf_json = {}
        sdrf_json["organism"] = sdrf_handler.get_organisms()
        sdrf_json["instrument"] = sdrf_handler.get_instruments()
        sdrf_json["diseases"] = sdrf_handler.get_diseases()
        sdrf_json["enzymes"] = sdrf_handler.get_enzymes()
        sdrf_json["cell_lines"] = sdrf_handler.get_cell_lines()
        sdrf_json["method"] = sdrf_handler.get_acquisition_properties()
        sdrf_json["experiment_type"] = sdrf_handler.get_experiment_type_from_sdrf()
        sdrf_json["feature_msg"] = sdrf_handler.extract_feature_properties().to_json(orient="values")
        sdrf_json["sample_map"] = sdrf_handler.get_sample_map()
        fix_m, var_m = sdrf_handler.get_mods()

        sdrf_json["modifications"] = {"Fixed": fix_m, "Variable": var_m}
        b = json.dumps(sdrf_json)
        f = open(json_path, "w")
        f.write(b)
        f.close()
