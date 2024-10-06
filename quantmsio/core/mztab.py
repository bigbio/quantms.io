import codecs
import os
import pandas as pd
from quantmsio.utils.pride_utils import get_quantmsio_modifications


def generate_modification_list(modification_str: str, modifications):

    if pd.isna(modification_str):
        return None
    modifications = get_quantmsio_modifications(modification_str, modifications)
    modifications_string = ""
    for key, value in modifications.items():
        modifications_string += "|".join(map(str, value["position"]))
        modifications_string = modifications_string + "-" + value["unimod_accession"] + ","
    modifications_string = modifications_string[:-1]  # Remove last comma
    modification_list = modifications_string.split(",")

    return modification_list


def fetch_modifications_from_mztab_line(line: str, _modifications: dict) -> dict:
    """
    get the modifications from a mztab line. An mzTab modification could be a fixed or variable modification.
    The structure of a fixed is the following:
      MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
      MTD	fixed_mod[1]-site	C
      MTD	fixed_mod[1]-position	Anywhere
    while the structure of a variable modification is the following:
      MTD	var_mod[1]	[UNIMOD, UNIMOD:21, Phospho, ]
      MTD	var_mod[1]-site	S
      MTD   var_mod[1]-position	Anywhere

    :param line: mztab line
    :param _modifications: modifications dictionary
    :return: modification dictionary
    """
    line = line.strip()
    line_parts = line.split("\t")
    if line_parts[0] == "MTD" and "_mod[" in line_parts[1]:
        if "site" not in line_parts[1] and "position" not in line_parts[1]:
            values = line_parts[2].replace("[", "").replace("]", "").split(",")
            accession = values[1].strip()
            name = values[2].strip()
            index = line_parts[1].split("[")[1].split("]")[0]
            _modifications[accession] = [name, index, None, None]
        elif "site" in line_parts[1]:
            index = line_parts[1].split("[")[1].split("]")[0]
            accession = None
            for (
                key,
                value,
            ) in _modifications.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
                if value[1] == index:
                    accession = key
            if accession is None:
                raise Exception("The accession for the modification is None")
            _modifications[accession][2] = line_parts[2]
        elif "position" in line_parts[1]:
            index = line_parts[1].split("[")[1].split("]")[0]
            accession = None
            for key, value in _modifications.items():
                if value[1] == index:
                    accession = key
            if accession is None:
                raise Exception("The accession for the modification is None")
            _modifications[accession][3] = line_parts[2]
    return _modifications


class MzTab:
    def __init__(self, mzTab_path: str) -> None:
        self.mztab_path = mzTab_path
        # psm pos
        self._psm_pos = None
        # psm len
        self._psm_len = None
        # pep pos
        self._pep_pos = None
        # pep len
        self._pep_len = None
        # prt pos
        self._prt_pos = None
        # prt len
        self._prt_len = None
        # load psms columns
        self._psms_columns = None
        # load pep columns
        self._pep_columns = None

    def __get_pos(self, header):
        if header == "PSH" and self._pep_pos is not None:
            return self._pep_pos + self._pep_len - 1
        elif header == "PEH" and self._prt_pos is not None:
            return self._prt_pos + self._prt_len - 1
        else:
            return 0

    def __extract_len(self, header):
        map_tag = {"PSH": "PSM", "PEH": "PEP", "PRH": "PRT"}
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = open(self.mztab_path)
        pos = self.__get_pos(header)
        f.seek(pos)
        line = f.readline()
        while not line.startswith(header):
            pos = f.tell()
            line = f.readline()

        if header == "PSH":
            self._psms_columns = line.split("\n")[0].split("\t")
        if header == "PEH":
            self._pep_columns = line.split("\n")[0].split("\t")

        line = f.readline()
        fle_len = 0
        while line.startswith(map_tag[header]):
            fle_len += 1
            line = f.readline()
        f.close()
        return fle_len, pos

    def __load_second(self, header, **kwargs):
        f = open(self.mztab_path)
        if header == "PSH":
            f.seek(self._psm_pos)
            return pd.read_csv(f, nrows=self._psm_len, **kwargs)
        elif header == "PEH":
            f.seek(self._pep_pos)
            return pd.read_csv(f, nrows=self._pep_len, **kwargs)
        else:
            f.seek(self._prt_pos)
            return pd.read_csv(f, nrows=self._prt_len, **kwargs)

    def __set_table_config(self, header, length, pos):
        if header == "PSH":
            self._psm_pos = pos
            self._psm_len = length
        elif header == "PEH":
            self._pep_pos = pos
            self._pep_len = length
        else:
            self._prt_pos = pos
            self._prt_len = length

    def skip_and_load_csv(self, header, **kwargs):
        if self._psm_pos is not None and header == "PSH":
            return self.__load_second(header, **kwargs)
        if self._pep_pos is not None and header == "PEH":
            return self.__load_second(header, **kwargs)
        if self._prt_pos is not None and header == "PRH":
            return self.__load_second(header, **kwargs)
        fle_len, pos = self.__extract_len(header)
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = open(self.mztab_path)
        f.seek(pos)
        self.__set_table_config(header, fle_len, pos)
        return pd.read_csv(f, nrows=fle_len, sep="\t", **kwargs)

    def extract_ms_runs(self):
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(self.mztab_path, "r", "utf-8")
        line = f.readline()
        ms_runs = {}
        while line.split("\t")[0] == "MTD":
            if line.split("\t")[1].split("-")[-1] == "location":
                ms_runs[line.split("\t")[1].split("-")[0]] = line.split("\t")[2].split("//")[-1].split(".")[0]
            line = f.readline()
        f.close()
        return ms_runs

    def get_protein_map(self, protein_str=None):
        """
        return: a dict about protein score
        """
        prt = self.skip_and_load_csv(
            "PRH",
            usecols=["ambiguity_members", "best_search_engine_score[1]"],
        )
        if protein_str:
            prt = prt[prt["ambiguity_members"].str.contains(f"{protein_str}", na=False)]
        prt_score = prt.groupby("ambiguity_members").min()
        protein_map = prt_score.to_dict()["best_search_engine_score[1]"]
        return protein_map

    def get_score_names(self):
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(self.mztab_path, "r", "utf-8")
        line = f.readline()
        score_names = {}
        while line.split("\t")[0] == "MTD":
            if "psm_search_engine_score" in line:
                msgs = line.split("\t")
                score_values = msgs[2].replace("[", "").replace("]", "").split(",")
                score_name = score_values[2].strip()
                if ":" in score_name:
                    score_name = score_name.split(":")[0]
                score_names[score_name] = msgs[1].replace("psm_", "")
            line = f.readline()
        f.close()
        return score_names

    def generate_positions(self, start, end) -> list:
        start = start.split(",")
        end = end.split(",")
        return [start + ":" + end for start, end in zip(start, end)]

    def get_modifications(self):
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = codecs.open(self.mztab_path, "r", "utf-8")
        line = f.readline()
        mod_dict = {}
        while line.split("\t")[0] == "MTD":
            if "_mod[" in line:
                mod_dict = fetch_modifications_from_mztab_line(line, mod_dict)
            line = f.readline()
        f.close()
        return mod_dict
