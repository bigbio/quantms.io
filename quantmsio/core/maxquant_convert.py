import pandas as pd
import pyarrow.parquet as pq
import re
import pyarrow as pa
import zipfile
from pathlib import Path

#from quantmsio.core.diann_convert import find_modification
from quantmsio.utils.pride_utils import get_quantmsio_modifications
from quantmsio.utils.pride_utils import get_peptidoform_proforma_version_in_mztab

import logging
#format the log entries
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

MODIFICATION_PATTERN = re.compile(r"\((.*?\))\)")
FEATURE_FIELDS = [
    pa.field(
        "sequence",
        pa.string(),
        metadata={"description": "Peptide sequence of the feature"},
    ),
    pa.field(
        "protein_accessions",
        pa.list_(pa.string()),
        metadata={"description": "accessions of associated proteins"},
    ),
    pa.field(
        "unique",
        pa.int32(),
        metadata={"description": "if the peptide is unique to a particular protein"},
    ),
    pa.field(
        "modifications",
        pa.list_(pa.string()),
        metadata={"description": "peptide modifications"},
    ),
    pa.field(
        "charge",
        pa.int32(),
        metadata={"description": "charge state of the feature"},
    ),
    pa.field(
        "peptidoform",
        pa.string(),
        metadata={"description": "peptidoform in proforma notation"},
    ),
    pa.field(
        "posterior_error_probability",
        pa.float64(),
        metadata={"description": "posterior error probability"},
    ),
    pa.field("intensity", pa.float64(), metadata={"description": "intensity value"}),
    pa.field(
        "sample_accession",
        pa.string(),
        metadata={"description": "accession of the associated sample"},
    ),
    pa.field(
        "condition",
        pa.string(),
        metadata={"description": "experimental condition, value of the experimental factor"},
    ),
    pa.field("fraction", pa.string(), metadata={"description": "fraction information"}),
    pa.field(
        "biological_replicate",
        pa.string(),
        metadata={"description": "biological replicate information"},
    ),
    pa.field(
        "fragment_ion",
        pa.string(),
        metadata={"description": "fragment ion information"},
    ),
    pa.field(
        "isotope_label_type",
        pa.string(),
        metadata={"description": "type of isotope label"},
    ),
    pa.field("run", pa.string(), metadata={"description": "experimental run information"}),
    pa.field(
        "channel",
        pa.string(),
        metadata={"description": "experimental channel information"},
    ),
    pa.field(
        "reference_file_name",
        pa.string(),
        metadata={"description": "file name of the reference file for the feature"},
    ),
]
SCHEMA =  pa.schema(
            FEATURE_FIELDS,
            metadata={"description": "Feature file in quantms.io format"},
        )

def find_modification(peptide):
    """
    Identify the modification site based on the peptide containing modifications.

    :param peptide: Sequences of peptides
    :type peptide: str
    :return: Modification sites
    :rtype: str

    Examples:
    >>> find_modification("PEPM(UNIMOD:35)IDE")
    '4-UNIMOD:35'
    >>> find_modification("SM(UNIMOD:35)EWEIRDS(UNIMOD:21)EPTIDEK")
    '2-UNIMOD:35,9-UNIMOD:21'
    """
    peptide = str(peptide)
    original_mods = MODIFICATION_PATTERN.findall(peptide)
    peptide = MODIFICATION_PATTERN.sub(".", peptide)
    position = [i for i, x in enumerate(peptide) if x == "."]
    for j in range(1, len(position)):
        position[j] -= j

    for k in range(0, len(original_mods)):
        original_mods[k] = str(position[k]) + "-" + original_mods[k].upper()

    original_mods = ",".join(str(i) for i in original_mods) if len(original_mods) > 0 else "null"

    return original_mods

def get_mods(sdrf_path):
    sdrf = pd.read_csv(sdrf_path, sep="\t", nrows=1)
    mod_cols = [col for col in sdrf.columns if col.startswith("comment[modification parameter]")]
    if not mod_cols:
        mod_cols = [col for col in sdrf.columns if col.startswith("comment[modification parameters]")]
    mods = {}
    for i,col in enumerate(mod_cols):
        mod_msg = sdrf[col].values[0].split(";")
        mod_dict = {k.split("=")[0]: k.split("=")[1] for k in mod_msg}
        mod = [mod_dict['NT'],str(i),mod_dict['TA'] if "TA" in mod_dict else 'X',mod_dict['PP'] if 'PP' in mod_dict else 'Anywhere']
        mods[mod_dict['AC']] = mod

    return mods

def get_mod_map(sdrf_path):
    sdrf = pd.read_csv(sdrf_path, sep="\t", nrows=1)
    mod_cols = [col for col in sdrf.columns if col.startswith("comment[modification parameters]")]
    mod_map = {}
    for col in mod_cols:
        mod_msg = sdrf[col].values[0].split(";")
        mod_dict = {k.split("=")[0]: k.split("=")[1] for k in mod_msg}
        mod = f"{mod_dict['NT']} ({mod_dict['TA']})" if "TA" in mod_dict else f"{mod_dict['NT']} ({mod_dict['PP']})"
        mod_map[mod] = mod_dict['AC']

    return mod_map

def generate_mods(row, mod_map):
    mod_seq = row['Modified sequence'].replace("_","")
    mod_p = find_modification(mod_seq)
    if mod_p== 'null' or mod_p==None:
        return None
    for mod in row['Modifications'].split(','):
        mod = re.search(r"[A-Za-z]+.*\)$", mod)
        if mod:
            mod = mod.group()
            if mod in mod_map.keys():
                if '(' in mod_p:
                    mod_p = mod_p.replace(mod.upper(), mod_map[mod])
                else:
                    mod_p = mod_p.replace(mod[:2].upper(), mod_map[mod])
    return mod_p

class MaxquantConvert:
    def __init__(self):
        self.map_column_names = {
            'Sequence': 'sequence',
            'Proteins': 'protein_accessions',
            'PEP': 'posterior_error_probability',
            'Modifications': 'modifications',
            'Charge': 'charge',
            'Modified sequence': 'peptidoform',
            'Raw file': 'reference_file_name',
            'Intensity': "intensity",
        }
        self.use_columns = list(self.map_column_names.keys())
        self.sdrf_map = {
            "source name": "sample_accession",
            "comment[fraction identifier]": "fraction",
            "characteristics[biological replicate]": "biological_replicate",
            "comment[label]": "channel",
        }
        self._modifications = None

    def _generate_modification_list(self, modification_str: str):
        if pd.isna(modification_str):
            return pd.NA
        modifications = get_quantmsio_modifications(modification_str, self._modifications)
        modifications_string = ""
        for key, value in modifications.items():
            modifications_string += "|".join(map(str, value["position"]))
            modifications_string = modifications_string + "-" + value["unimod_accession"] + ","
        modifications_string = modifications_string[:-1]  # Remove last comma
        modification_list = modifications_string.split(",")

        return modification_list
    
    def main_operate(self,df:pd.DataFrame,mods_map:dict):
        df.loc[:,'Modifications'] = df[['Modified sequence','Modifications']].apply(lambda row: generate_mods(row, mods_map),axis=1)
        df = df.query('`Potential contaminant`!="+"')
        df = df.dropna(subset=['Intensity','Proteins'])
        df = df.drop('Potential contaminant', axis=1)
        df = df[df['PEP']<0.01]
        df = df.rename(columns=self.map_column_names)
        df.loc[:, "peptidoform"] = df[["sequence", "modifications"]].apply(
            lambda row: get_peptidoform_proforma_version_in_mztab(
                row["sequence"], row["modifications"], self._modifications
            ),
            axis=1,
        )
        df['unique'] = df['protein_accessions'].apply(lambda x: "0" if ";" in str(x) else "1")
        df["isotope_label_type"] = 'L'
        df["fragment_ion"] = None
        return df
    
    def iter_batch(self, file_path:str, mods_map:dict, chunksize: int=None):
        for df in pd.read_csv(file_path, sep='\t', usecols=self.use_columns + ['Potential contaminant'], low_memory=False, chunksize=chunksize):
            df = self.main_operate(df, mods_map)
            yield df

    def open_from_zip_archive(self, zip_file, file_name):
        """Open file from zip archive."""
        with zipfile.ZipFile(zip_file) as z:
            with z.open(file_name) as f:
                df = pd.read_csv(f, sep='\t', usecols=self.use_columns + ['Potential contaminant'], low_memory=False)
        return df

    def read_zip_file(self, zip_path:str):
        filepath = Path(zip_path)
        df = self.open_from_zip_archive(zip_path,f'{filepath.stem}/evidence.txt')
        return df

    def merge_sdrf(self, data:pd.DataFrame, sdrf_path:str):
        sdrf = pd.read_csv(sdrf_path, sep="\t")
        factor = "".join(filter(lambda x: x.startswith("factor"), sdrf.columns))
        sdrf = sdrf[
            [
                "comment[data file]",
                "source name",
                factor,
                "comment[fraction identifier]",
                "comment[label]",
                "comment[technical replicate]",
                "characteristics[biological replicate]"
            ]
        ]
        sdrf["comment[data file]"] = sdrf["comment[data file]"].apply(lambda x: x.split(".")[0])
        res = pd.merge(
                data,
                sdrf,
                left_on=["reference_file_name"],
                right_on=["comment[data file]"],
                how="left",
            )
        samples = sdrf["source name"].unique()
        mixed_map = dict(zip(samples, range(1, len(samples) + 1)))
        res.loc[:, "run"] = res[
            [
                "source name",
                "comment[technical replicate]",
                "comment[fraction identifier]",
            ]
        ].apply(
            lambda row: str(mixed_map[row["source name"]])
            + "_"
            + str(row["comment[technical replicate]"])
            + "_"
            + str(row["comment[fraction identifier]"]),
            axis=1,
        )
        res.drop(
            [
                "comment[data file]",
                "comment[technical replicate]",
            ],
            axis=1,
            inplace=True,
        )
        res.rename(columns=self.sdrf_map, inplace=True)
        res.rename(columns={factor: "condition"}, inplace=True)
        return res

    def convert_to_parquet(self, file_path:str, sdrf_path:str, output_path:str, chunksize: int=None):
        self._modifications = get_mods(sdrf_path)
        mods_map = get_mod_map(sdrf_path)
        pqwriter = None
        for df in self.iter_batch(file_path,mods_map,chunksize=chunksize):
            df = self.merge_sdrf(df,sdrf_path)
            df = self.format_to_parquet(df)
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, df.schema)
            pqwriter.write_table(df)
        if pqwriter:
            pqwriter.close()

    def convert_zip_to_parquet(self, files: list, sdrf_path:str,output_path:str):
        self._modifications = get_mods(sdrf_path)
        mods_map = get_mod_map(sdrf_path)
        pqwriter = None
        for file in files:
            try:
                df = self.read_zip_file(file)
                df = self.main_operate(df,mods_map)
                df.loc[:,'reference_file_name'] = file.split('.')[0]
                df = self.merge_sdrf(df,sdrf_path)
                df = self.format_to_parquet(df)
                logging.log(logging.INFO, f"Processing file {file}")
            except Exception as e:
                logging.log(logging.ERROR, f"Error processing file {file}: {e}")
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, df.schema)
            pqwriter.write_table(df)
        if pqwriter:
            pqwriter.close()

    def format_to_parquet(self,res:pd.DataFrame):

        res["sequence"] = res["sequence"].astype(str)
        res["protein_accessions"] = res["protein_accessions"].str.split(";")
        res["unique"] = res["unique"].map(lambda x: None if pd.isna(x) else int(x)).astype("Int32")
        res["modifications"] = res["modifications"].apply(lambda x: self._generate_modification_list(x))
        res["charge"] = res["charge"].map(lambda x: None if pd.isna(x) else int(x)).astype("Int32")
        res["posterior_error_probability"] = res["posterior_error_probability"].astype(float)
        res["intensity"] = res["intensity"].astype(float)
        res["fraction"] = res["fraction"].astype(int).astype(str)
        res["biological_replicate"] = res["biological_replicate"].astype(str)
        res["fragment_ion"] = res["fragment_ion"].astype(str)
        res["run"] = res["run"].astype(str)
        res["condition"] = res["condition"].astype(str)

        return pa.Table.from_pandas(res, schema=SCHEMA)

    def add_additional_msg(self, df:pd.DataFrame):
        # ?????
        df.loc[:,'best_psm_reference_file_name'] = None
        df.loc[:,'best_psm_scan_number'] = None
        df.loc[:,'exp_mass_to_charge'] = None
        df.loc[:,'gene_accessions'] = None
        df.loc[:,'gene_names'] = None
        df.loc[:,'global_qvalue'] = None
        df.loc[:,'id_scores'] = None
        df.loc[:,'intensity_array'] = None
        df.loc[:,'is_decoy'] = None
        df.loc[:,'mz_array'] = None
        df.loc[:,'num_peaks'] = None
        df.loc[:,'protein_end_positions'] = None
        df.loc[:,'protein_global_qvalue'] = None
        df.loc[:,'protein_start_positions'] = None
        df.loc[:,'retention_time'] = None
        df.loc[:,'spectral_count'] = None
        return df
    