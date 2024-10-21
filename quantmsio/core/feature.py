import pandas as pd
import os
import pyarrow as pa
import pyarrow.parquet as pq
from quantmsio.utils.file_utils import extract_protein_list
from quantmsio.core.mztab import MzTab, generate_modification_list
from quantmsio.core.psm import Psm
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.utils.pride_utils import (
    clean_peptidoform_sequence,
    get_petidoform_msstats_notation,
    generate_scan_number,
    get_peptidoform_proforma_version_in_mztab,
)
from quantmsio.utils.constants import ITRAQ_CHANNEL, TMT_CHANNELS
from quantmsio.core.common import MSSTATS_MAP, MSSTATS_USECOLS, SDRF_USECOLS, SDRF_MAP, QUANTMSIO_VERSION
from quantmsio.core.format import FEATURE_FIELDS

FEATURE_SCHEMA = pa.schema(
    FEATURE_FIELDS,
    metadata={"description": "feature file in quantms.io format"},
)


class Feature(MzTab):
    def __init__(self, mzTab_path, sdrf_path, msstats_in_path):
        super(Feature, self).__init__(mzTab_path)
        self._msstats_in = msstats_in_path
        self._sdrf_path = sdrf_path
        self._ms_runs = self.extract_ms_runs()
        self._protein_global_qvalue_map = self.get_protein_map()
        self._modifications = self.get_modifications()
        self._score_names = self.get_score_names()
        self.experiment_type = SDRFHandler(sdrf_path).get_experiment_type_from_sdrf()

    def _extract_pep_columns(self):
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = open(self.mztab_path)
        line = f.readline()
        while not line.startswith("PEH"):
            line = f.readline()
        self._pep_columns = line.split("\n")[0].split("\t")

    def extract_from_pep(self):
        self._extract_pep_columns()
        pep_usecols = [
            "opt_global_cv_MS:1000889_peptidoform_sequence",
            "charge",
            "best_search_engine_score[1]",
            "spectra_ref",
        ]
        live_cols = [col for col in pep_usecols if col in self._pep_columns]
        not_cols = [col for col in pep_usecols if col not in live_cols]
        if "opt_global_cv_MS:1000889_peptidoform_sequence" in not_cols:
            if "sequence" in self._pep_columns and "modifications" in self._pep_columns:
                live_cols.append("sequence")
                live_cols.append("modifications")
            else:
                raise Exception("The peptide table don't have opt_global_cv_MS:1000889_peptidoform_sequence columns")
        if "charge" in not_cols or "best_search_engine_score[1]" in not_cols:
            raise Exception("The peptide table don't have best_search_engine_score[1] or charge columns")

        pep = self.skip_and_load_csv("PEH", usecols=live_cols)

        if "opt_global_cv_MS:1000889_peptidoform_sequence" not in pep.columns:
            pep.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = pep[["modifications", "sequence"]].apply(
                lambda row: get_petidoform_msstats_notation(row["sequence"], row["modifications"], self._modifications),
                axis=1,
            )

        # check spectra_ref
        if "spectra_ref" not in pep.columns:
            pep.loc[:, "scan_number"] = None
            pep.loc[:, "spectra_ref"] = None
        else:
            pep.loc[:, "scan_number"] = pep["spectra_ref"].apply(generate_scan_number)
            pep["spectra_ref"] = pep["spectra_ref"].apply(lambda x: self._ms_runs[x.split(":")[0]])
        pep_msg = pep.iloc[
            pep.groupby(["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"]).apply(
                lambda row: row["best_search_engine_score[1]"].idxmin()
            )
        ]
        pep_msg = pep_msg.set_index(["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"])

        pep_msg.loc[:, "pep_msg"] = pep_msg[["best_search_engine_score[1]", "spectra_ref", "scan_number"]].apply(
            lambda row: [
                row["best_search_engine_score[1]"],
                row["spectra_ref"],
                row["scan_number"],
            ],
            axis=1,
        )

        map_dict = pep_msg.to_dict()["pep_msg"]
        return map_dict

    def extract_psm_msg(self, chunksize=1000000, protein_str=None):
        P = Psm(self.mztab_path)
        pep_dict = self.extract_from_pep()
        map_dict = {}

        def merge_pep_msg(row):
            key = (row["opt_global_cv_MS:1000889_peptidoform_sequence"], row["precursor_charge"])
            if key in pep_dict:
                return pep_dict[key]
            else:
                return [None, None, None]

        for psm in P.iter_psm_table(chunksize=chunksize, protein_str=protein_str):
            P.transform_psm(psm)
            if "opt_global_cv_MS:1000889_peptidoform_sequence" not in psm.columns:
                psm.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = psm[["modifications", "sequence"]].apply(
                    lambda row: get_petidoform_msstats_notation(
                        row["sequence"], row["modifications"], self._modifications
                    ),
                    axis=1,
                )
            psm[["best_qvalue", "scan_reference_file_name", "scan"]] = psm[
                ["opt_global_cv_MS:1000889_peptidoform_sequence", "precursor_charge"]
            ].apply(
                merge_pep_msg,
                axis=1,
                result_type="expand",
            )
            for key, df in psm.groupby(["opt_global_cv_MS:1000889_peptidoform_sequence", "precursor_charge"]):
                df = df.reset_index(drop=True)
                if key not in map_dict:
                    map_dict[key] = [None for _ in range(10)]
                qvalue = None
                temp_df = None
                if len(df["best_qvalue"].unique()) > 1 or not pd.isna(df.loc[0, "best_qvalue"]):
                    temp_df = df.iloc[df["best_qvalue"].idxmin()]
                    qvalue = "best_qvalue"
                elif len(df["global_qvalue"].unique()) > 1 or not pd.isna(df.loc[0, "best_qvalue"]):
                    temp_df = df.iloc[df["global_qvalue"].idxmin()]
                    qvalue = "global_qvalue"
                if qvalue is not None:
                    best_qvalue = temp_df[qvalue]
                    if map_dict[key][0] is None or float(map_dict[key][0]) > float(best_qvalue):
                        map_dict[key][0] = temp_df[qvalue]
                        map_dict[key][1] = temp_df["scan_reference_file_name"]
                        map_dict[key][2] = temp_df["scan"]
                        map_dict[key][3] = temp_df["mp_accessions"]
                        map_dict[key][4] = temp_df["modifications"]
                        map_dict[key][5] = temp_df["posterior_error_probability"]
                        map_dict[key][6] = temp_df["is_decoy"]
                        map_dict[key][7] = temp_df["calculated_mz"]
                        map_dict[key][8] = temp_df["observed_mz"]
                        map_dict[key][9] = temp_df["additional_scores"]
        return map_dict

    def transform_msstats_in(self, chunksize=1000000, protein_str=None):
        cols = pd.read_csv(self._msstats_in, nrows=0).columns
        if self.experiment_type == "LFQ":
            MSSTATS_USECOLS.add("PrecursorCharge")
            MSSTATS_MAP["PrecursorCharge"] = "precursor_charge"
        else:
            MSSTATS_USECOLS.add("Charge")
            MSSTATS_MAP["Charge"] = "precursor_charge"
        nocols = MSSTATS_USECOLS - set(cols)
        for msstats in pd.read_csv(self._msstats_in, chunksize=chunksize, usecols=list(MSSTATS_USECOLS - set(nocols))):
            if protein_str:
                msstats = msstats[msstats["ProteinName"].str.contains(f"{protein_str}", na=False)]
            for col in nocols:
                if col == "Channel":
                    msstats.loc[:, col] = "LFQ"
                else:
                    msstats.loc[:, col] = None
            msstats["Reference"] = msstats["Reference"].apply(lambda x: x.split(".")[0])
            msstats.loc[:, "sequence"] = msstats["PeptideSequence"].apply(clean_peptidoform_sequence)
            if self.experiment_type != "LFQ":
                if "TMT" in self.experiment_type:
                    msstats["Channel"] = msstats["Channel"].apply(
                        lambda row: TMT_CHANNELS[self.experiment_type][row - 1]
                    )
                else:
                    msstats["Channel"] = msstats["Channel"].apply(
                        lambda row: ITRAQ_CHANNEL[self.experiment_type][row - 1]
                    )
            msstats.loc[:, "unique"] = msstats["ProteinName"].apply(lambda x: 0 if ";" in x else 1)
            msstats.rename(columns=MSSTATS_MAP, inplace=True)
            yield msstats

    @staticmethod
    def transform_sdrf(sdrf_path):
        sdrf = pd.read_csv(sdrf_path, sep="\t")
        factor = list(filter(lambda x: x.startswith("factor"), sdrf.columns))
        usecols = list(SDRF_USECOLS) + factor
        sdrf = sdrf[usecols]
        sdrf["comment[data file]"] = sdrf["comment[data file]"].apply(lambda x: x.split(".")[0])
        samples = sdrf["source name"].unique()
        mixed_map = dict(zip(samples, range(1, len(samples) + 1)))
        sdrf.loc[:, 'conditions'] = sdrf[factor].apply(lambda row:[str(row[col]) for col in factor],axis=1)
        sdrf.loc[:, "run"] = sdrf[
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
        sdrf.drop(
            [
                "comment[technical replicate]",
            ] + factor,
            axis=1,
            inplace=True,
        )
        sdrf.rename(columns=SDRF_MAP, inplace=True)
        return sdrf

    def merge_msstats_and_sdrf(self, msstats):
        sdrf = self.transform_sdrf(self._sdrf_path)
        if self.experiment_type != "LFQ":
            msstats = pd.merge(
                msstats,
                sdrf,
                left_on=["reference_file_name", "channel"],
                right_on=["reference", "label"],
                how="left",
            )
        else:
            msstats = pd.merge(
                msstats,
                sdrf,
                left_on=["reference_file_name"],
                right_on=["reference"],
                how="left",
            )
        msstats.drop(
            [
                "reference",
                "label",
            ],
            axis=1,
            inplace=True,
        )
        return msstats

    def merge_msstats_and_psm(self, msstats, map_dict):
        map_features = [
            "global_qvalue",
            "scan_reference_file_name",
            "scan",
            "mp_accessions",
            "modifications",
            "posterior_error_probability",
            "is_decoy",
            "calculated_mz",
            "observed_mz",
            "additional_scores",
        ]
        for i, feature in enumerate(map_features):
            msstats.loc[:, feature] = msstats[["peptidoform", "precursor_charge"]].apply(
                lambda row: map_dict[(row["peptidoform"], row["precursor_charge"])][i],
                axis=1,
            )
        return msstats

    def generate_feature(self, chunksize=1000000, protein_str=None):
        map_dict = self.extract_psm_msg(chunksize, protein_str)
        for msstats in self.transform_msstats_in(chunksize, protein_str):
            msstats = self.merge_msstats_and_sdrf(msstats)
            msstats = self.merge_msstats_and_psm(msstats, map_dict)
            self.add_additional_msg(msstats)
            self.convert_to_parquet_format(msstats, self._modifications)
            feature = self.transform_feature(msstats)
            yield feature

    @staticmethod
    def slice(df, partitions):
        cols = df.columns
        if not isinstance(partitions, list):
            raise Exception(f"{partitions} is not a list")
        if len(partitions) == 0:
            raise Exception(f"{partitions} is empty")
        for partion in partitions:
            if partion not in cols:
                raise Exception(f"{partion} does not exist")
        for key, df in df.groupby(partitions):
            yield key, df

    def generate_slice_feature(self, partitions, chunksize=1000000, protein_str=None):
        map_dict = self.extract_psm_msg(chunksize, protein_str)
        for msstats in self.transform_msstats_in(chunksize, protein_str):
            msstats = self.merge_msstats_and_sdrf(msstats)
            msstats = self.merge_msstats_and_psm(msstats, map_dict)
            self.add_additional_msg(msstats)
            self.convert_to_parquet_format(msstats, self._modifications)
            for key, df in self.slice(msstats, partitions):
                feature = self.transform_feature(df)
                yield key, feature

    @staticmethod
    def transform_feature(df):
        return pa.Table.from_pandas(df, schema=FEATURE_SCHEMA)

    def write_feature_to_file(
        self,
        output_path,
        chunksize=1000000,
        protein_file=None,
    ):
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None
        pqwriter = None
        for feature in self.generate_feature(chunksize, protein_str):
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, feature.schema)
            pqwriter.write_table(feature)
        if pqwriter:
            pqwriter.close()

    def write_features_to_file(self, output_folder, filename, partitions, chunksize=1000000, protein_file=None):
        pqwriters = {}
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None
        for key, feature in self.generate_slice_feature(partitions, chunksize, protein_str):
            folder = [output_folder] + [str(col) for col in key]
            folder = os.path.join(*folder)
            if not os.path.exists(folder):
                os.makedirs(folder, exist_ok=True)
            save_path = os.path.join(*[folder, filename])
            if not os.path.exists(save_path):
                pqwriter = pq.ParquetWriter(save_path, feature.schema)
                pqwriters[key] = pqwriter
            pqwriters[key].write_table(feature)

        for pqwriter in pqwriters.values():
            pqwriter.close()

    def add_additional_msg(self, msstats):
        msstats.loc[:, "pg_global_qvalue"] = msstats["pg_accessions"].map(self._protein_global_qvalue_map)
        msstats.loc[:, "peptidoform"] = msstats[["modifications", "sequence"]].apply(
            lambda row: get_peptidoform_proforma_version_in_mztab(
                row["sequence"], row["modifications"], self._modifications
            ),
            axis=1,
        )
        msstats.loc[:, "additional_intensities"] = None
        msstats.loc[:, "best_id_score"] = None
        msstats.loc[:, "modification_details"] = None
        msstats.loc[:, "predicted_rt"] = None
        msstats.loc[:, "gg_accessions"] = None
        msstats.loc[:, "gg_names"] = None
        msstats.loc[:, "cv_params"] = None
        msstats.loc[:, "rt_start"] = None
        msstats.loc[:, "rt_stop"] = None

    @staticmethod
    def convert_to_parquet_format(res, modifications):
        res["pg_accessions"] = res["pg_accessions"].str.split(";")
        res["mp_accessions"] = res["mp_accessions"].str.split(";")
        res["pg_global_qvalue"] = res["pg_global_qvalue"].astype(float)
        res["unique"] = res["unique"].astype("Int32")
        res["modifications"] = res["modifications"].apply(lambda x: generate_modification_list(x, modifications))
        res["precursor_charge"] = res["precursor_charge"].map(lambda x: None if pd.isna(x) else int(x)).astype("Int32")
        res["calculated_mz"] = res["calculated_mz"].astype(float)
        res["observed_mz"] = res["observed_mz"].astype(float)
        res["posterior_error_probability"] = res["posterior_error_probability"].astype(float)
        res["global_qvalue"] = res["global_qvalue"].astype(float)
        res["is_decoy"] = res["is_decoy"].map(lambda x: None if pd.isna(x) else int(x)).astype("Int32")
        res["fraction"] = res["fraction"].astype(int).astype(str)
        res["biological_replicate"] = res["biological_replicate"].astype(str)
        res["scan"] = res["scan"].astype(str)
        res["scan_reference_file_name"] = res["scan_reference_file_name"].astype(str)
        if "rt" in res.columns:
            res["rt"] = res["rt"].astype(float)
        else:
            res.loc[:, "rt"] = None
        # return pa.Table.from_pandas(res, schema=FEATURE_SCHEMA)
