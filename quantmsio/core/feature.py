import pandas as pd
import os
from quantmsio.core.mzTab import MzTab
from quantmsio.core.psm import Psm
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.utils.pride_utils import clean_peptidoform_sequence,get_petidoform_msstats_notation,generate_scan_number
from quantmsio.utils.constants import ITRAQ_CHANNEL,TMT_CHANNELS
from quantmsio.core.common import MSSTATS_MAP,MSSTATS_USECOLS,SDRF_USECOLS,SDRF_MAP
    
class Feature(MzTab):
    def __init__(self,mzTab_path,sdrf_path,msstats_in_path):
        super(Feature,self).__init__(mzTab_path)
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
            pep.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = pep[
                ["modifications", "sequence"]
            ].apply(
                lambda row: get_petidoform_msstats_notation(row["sequence"], row["modifications"], self._modifications),
                axis=1,
            )

        # check spectra_ref
        if "spectra_ref" not in pep.columns:
            pep.loc[:, "scan_number"] = None
            pep.loc[:, "spectra_ref"] = None
        else:
            pep.loc[:, "scan_number"] = pep["spectra_ref"].apply(lambda x: generate_scan_number(x))
            pep["spectra_ref"] = pep["spectra_ref"].apply(lambda x: self._ms_runs[x.split(":")[0]])
        pep_msg = pep.iloc[
            pep.groupby(["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"]).apply(
                lambda row: row["best_search_engine_score[1]"].idxmin(),include_groups=False
            )
        ]
        pep_msg = pep_msg.set_index(["opt_global_cv_MS:1000889_peptidoform_sequence", "charge"])

        pep_msg.loc[:, "pep_msg"] = pep_msg[
            ["best_search_engine_score[1]", "spectra_ref", "scan_number"]
        ].apply(
            lambda row: [
                row["best_search_engine_score[1]"],
                row["spectra_ref"],
                row["scan_number"],
            ],
            axis=1,
        )

        map_dict = pep_msg.to_dict()["pep_msg"]
        return map_dict
    
    def extract_psm_msg(self,chunksize=1000000,protein_str=None):
        P = Psm(self.mztab_path)
        pep_dict = self.extract_from_pep()
        map_dict = {}
        def merge_pep_msg(row):
            key = (row["opt_global_cv_MS:1000889_peptidoform_sequence"],row["precursor_charge"])
            if(key in pep_dict):
                return pep_dict[key]
            else:
                return [None,None,None]
        for psm in P.iter_psm_table(chunksize=chunksize,protein_str=protein_str):
            P.transform_psm(psm)
            if "opt_global_cv_MS:1000889_peptidoform_sequence" not in psm.columns:
                psm.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = psm[
                    ["modifications", "sequence"]
                ].apply(
                    lambda row: get_petidoform_msstats_notation(
                        row["sequence"], row["modifications"], self._modifications
                    ),
                    axis=1,
                )
            psm[["best_qvalue","psm_reference_file_name","psm_scan_number"]] = psm[["opt_global_cv_MS:1000889_peptidoform_sequence", "precursor_charge"]].apply(
                lambda row: merge_pep_msg(row),
                axis=1,
                result_type="expand",
            )
            for key, df in psm.groupby(["opt_global_cv_MS:1000889_peptidoform_sequence", "precursor_charge"]):
                df = df.reset_index(drop=True)
                if key not in map_dict:
                    map_dict[key] = [None for i in range(10)]
                qvalue = None
                temp_df = None
                if(len(df["best_qvalue"].unique()) >1 or not pd.isna(df.loc[0,"best_qvalue"])):
                    temp_df = df.iloc[df["best_qvalue"].idxmin()]
                    qvalue = "best_qvalue"
                elif(len(df["global_qvalue"].unique()) > 1 or not pd.isna(df.loc[0,"best_qvalue"])):
                    temp_df = df.iloc[df["global_qvalue"].idxmin()]
                    qvalue = "global_qvalue"
                #print(temp_df)
                if(qvalue != None):
                    best_qvalue = temp_df[qvalue]
                    if map_dict[key][0] == None or float(map_dict[key][0]) > float(best_qvalue):
                        map_dict[key][0] = temp_df[qvalue]
                        map_dict[key][1] = temp_df["psm_reference_file_name"]
                        map_dict[key][2] = temp_df["psm_scan_number"]
                        map_dict[key][3] = temp_df["pg_positions"]
                        map_dict[key][4] = temp_df["modifications"]
                        map_dict[key][5] = temp_df["posterior_error_probability"]
                        map_dict[key][6] = temp_df["is_decoy"]
                        map_dict[key][7] = temp_df["calculated_mz"]
                        map_dict[key][8] = temp_df["observed_mz"]
                        map_dict[key][9] = temp_df["additional_scores"]
        return map_dict
        
  
    def transform_msstats_in(self,chunksize=1000000,protein_str=None):
        cols = pd.read_csv(self._msstats_in,nrows=0).columns
        if self.experiment_type == 'LFQ':
            MSSTATS_USECOLS.add('PrecursorCharge')
            MSSTATS_MAP['PrecursorCharge'] = 'precursor_charge'
        else:
            MSSTATS_USECOLS.add('Charge')
            MSSTATS_MAP['Charge'] = 'precursor_charge'
        nocols = MSSTATS_USECOLS - set(cols)
        for msstats in pd.read_csv(self._msstats_in,chunksize=chunksize,usecols=list(MSSTATS_USECOLS - set(nocols))):
            if protein_str:
                msstats_in = msstats_in[msstats_in["ProteinName"].str.contains(f"{protein_str}", na=False)]
            for col in nocols:
                if col == "Channel":
                    msstats.loc[:, col] = "LFQ"
                else:
                    msstats.loc[:, col] = None
            msstats["Reference"] = msstats["Reference"].apply(lambda x: x.split(".")[0])
            msstats.loc[:, "sequence"] = msstats["PeptideSequence"].apply(
                lambda x: clean_peptidoform_sequence(x)
            )
            if self.experiment_type != 'LFQ':
                if "TMT" in self.experiment_type:
                    msstats["Channel"] = msstats["Channel"].apply(
                        lambda row: TMT_CHANNELS[self.experiment_type][row - 1]
                    )
                else:
                    msstats["Channel"] = msstats["Channel"].apply(
                        lambda row: ITRAQ_CHANNEL[self.experiment_type][row - 1]
                    )
            msstats.loc[:,"unique"] = msstats['ProteinName'].apply(lambda x: 0 if ';' in x else 1)
            msstats.rename(columns=MSSTATS_MAP,inplace=True)
            yield msstats
    
    def transform_sdrf(self):
        sdrf = pd.read_csv(self._sdrf_path, sep="\t")
        factor = "".join(filter(lambda x: x.startswith("factor"), sdrf.columns))
        SDRF_USECOLS.add(factor)
        sdrf = sdrf[list(SDRF_USECOLS)]
        sdrf["comment[data file]"] = sdrf["comment[data file]"].apply(lambda x: x.split(".")[0])
        samples = sdrf["source name"].unique()
        mixed_map = dict(zip(samples, range(1, len(samples) + 1)))
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
            ],
            axis=1,
            inplace=True,
        )
        sdrf.rename(columns=SDRF_MAP,inplace=True)
        sdrf.rename(columns={factor: "condition"}, inplace=True)
        return sdrf

    def merge_msstats_and_sdrf(self, msstats):
        sdrf = self.transform_sdrf()
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
    
    def merge_msstats_and_psm(self,msstats,map_dict):
        map_features = [
            "global_qvalue",
            "psm_reference_file_name",
            "psm_scan_number",
            "pg_positions",
            "modifications",
            "posterior_error_probability",
            "is_decoy",
            "calculated_mz",
            "observed_mz",
            "additional_scores"
        ]
        for i, feature in enumerate(map_features):
            msstats.loc[:, feature] = msstats[["peptidoform", "precursor_charge"]].apply(
                lambda row: map_dict[(row["peptidoform"], row["precursor_charge"])][i],
                axis=1,
            )
        return msstats