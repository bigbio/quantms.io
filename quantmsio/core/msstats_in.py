from quantmsio.core.duckdb import DuckDB
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.core.common import MSSTATS_USECOLS, MSSTATS_MAP
from quantmsio.utils.constants import ITRAQ_CHANNEL, TMT_CHANNELS
from quantmsio.utils.pride_utils import clean_peptidoform_sequence
from quantmsio.operate.tools import get_protein_accession


class MsstatsIN(DuckDB):
    def __init__(self, report_path, sdrf_path, duckdb_max_memory="16GB", duckdb_threads=4):
        super(MsstatsIN, self).__init__(report_path, duckdb_max_memory, duckdb_threads)
        self._sdrf = SDRFHandler(sdrf_path)
        self.experiment_type = self._sdrf.get_experiment_type_from_sdrf()
        self._sample_map = self._sdrf.get_sample_map_run()

    def get_runs(self):
        references = self._duckdb.sql(f"SELECT Reference FROM report").df()
        references = references["Reference"].str.split(".").str[0]
        return list(set(references))

    def iter_runs(self, file_num=10, columns: list = None):
        references = self.get_runs()
        ref_list = [references[i : i + file_num] for i in range(0, len(references), file_num)]
        for refs in ref_list:
            batch_df = self.query_field("Reference", refs, columns)
            yield batch_df

    def generate_msstats_in(self, file_num=10, protein_str=None):
        msstats_map = MSSTATS_MAP.copy()
        usecols = list(MSSTATS_USECOLS)
        if self.experiment_type == "LFQ":
            usecols.remove("Channel")
            usecols.remove("RetentionTime")
            usecols += ["PrecursorCharge"]
            msstats_map["PrecursorCharge"] = "precursor_charge"
        else:
            usecols += ["Charge"]
            msstats_map["Charge"] = "precursor_charge"
        for msstats in self.iter_runs(file_num=file_num, columns=usecols):
            if self.experiment_type == "LFQ":
                msstats.loc[:, "Channel"] = "LFQ"
                msstats.loc[:, "RetentionTime"] = None
            if protein_str:
                msstats = msstats[msstats["ProteinName"].str.contains(f"{protein_str}", na=False)]
            msstats.rename(columns=msstats_map, inplace=True)
            self.transform_msstats_in(msstats)
            self.transform_experiment(msstats)
            yield msstats

    def transform_msstats_in(self, msstats):
        msstats["reference_file_name"] = msstats["reference_file_name"].str.split(".").str[0]
        msstats.loc[:, "sequence"] = msstats["peptidoform"].apply(clean_peptidoform_sequence)
        if self.experiment_type != "LFQ":
            if "TMT" in self.experiment_type:
                msstats["channel"] = msstats["channel"].apply(lambda row: TMT_CHANNELS[self.experiment_type][row - 1])
            else:
                msstats["channel"] = msstats["channel"].apply(lambda row: ITRAQ_CHANNEL[self.experiment_type][row - 1])
        msstats.loc[:, "unique"] = msstats["pg_accessions"].apply(lambda x: 0 if ";" in x or "," in x else 1)
        msstats["pg_accessions"] = msstats["pg_accessions"].apply(get_protein_accession)
        msstats.loc[:, "anchor_protein"] = msstats["pg_accessions"].str[0]

    def transform_experiment(self, msstats):
        intensities_map = {}
        select_cols = ["map", "reference_file_name", "peptidoform", "precursor_charge", "channel", "intensity"]
        if self.experiment_type != "LFQ":
            msstats.loc[:, "map"] = (
                msstats["reference_file_name"] + msstats["peptidoform"] + msstats["precursor_charge"].astype(str)
            )

            def get_intensities_map(rows):
                key = rows["map"]
                sample_key = rows["reference_file_name"] + "-" + rows["channel"]
                if key not in intensities_map:
                    intensities_map[key] = []
                intensities_map[key].append(
                    {
                        "sample_accession": self._sample_map[sample_key],
                        "channel": rows["channel"],
                        "intensity": rows["intensity"],
                    }
                )

            msstats[select_cols].apply(get_intensities_map, axis=1)
            msstats.drop_duplicates(subset=["reference_file_name", "peptidoform", "precursor_charge"], inplace=True)
            msstats.reset_index(inplace=True, drop=True)
            msstats.loc[:, "intensities"] = msstats["map"].map(intensities_map)
            msstats.drop(["map"], inplace=True, axis=1)
        else:
            msstats.loc[:, "intensities"] = msstats[["reference_file_name", "channel", "intensity"]].apply(
                lambda rows: [
                    {
                        "sample_accession": self._sample_map[rows["reference_file_name"] + "-" + rows["channel"]],
                        "channel": rows["channel"],
                        "intensity": rows["intensity"],
                    }
                ],
                axis=1,
            )
