import pyopenms as oms
from scipy.linalg._solve_toeplitz import float64


class ConsensusXMLHandler:

    def __init__(self) -> None:
        self._consensus_xml_path = None

    def get_intensity_map(self, consensusxml_path: str, experiment_type: str = None) -> dict:
        """
        Get the intensity map from a consensusxml file. The intensity map is a dictionary with the following structure:
        - key: peptide sequence + ":_:" + charge + ":_:" + reference file
        - value: dictionary with the following structure:
          - rt: retention time
          - mz: mass to charge ratio
          - intensity: intensity
        :param consensusxml_path: path to the consensusxml file
        :param experiment_type: experiment type (e.g. lfq, tmt, etc.)
        :return: intensity map
        """
        self._consensus_xml_path = consensusxml_path
        consensus_map = oms.ConsensusMap()
        oms.ConsensusXMLFile().load(self._consensus_xml_path, consensus_map)

        df = consensus_map.get_df()
        df = df[df.sequence != "None"]

        if experiment_type is not None and "LABEL FREE" in experiment_type.upper():
            return self._get_intensity_map_lfq(df)
        elif experiment_type is not None and "TMT" in experiment_type.upper():
            return self._get_intensity_map_tmt(df)

        return self._get_intensity_map_lfq(df)  # If not experiment type is provided, we assume it is label free

    @staticmethod
    def _get_intensity_map_lfq(df):
        """
        Get the intensity map for label free experiments
        :param df: pandas dataframe with the consensusxml data
        :return: intensity map
        """
        peptide_columns = ["sequence", "charge", "RT", "mz", "quality"]
        intensity_columns = [column for column in df.columns if column not in peptide_columns]
        intensity_map = {}
        for index, row in df.iterrows():
            for column in intensity_columns:
                if float64(row[f'{column}']) > 0.0:
                    reference_file = column.split(".")[0]
                    key = row.sequence + ":_:" + str(row.charge) + ":_:" + reference_file
                    if key not in intensity_map:
                        intensity_map[key] = {"rt": row.RT, "mz": row.mz, "intensity": row[column]}
                    else:
                        if row[column] > intensity_map[key]["intensity"]:
                            intensity_map[key] = {"rt": row.RT, "mz": row.mz, "intensity": row[column]}
        return intensity_map

    @staticmethod
    def _get_intensity_map_tmt(df):
        """
        Get the intensity map for TMT experiments
        :param df: pandas dataframe with the consensusxml data
        :return: intensity map
        """
        peptide_columns = ["sequence", "charge", "RT", "mz", "quality", "file"]
        intensity_columns = [column for column in df.columns if column not in peptide_columns]
        intensity_map = {}
        for index, row in df.iterrows():
            for column in intensity_columns:
                if float64(row[f'{column}']) > 0.0:
                    reference_file = row.file.split(".")[0]
                    channel = "TMT" + column.split("_")[1]  # A TMT channel has in consesusXML the following format:
                    # tmt10plex_129N -> TMT129N
                    key = row.sequence + ":_:" + str(row.charge) + ":_:" + reference_file + ":_:" + channel
                    if key not in intensity_map:
                        intensity_map[key] = {"rt": row.RT, "mz": row.mz, "intensity": row[column], "channel": channel}
                    else:
                        if row[column] > intensity_map[key]["intensity"]:
                            intensity_map[key] = {"rt": row.RT, "mz": row.mz, "intensity": row[column],
                                                  "channel": channel}
        return intensity_map
