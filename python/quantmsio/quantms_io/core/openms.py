import pyopenms as oms
from scipy.linalg._solve_toeplitz import float64


class ConsensusXMLHandler:

    def __init__(self) -> None:
        self._consensus_xml_path = None

    def get_intensity_map(self, consensusxml_path: str) -> dict:
        """
        Get the intensity map from a consensusxml file. The intensity map is a dictionary with the following structure:
        - key: peptide sequence + ":_:" + charge + ":_:" + reference file
        - value: dictionary with the following structure:
          - rt: retention time
          - mz: mass to charge ratio
          - intensity: intensity
        :param consensusxml_path: path to the consensusxml file
        :return: intensity map
        """
        self._consensus_xml_path = consensusxml_path
        consensus_map = oms.ConsensusMap()
        oms.ConsensusXMLFile().load(self._consensus_xml_path, consensus_map)

        df = consensus_map.get_df()
        df = df[df.sequence != "None"]
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
