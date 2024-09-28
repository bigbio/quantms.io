import warnings
from typing import Any
from typing import Tuple

import numpy as np
import pyopenms as oms
from pyopenms import SpectrumLookup


class OpenMSHandler:
    def __init__(self) -> None:
        self._mzml_exp = None
        self._consensus_xml_path = None
        self._spec_lookup = None

    def get_spectrum_from_scan(self, mzml_path: str, scan_number: int) -> Tuple[Any, Any]:
        """
        Get a spectrum from a mzML file using the scan number
        :param mzml_path: path to the mzML file
        :param scan_number: scan number
        :return: spectrum
        """
        if self._mzml_exp is None:
            self._mzml_exp = oms.MSExperiment()
            oms.MzMLFile().load(mzml_path, self._mzml_exp)
            self._spec_lookup = SpectrumLookup()
            self._spec_lookup.readSpectra(self._mzml_exp, "scan=(?<SCAN>\\d+)")
        try:
            index = self._spec_lookup.findByScanNumber(scan_number)
        except IndexError:
            message = "scan_number" + str(scan_number) + "not found in file: " + mzml_path
            warnings.warn(message, category=None, stacklevel=1, source=None)
            return [], []
        spectrum = self._mzml_exp.getSpectrum(index)
        spectrum_mz, spectrum_intensities = spectrum.get_peaks()
        return spectrum_mz, spectrum_intensities

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
            return self._get_intensity_map_tmt_or_itraq(df, experiment_type)
        elif experiment_type is not None and "ITRAQ" in experiment_type.upper():
            return self._get_intensity_map_tmt_or_itraq(df, experiment_type)
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
                if np.float64(row[f"{column}"]) > 0.0:
                    reference_file = column.split(".")[0]
                    key = row.sequence + ":_:" + str(row.charge) + ":_:" + reference_file
                    if key not in intensity_map:
                        intensity_map[key] = {
                            "rt": row.RT,
                            "mz": row.mz,
                            "intensity": row[column],
                        }
                    else:
                        if row[column] > intensity_map[key]["intensity"]:
                            intensity_map[key] = {
                                "rt": row.RT,
                                "mz": row.mz,
                                "intensity": row[column],
                            }
        return intensity_map

    @staticmethod
    def _get_intensity_map_tmt_or_itraq(df, experiment_type):
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
                if np.float64(row[f"{column}"]) > 0.0:
                    reference_file = row.file.split(".")[0]
                    if "TMT" in experiment_type.upper():
                        channel = "TMT" + column.split("_")[1]  # A TMT channel has in consesusXML the following format:
                        # tmt10plex_129N -> TMT129N
                    else:
                        channel = "ITRAQ" + column.split("_")[1]
                    key = row.sequence + ":_:" + str(row.charge) + ":_:" + reference_file + ":_:" + channel
                    if key not in intensity_map:
                        intensity_map[key] = {
                            "rt": row.RT,
                            "mz": row.mz,
                            "intensity": row[column],
                            "channel": channel,
                        }
                    else:
                        if row[column] > intensity_map[key]["intensity"]:
                            intensity_map[key] = {
                                "rt": row.RT,
                                "mz": row.mz,
                                "intensity": row[column],
                                "channel": channel,
                            }
        return intensity_map
