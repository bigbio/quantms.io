"""MzSpectra data structure — mass spectrometry spectral data."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema

MzSchema = load_schema("mz")


class MzSpectra(BaseStructure):
    """Mass spectrometry spectral data (scan-level)."""

    _schema_class = MzSchema

    def ms1(self) -> "MzSpectra":
        """Filter to MS1 scans only."""
        return self.filter("ms_level = 1")

    def ms2(self) -> "MzSpectra":
        """Filter to MS2 scans only."""
        return self.filter("ms_level = 2")

    def by_rt_range(self, rt_start: float, rt_end: float) -> "MzSpectra":
        """Filter scans by retention time range (minutes)."""
        return self.filter(f"scan_start_time >= {rt_start} AND scan_start_time <= {rt_end}")
