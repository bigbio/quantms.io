"""DIA-NN m/z adapter — generate the ``mz`` view from mzML-derived MS2 spectra."""

from __future__ import annotations

import logging
from pathlib import Path

import duckdb
import numpy as np

from qpx.core.spectra_codec import encode_intensity, encode_mz
from qpx.core.sql import escape_path, sql_build
from qpx.writers.mz import MzWriter

logger = logging.getLogger(__name__)


class DiannMzAdapter:
    """Generate the ``mz`` view from mzML-statistics MS2 spectra.

    Reads per-run ``*_ms2_info.parquet`` files (``scan``, ``ms_level``,
    ``mz_array``, ``intensity_array``) produced by the quantmsdiann workflow,
    deduplicates by ``(run_file_name, scan)`` (the run name comes from the file
    stem), losslessly encodes the peak arrays into the ``mz`` / ``intensity``
    binary columns, and writes the ``mz`` view. Only the feature-matched MS2
    scans are kept when ``matched_scans`` is provided.
    """

    def __init__(self, compression: str = "zstd", batch_size: int = 100_000):
        self._compression = compression
        self._batch_size = batch_size

    def convert(
        self,
        mzml_statistics_folder: str,
        output_path: str,
        matched_scans: set[tuple[str, int]] | None = None,
    ) -> int:
        """Write ``{output_path}`` from a folder of ``*_ms2_info.parquet`` files.

        Returns the number of unique spectra written.
        """
        ms2_files = sorted(Path(mzml_statistics_folder).glob("*_ms2_info.parquet"))
        if not ms2_files:
            raise FileNotFoundError(f"no *_ms2_info.parquet files in {mzml_statistics_folder}")

        seen: set[tuple[str, int]] = set()
        pending: list[dict] = []
        written = 0

        con = duckdb.connect()
        writer = MzWriter(output_path, compression=self._compression)
        try:
            for ms2_file in ms2_files:
                run_name = ms2_file.name.removesuffix("_ms2_info.parquet")
                safe_path = escape_path(str(ms2_file))
                cursor = con.execute(
                    sql_build(
                        "SELECT scan, ms_level, mz_array, intensity_array FROM parquet_scan('$path')",
                        path=safe_path,
                    )
                )
                for scan, ms_level, mz, intensity in cursor.fetchall():
                    scan0 = int(scan)
                    key = (run_name, scan0)
                    if key in seen:
                        continue
                    if matched_scans is not None and key not in matched_scans:
                        continue
                    seen.add(key)
                    pending.append(self.build_record(run_name, scan0, ms_level, mz, intensity))
                    if len(pending) >= self._batch_size:
                        writer.write_batch(pending)
                        written += len(pending)
                        pending = []
            if pending:
                writer.write_batch(pending)
                written += len(pending)
            writer.close()
        finally:
            con.close()
        logger.info("wrote %s unique spectra to %s", written, output_path)
        return written

    @staticmethod
    def build_record(run_name: str, scan: int, ms_level: int, mz, intensity) -> dict:
        """Build one ``mz``-view record with losslessly-encoded peak arrays."""
        mz = np.asarray(mz, dtype=np.float32)
        intensity = np.asarray(intensity, dtype=np.float32)
        order = np.argsort(mz, kind="stable")
        mz = mz[order]
        intensity = intensity[order]
        return {
            "id": f"{run_name}:scan={scan}",
            "run_file_name": run_name,
            "scan": scan,
            "ms_level": int(ms_level) if ms_level is not None else 2,
            "centroid": True,
            "scan_start_time": None,
            "inverse_ion_mobility": None,
            "ion_injection_time": None,
            "total_ion_current": None,
            "precursors": [],
            "mz": encode_mz(mz),
            "intensity": encode_intensity(intensity),
            "cv_params": [],
        }
