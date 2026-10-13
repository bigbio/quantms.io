"""DIA-NN m/z adapter — generate the ``mz`` view from feature-matched MS2 spectra."""

from __future__ import annotations

import logging
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from qpx.converters.diann.feature_adapter import _ms_info_scans
from qpx.core.sql import escape_path, sql_build
from qpx.writers.mz import MzWriter

logger = logging.getLogger(__name__)

# Rows fetched from DuckDB at a time, so only a bounded slice of peak arrays is in Python.
_FETCH_ROWS = 500


def _scan_numbers(con: duckdb.DuckDBPyConnection, path: Path) -> pd.DataFrame:
    """Distinct ``scan`` values of a stats table (``raw``) with their scan numbers (``scan``).

    Scan numbers are resolved as the feature view resolves them; IDs without a
    single scan number are dropped.
    """
    raw = con.execute(sql_build("SELECT DISTINCT scan FROM parquet_scan('$path')", path=escape_path(str(path)))).df()["scan"]
    scans = pd.DataFrame({"raw": raw, "scan": _ms_info_scans(raw)}).dropna(subset=["scan"])
    return scans.astype({"scan": "int64"})


def _rt_lookup(con: duckdb.DuckDBPyConnection, stats: Path) -> dict[str, dict[int, float]]:
    """Map run -> {scan: rt} from the per-run ``*_ms_info.parquet`` tables."""
    lookup: dict[str, dict[int, float]] = {}
    for path in sorted(stats.glob("*_ms_info.parquet")):
        table = con.execute(
            sql_build("SELECT scan, rt FROM parquet_scan('$path') WHERE rt IS NOT NULL", path=escape_path(str(path)))
        ).df()
        scans = _ms_info_scans(table["scan"])
        known = scans.notna()
        run_name = path.name.removesuffix("_ms_info.parquet")
        lookup[run_name] = dict(zip(scans[known].astype("int64").tolist(), table.loc[known, "rt"].astype(float).tolist()))
    return lookup


def _matched_spectra(con: duckdb.DuckDBPyConnection, ms2_file: Path, run_name: str, matched_scans):
    """Yield ``(scan, ms_level, mz_array, intensity_array)`` once per matched scan of one run.

    DuckDB applies the scan filter and rows are fetched in chunks, so only the
    matched spectra reach Python. ``matched_scans`` of ``None`` keeps every scan.
    """
    keep = _scan_numbers(con, ms2_file)
    if matched_scans is not None:
        keep = keep[keep["scan"].isin({scan for run, scan in matched_scans if run == run_name})]
    con.register("keep_scans", keep.drop_duplicates("scan"))
    try:
        cursor = con.execute(
            sql_build(
                "SELECT k.scan, m.ms_level, m.mz_array, m.intensity_array "
                "FROM parquet_scan('$path') m JOIN keep_scans k ON m.scan = k.raw",
                path=escape_path(str(ms2_file)),
            )
        )
        seen: set[int] = set()
        while rows := cursor.fetchmany(_FETCH_ROWS):
            for scan, ms_level, mz, intensity in rows:
                if scan not in seen:
                    seen.add(scan)
                    yield scan, ms_level, mz, intensity
    finally:
        con.unregister("keep_scans")


class DiannMzAdapter:
    """Generate the ``mz`` view from mzML-statistics MS2 spectra.

    Reads per-run ``*_ms2_info.parquet`` files (``scan``, ``ms_level``,
    ``mz_array``, ``intensity_array``) produced by the quantmsdiann workflow,
    keeping only the feature-matched MS2 scans (``matched_scans``) and writing
    their peak arrays to the ``mz`` / ``intensity`` float32 list columns.
    ``scan_start_time`` is filled from the per-run ``*_ms_info.parquet`` ``rt``.
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
        """Write ``{output_path}`` from ``*_ms2_info.parquet`` / ``*_ms_info.parquet``.

        Returns the number of unique feature-matched spectra written.
        """
        stats = Path(mzml_statistics_folder)
        ms2_files = sorted(stats.glob("*_ms2_info.parquet"))
        if not ms2_files:
            raise FileNotFoundError(f"no *_ms2_info.parquet files in {mzml_statistics_folder}")

        con = duckdb.connect()
        rt_lookup = _rt_lookup(con, stats)

        pending: list[dict] = []
        written = 0
        writer = MzWriter(output_path, compression=self._compression, batch_size=self._batch_size)
        try:
            for ms2_file in ms2_files:
                run_name = ms2_file.name.removesuffix("_ms2_info.parquet")
                run_rt = rt_lookup.get(run_name, {})
                for scan, ms_level, mz, intensity in _matched_spectra(con, ms2_file, run_name, matched_scans):
                    pending.append(self.build_record(run_name, scan, ms_level, mz, intensity, run_rt.get(scan)))
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
        logger.info("wrote %s unique feature-matched spectra to %s", written, output_path)
        return written

    @staticmethod
    def build_record(run_name: str, scan: int, ms_level: int, mz, intensity, rt: float | None = None) -> dict:
        """Build one ``mz``-view record with float32 peak arrays."""
        mz = np.asarray(mz, dtype=np.float32)
        intensity = np.asarray(intensity, dtype=np.float32)
        return {
            "id": f"{run_name}:scan={scan}",
            "run_file_name": run_name,
            "scan": [scan],
            "ms_level": int(ms_level) if ms_level is not None else 2,
            "centroid": True,
            "scan_start_time": (float(rt) / 60.0) if rt is not None else 0.0,
            "inverse_ion_mobility": None,
            "ion_injection_time": 0.0,
            "total_ion_current": float(intensity.sum()) if len(intensity) else 0.0,
            "precursors": [],
            "mz": mz,
            "intensity": intensity,
            "cv_params": [],
        }
