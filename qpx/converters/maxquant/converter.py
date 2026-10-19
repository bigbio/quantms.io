"""MaxQuant orchestrator — composes PSM, Feature, and PG adapters."""

import logging
from pathlib import Path

from qpx._version import __version__
from qpx.converters.mappings import get_tool_meta
from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter
from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM, RUN, SAMPLE
from qpx.core.scores import field_ontology_entries, score_ontology_entries

logger = logging.getLogger(__name__)


class MaxQuantConverter(BaseOrchestrator):
    """Orchestrate full MaxQuant conversion to QPX format."""

    def __init__(self, max_memory=None, max_cpus=None, compression: str = "zstd"):
        self._memory = max_memory or "16GB"
        self._cpus = max_cpus or 4
        self._compression = compression
        self._resolved_mappings_by_view: dict[str, dict] = {}

    def convert(
        self,
        output_folder,
        msms_file=None,
        evidence_file=None,
        protein_groups_file=None,
        sdrf_file=None,
        output_prefix=None,
        structures=None,
        batch_size=100_000,
        n_workers=None,
        spectral_data=False,
        project_accession=None,
        fixed_mod_only: bool = False,
    ):
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        prefix = output_prefix or "maxquant"

        if structures is None:
            structures = []
            if msms_file:
                structures.append(PSM)
            if evidence_file:
                structures.append(FEATURE)
            if protein_groups_file:
                structures.append(PG)

        ontology_entries: list[dict] = []

        # Convert SDRF for sample/run metadata (if SDRF has required columns)
        if sdrf_file:
            try:
                from qpx.converters.sdrf import SdrfConverter

                with SdrfConverter() as sdrf_conv:
                    sdrf_conv.convert(
                        sdrf_path=str(sdrf_file),
                        sample_output=str(output_folder / f"{prefix}.sample.parquet"),
                        run_output=str(output_folder / f"{prefix}.run.parquet"),
                    )
                    ontology_entries.extend(sdrf_conv.run_ontology_entries())
                logger.info("SDRF conversion complete")
            except (ValueError, KeyError, OSError, TypeError) as exc:
                logger.warning("SDRF conversion skipped (incomplete SDRF?): %s", exc)
                # Clean up any corrupt partial files left by the failed write
                for suffix in (".sample.parquet", ".run.parquet"):
                    corrupt = output_folder / f"{prefix}{suffix}"
                    if corrupt.exists():
                        corrupt.unlink()
                        logger.debug("Removed corrupt %s", corrupt)

        duckdb_threads = n_workers if n_workers else self._cpus

        if PSM in structures and msms_file:
            with MaxQuantPsmAdapter(
                duckdb_memory=self._memory,
                duckdb_threads=duckdb_threads,
                compression=self._compression,
            ) as adapter:
                adapter.convert(
                    msms_path=str(msms_file),
                    output_path=str(output_folder / f"{prefix}.psm.parquet"),
                    chunksize=batch_size,
                    spectral_data=spectral_data,
                )
                ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PSM))
                self._resolved_mappings_by_view[PSM] = adapter.get_resolved_columns()
            logger.info("MaxQuant PSM conversion complete")

        if FEATURE in structures and evidence_file:
            with MaxQuantFeatureAdapter(
                duckdb_memory=self._memory,
                duckdb_threads=duckdb_threads,
                compression=self._compression,
            ) as adapter:
                adapter.convert(
                    evidence_path=str(evidence_file),
                    output_path=str(output_folder / f"{prefix}.feature.parquet"),
                    sdrf_path=str(sdrf_file) if sdrf_file else None,
                    protein_groups_path=(str(protein_groups_file) if protein_groups_file else None),
                    chunksize=batch_size,
                    fixed_mod_only=fixed_mod_only,
                )
                ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=FEATURE))
                self._resolved_mappings_by_view[FEATURE] = adapter.get_resolved_columns()
            logger.info("MaxQuant feature conversion complete")

        if PG in structures and protein_groups_file:
            with MaxQuantPgAdapter(
                duckdb_memory=self._memory,
                duckdb_threads=duckdb_threads,
                compression=self._compression,
            ) as adapter:
                adapter.convert(
                    protein_groups_path=str(protein_groups_file),
                    output_path=str(output_folder / f"{prefix}.pg.parquet"),
                    sdrf_path=str(sdrf_file) if sdrf_file else None,
                    evidence_path=str(evidence_file) if evidence_file else None,
                    chunksize=batch_size,
                )
                ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PG))
                self._resolved_mappings_by_view[PG] = adapter.get_resolved_columns()
            logger.info("MaxQuant PG conversion complete")

        for view_name, mappings in self._resolved_mappings_by_view.items():
            ontology_entries.extend(
                field_ontology_entries(
                    view=view_name,
                    resolved_mappings=mappings,
                    tool_name=get_tool_meta("maxquant")["tool_name"],
                )
            )

        self._write_ontology(output_folder, prefix, ontology_entries)
        self._write_provenance(output_folder, prefix, structures)
        self._write_dataset(
            output_folder,
            prefix,
            project_accession,
            software_name="MaxQuant",
            software_version=None,
        )

    def _write_provenance(self, output_folder: Path, prefix: str, structures: list[str]) -> None:
        """Write provenance.parquet with MaxQuant + QPX conversion steps."""
        records = [
            {
                "step_order": 1,
                "step_category": "database_search",
                "step_name": "spectrum_identification",
                "tool_name": "MaxQuant",
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures,
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "maxquant_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures + [SAMPLE, RUN, ONTOLOGY],
            },
        ]
        super()._write_provenance(output_folder, prefix, records)
