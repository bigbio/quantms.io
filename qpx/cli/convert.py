"""
Convert subcommands — convert external tool outputs to QPX format.

Each command wraps a converter from ``qpx.converters.*``.

Subcommands:
    qpxc convert diann      — DIA-NN report to QPX
    qpxc convert maxquant   — MaxQuant output to QPX
    qpxc convert fragpipe   — FragPipe output to QPX
    qpxc convert mzidentml  — mzIdentML (incl. XL-MS 1.3) to QPX
    qpxc convert openms     — OpenMS -out_qpx enrichment to full QPX
    qpxc convert spectronaut — Spectronaut report to QPX
    qpxc convert cdap       — CPTAC CDAP .psm files to QPX
    qpxc convert sdrf       — SDRF to sample.parquet + run.parquet
    qpxc convert quantms-msstats — QuantMS MSstats + SDRF to QPX
"""

from __future__ import annotations

import functools
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import click

from qpx.converters.openms import OpenMSConverter

logger = logging.getLogger("qpx.cli.convert")


def _warn(message: str) -> None:
    """Print a warning to stderr with the shared ``WARNING:`` prefix."""
    click.echo(f"WARNING: {message}", err=True)


def _log_summary(output_folder) -> None:
    """Log a conversion summary for a just-written output folder.

    Shared by every convert subcommand. Guarded so a summary failure NEVER fails
    an otherwise-successful conversion (a warning is logged instead).
    """
    from qpx.converters.summary import log_conversion_summary

    log_conversion_summary(output_folder, logger=logger)


def _write_mudata(output_folder: Path, prefix: str, enabled: bool) -> None:
    """Refresh or remove the dataset's MuData view after a conversion.

    Uses qpx's own writer, which refuses a MuData missing a required
    quantification modality rather than writing a partial view, and is
    best-effort: the parquet views are the dataset's source of truth, so a view
    that cannot be built is reported and the conversion still succeeds.
    """
    if not enabled:
        try:
            (Path(output_folder) / f"{prefix}.h5mu").unlink(missing_ok=True)
        except OSError as exc:
            logger.warning("Could not remove stale MuData for %s: %s", prefix, exc)
        return
    from qpx.mudata import write_dataset_mudata

    written = write_dataset_mudata(Path(output_folder), prefix)
    if written is None:
        _warn("no MuData view was written (see the log); the Parquet views are complete")
    else:
        click.echo(f"MuData view: {written.name}")


def _annotate_protein_properties(output_folder: Path, fasta: Optional[Path], prefix: Optional[str] = None) -> None:
    """Fill null protein properties from an optional FASTA after a conversion.

    The conversion has already succeeded, so a dataset this step cannot handle
    (e.g. partitioned views) is reported and left as converted rather than failed.
    """
    if fasta is None:
        return
    from qpx.cli.transform import annotate_dataset_protein_properties

    click.echo(f"Filling protein properties from {fasta.name}")
    try:
        annotate_dataset_protein_properties(Path(output_folder), fasta, in_place=True, prefix=prefix)
    except click.ClickException as exc:
        _warn(f"protein properties from FASTA skipped: {exc.format_message()}")


def _maybe_enrich_pride(output_folder, project_accession: str | None, enrich: bool, prefix: str) -> None:
    """Optionally enrich a converted dataset with PRIDE metadata."""
    if not enrich:
        return
    if not project_accession:
        _warn("--enrich-pride requires --project-accession; skipping enrichment.")
        return
    try:
        from qpx.dataset import Dataset

        with Dataset(str(output_folder), file_prefix=prefix) as ds:
            ds.enrich_from_pride(project_accession)
        click.echo(f"Enriched dataset with PRIDE metadata for {project_accession}")
    except Exception as exc:
        _warn(f"Could not enrich from PRIDE: {exc}")


@dataclass(frozen=True)
class _DatasetMetadata:
    """Validated dataset-metadata options of one convert command.

    ``config_json`` holds the workflow parameters as compact JSON, or None when no config file is given.
    """

    software_name: str | None
    software_version: str | None
    config_json: str | None
    enrich_pride: bool
    compute_integrity: bool


def _dataset_metadata_options(func):
    """Add the dataset-metadata options shared by the convert commands.

    The options set the workflow name and version, enrich from PRIDE and record
    file integrity. They are listed in the order declared here. The command receives
    them as one ``metadata`` argument.
    """
    options = (
        click.option(
            "--software-name",
            default=None,
            help=(
                "Workflow that produced the input (e.g. nf-core/mhcquant). "
                "Sets dataset.software_name/software_version and adds a 'workflow' provenance step."
            ),
        ),
        click.option("--software-version", default=None, help="Version of --software-name."),
        click.option(
            "--software-config",
            default=None,
            type=click.Path(exists=True, dir_okay=False, path_type=Path),
            help=(
                "JSON file with the workflow parameters (e.g. Nextflow params.json), "
                "stored in the workflow provenance step's config."
            ),
        ),
        click.option(
            "--enrich-pride",
            is_flag=True,
            default=False,
            help=(
                "Fetch project title, description and PubMed ID from the PRIDE API after "
                "conversion (needs --project-accession; best-effort)."
            ),
        ),
        click.option(
            "--compute-integrity",
            is_flag=True,
            default=False,
            help="Write file checksums, row counts and sizes for the output folder into dataset.parquet, as the last step.",
        ),
    )

    @functools.wraps(func)
    def wrapper(*, software_name, software_version, software_config, enrich_pride, compute_integrity, **kwargs):
        """Validate the workflow options, then call the command with them as ``metadata``."""
        metadata = _DatasetMetadata(
            software_name=software_name,
            software_version=software_version,
            config_json=_software_config_json(software_name, software_version, software_config),
            enrich_pride=enrich_pride,
            compute_integrity=compute_integrity,
        )
        return func(metadata=metadata, **kwargs)

    for option in reversed(options):
        wrapper = option(wrapper)
    return wrapper


def _software_config_json(software_name, software_version, software_config) -> str | None:
    """Validate the workflow options and return the parameters as compact, sorted JSON.

    Runs before the conversion, so an invalid option fails without writing output.
    """
    if (software_version is not None or software_config is not None) and not software_name:
        raise click.UsageError("--software-version and --software-config require --software-name.")
    if software_config is None:
        return None
    try:
        with open(software_config, encoding="utf-8") as handle:
            config = json.load(handle)
    except (OSError, ValueError) as exc:
        raise click.UsageError(f"--software-config {software_config} is not a readable JSON file: {exc}") from exc
    return json.dumps(config, separators=(",", ":"), sort_keys=True)


def _apply_software(
    output_folder,
    prefix: str,
    software_name: str | None,
    software_version: str | None,
    config_json: str | None,
) -> None:
    """Record the workflow that produced the dataset in dataset and provenance parquet.

    Sets the software fields of dataset.parquet and prepends a 'workflow' step to
    provenance.parquet. Skipped with a warning when either metadata file is missing.
    """
    if not software_name:
        return
    import pandas as pd

    from qpx.converters.orchestrator import build_provenance_step
    from qpx.dataset import Dataset

    with Dataset(str(output_folder), file_prefix=prefix) as ds:
        if ds.dataset_meta is None or ds.provenance is None:
            _warn("no dataset/provenance metadata found; --software-name skipped")
            return
        record = ds.dataset_meta.to_df().iloc[0].to_dict()
        record["software_name"] = software_name
        record["software_version"] = software_version
        ds.save_structure([record], "dataset", prefix=prefix)

        workflow = build_provenance_step(
            step_category="workflow",
            step_name="pipeline_execution",
            tool_name=software_name,
            tool_version=software_version,
            config=config_json,
        )
        frame = pd.concat([pd.DataFrame([workflow]), ds.provenance.to_df()], ignore_index=True)
        frame["step_order"] = range(1, len(frame) + 1)
        ds.save_structure(frame, "provenance", prefix=prefix)


def _maybe_compute_integrity(output_folder, prefix: str, enabled: bool) -> None:
    """Record checksums, row counts and sizes of the output folder in dataset parquet.

    Call this after the last file write of a conversion, so the record covers every Parquet and h5ad file.
    """
    if not enabled:
        return
    from qpx.dataset import Dataset

    with Dataset(str(output_folder), file_prefix=prefix) as ds:
        if ds.dataset_meta is None:
            _warn("no dataset.parquet found; --compute-integrity skipped")
            return
        integrity = ds.compute_integrity()
        record = ds.dataset_meta.to_df().iloc[0].to_dict()
        record.update(integrity)
        ds.save_structure([record], "dataset", prefix=prefix)
    click.echo(f"Integrity recorded for {integrity['total_structures']} files")


def _finish_conversion(
    output_folder: Path,
    prefix: str,
    project_accession: str | None,
    metadata: _DatasetMetadata,
    *,
    fasta: Path | None = None,
    mudata: bool | None = None,
    written: bool = True,
) -> None:
    """Run the post-conversion steps shared by the convert commands, in order.

    Dataset-changing steps run only when ``written`` is true. The MuData view is
    refreshed when ``mudata`` is not None. The integrity record runs last, so it covers every
    Parquet and h5ad file (not the ``.h5mu`` view).
    """
    if written:
        _annotate_protein_properties(output_folder, fasta, prefix)
        _apply_software(output_folder, prefix, metadata.software_name, metadata.software_version, metadata.config_json)
        _maybe_enrich_pride(output_folder, project_accession, metadata.enrich_pride, prefix)
    if mudata is not None:
        _write_mudata(output_folder, prefix, mudata and written)
    if written:
        _maybe_compute_integrity(output_folder, prefix, metadata.compute_integrity)
    _log_summary(output_folder)


# ---------------------------------------------------------------------------
# Convert group
# ---------------------------------------------------------------------------


@click.group()
def convert():
    """Convert external tool outputs to QPX format."""
    pass


# ---------------------------------------------------------------------------
# DIA-NN
# ---------------------------------------------------------------------------


@convert.command("diann")
@click.option(
    "--report-path",
    help="DIA-NN report file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--mzml-info-folder",
    help="Folder containing mzML info files (optional; scan/mz fields left empty if omitted)",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--qvalue-threshold",
    help=(
        "Optional q-value threshold. Unset (the default) converts the DIA-NN "
        "report as reported — no filtering — since DIA-NN already FDR-filters its "
        "main report and every per-row q-value column is carried through for "
        "downstream filtering. When given, the feature and pg views keep only precursors "
        "with Q.Value below it; protein-group q-values are carried, not filtered."
    ),
    default=None,
    type=float,
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option("--output-prefix", help="Prefix for output file names")
@click.option(
    "--pg-matrix-path",
    help="DIA-NN protein quantities matrix file (enables PG conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--protein-file",
    help="Protein file for filtering",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--partitions",
    help="Field(s) for splitting output files (comma-separated)",
)
@click.option(
    "--max-memory",
    help="Maximum memory limit (e.g. '16GB')",
)
@click.option(
    "--max-cpus",
    help="Maximum number of threads",
    type=int,
)
@click.option(
    "--batch-size",
    help="Number of files to process simultaneously",
    default=100,
    type=int,
)
@click.option(
    "--standardized-intensities",
    help="Calculate standardized intensity metrics for PG output",
    is_flag=True,
    default=False,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD020192)",
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option(
    "--diann-log",
    help="DIA-NN summary log file (version auto-detected from first line)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--fasta",
    help=(
        "Optional protein FASTA used for the search. When given, fills null pg.sequence_coverage, "
        "pg.molecular_weight and feature.pg_positions for target rows after conversion; proteins "
        "absent from it (e.g. DIA-NN internal decoys) stay null. See 'qpxc transform protein-properties'."
    ),
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--mudata/--no-mudata",
    default=True,
    show_default=True,
    help=(
        "Write the dataset's MuData (.h5mu) view after conversion. Best-effort: a view that "
        "cannot be built is reported and the conversion still succeeds. Use --no-mudata to skip "
        "the build on very large datasets."
    ),
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@_dataset_metadata_options
def convert_diann_cmd(
    report_path: Path,
    sdrf_file: Path,
    mzml_info_folder: Path,
    qvalue_threshold: Optional[float],
    output_folder: Path,
    output_prefix: Optional[str],
    pg_matrix_path: Optional[Path],
    protein_file: Optional[Path],
    partitions: Optional[str],
    max_memory: Optional[str],
    max_cpus: Optional[int],
    batch_size: int,
    standardized_intensities: bool,
    project_accession: Optional[str],
    compression: str,
    diann_log: Optional[Path],
    fasta: Optional[Path],
    mudata: bool,
    verbose: bool,
    metadata: _DatasetMetadata,
):
    """Convert DIA-NN report to QPX format.

    Reads a DIA-NN report.tsv file and converts feature-level quantification
    data into QPX Parquet format. When --pg-matrix-path is provided, also
    produces protein group output.

    \b
    Examples:
        # Feature conversion
        qpxc convert diann \\
            --report-path report.tsv \\
            --sdrf-file data.sdrf.tsv \\
            --mzml-info-folder ./mzml_info \\
            --output-folder ./qpx_output

        # Feature + protein groups
        qpxc convert diann \\
            --report-path report.tsv \\
            --sdrf-file data.sdrf.tsv \\
            --mzml-info-folder ./mzml_info \\
            --pg-matrix-path report.pg_matrix.tsv \\
            --output-folder ./qpx_output \\
            --standardized-intensities
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)
    prefix = output_prefix or "diann"

    from qpx.converters.diann import DiaNNConverter

    converter = DiaNNConverter(
        report_path=report_path,
        sdrf_path=sdrf_file,
        max_memory=max_memory,
        max_cpus=max_cpus,
        compression=compression,
        diann_log=str(diann_log) if diann_log else None,
    )
    converter.convert_features(
        mzml_info_folder=mzml_info_folder,
        qvalue_threshold=qvalue_threshold,
        output_folder=output_folder,
        output_prefix=output_prefix,
        protein_file=protein_file,
        partitions=partitions,
        batch_size=batch_size,
    )
    if pg_matrix_path:
        converter.convert_pg(
            pg_matrix_path=pg_matrix_path,
            output_folder=output_folder,
            output_prefix=output_prefix,
            batch_size=batch_size,
            standardized_intensities=standardized_intensities,
            qvalue_threshold=qvalue_threshold,
        )

    converter.convert_mz(
        mzml_statistics_folder=mzml_info_folder,
        output_folder=output_folder,
        output_prefix=output_prefix,
        feature_parquet=output_folder / f"{prefix}.feature.parquet",
    )
    converter.convert_sdrf(output_folder=output_folder, prefix=prefix)
    converter.write_ontology(output_folder, prefix=prefix)
    converter.write_provenance(output_folder, prefix=prefix)
    converter.write_dataset(output_folder, prefix=prefix, project_accession=project_accession)

    _finish_conversion(output_folder, prefix, project_accession, metadata, fasta=fasta, mudata=mudata)

    click.echo(f"DIA-NN conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# MaxQuant
# ---------------------------------------------------------------------------


@convert.command("maxquant")
@click.option(
    "--msms-file",
    help="MaxQuant msms.txt file (for PSM conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--evidence-file",
    help="MaxQuant evidence.txt file (for feature conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--protein-groups-file",
    help="MaxQuant proteinGroups.txt file (for PG conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file (required for feature and PG)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
)
@click.option(
    "--structures",
    help="Comma-separated list of structures to produce (psm, feature, pg). Default: all available.",
    default=None,
)
@click.option(
    "--protein-file",
    help="Protein list file for filtering feature output",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--batch-size",
    help="Processing batch size",
    default=100_000,
    type=int,
)
@click.option(
    "--max-cpus",
    help="Maximum number of threads",
    type=int,
)
@click.option(
    "--max-memory",
    help="Maximum memory limit (e.g. '16GB')",
)
@click.option(
    "--spectral-data",
    help="Include spectral data fields in PSM output",
    is_flag=True,
)
@click.option(
    "--standardized-intensities",
    help="Calculate standardized intensity metrics for PG output",
    is_flag=True,
    default=False,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD020192)",
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@click.option(
    "--fixed-mod-only",
    help="Only include features with fixed modifications (Unmodified or Carbamidomethyl C)",
    is_flag=True,
    default=False,
)
@_dataset_metadata_options
def convert_maxquant_cmd(
    msms_file: Optional[Path],
    evidence_file: Optional[Path],
    protein_groups_file: Optional[Path],
    sdrf_file: Optional[Path],
    output_folder: Path,
    output_prefix: Optional[str],
    structures: Optional[str],
    protein_file: Optional[Path],
    batch_size: int,
    max_cpus: Optional[int],
    max_memory: Optional[str],
    spectral_data: bool,
    standardized_intensities: bool,
    project_accession: Optional[str],
    compression: str,
    verbose: bool,
    fixed_mod_only: bool,
    metadata: _DatasetMetadata,
):
    """Convert MaxQuant output to QPX format.

    Reads MaxQuant result files (msms.txt, evidence.txt, proteinGroups.txt)
    and writes corresponding QPX Parquet files.

    \b
    Examples:
        # Convert everything
        qpxc convert maxquant \\
            --msms-file msms.txt \\
            --evidence-file evidence.txt \\
            --protein-groups-file proteinGroups.txt \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_output

        # Convert PSMs only
        qpxc convert maxquant \\
            --msms-file msms.txt \\
            --output-folder ./qpx_output \\
            --structures psm
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)

    # Determine which structures to produce
    if structures:
        requested = {s.strip() for s in structures.split(",")}
    else:
        # Auto-detect from provided files
        requested = set()
        if msms_file:
            requested.add("psm")
        if evidence_file:
            requested.add("feature")
        if protein_groups_file:
            requested.add("pg")
    if not requested:
        raise click.ClickException(
            "No input files provided. Supply at least one of --msms-file, --evidence-file, or --protein-groups-file."
        )

    from qpx.converters.maxquant import MaxQuantConverter

    converter = MaxQuantConverter(max_memory=max_memory, max_cpus=max_cpus, compression=compression)
    converter.convert(
        msms_file=msms_file,
        evidence_file=evidence_file,
        protein_groups_file=protein_groups_file,
        sdrf_file=sdrf_file,
        output_folder=output_folder,
        output_prefix=output_prefix,
        structures=list(requested),
        protein_file=protein_file,
        batch_size=batch_size,
        n_workers=max_cpus,
        spectral_data=spectral_data,
        standardized_intensities=standardized_intensities,
        project_accession=project_accession,
        fixed_mod_only=fixed_mod_only,
    )

    prefix = output_prefix or "maxquant"
    _finish_conversion(output_folder, prefix, project_accession, metadata)

    click.echo(f"MaxQuant conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# FragPipe
# ---------------------------------------------------------------------------


@convert.command("fragpipe")
@click.option(
    "--psm-file",
    help="FragPipe psm.tsv file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--ion-file",
    help="FragPipe combined_ion.tsv file (for feature conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--pg-file",
    help="FragPipe combined_protein.tsv file (for PG conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file (for sample/run conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--experiment-annotation-file",
    help=("FragPipe experiment_annotation.tsv mapping experiments to member raw files"),
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
)
@click.option(
    "--batch-size",
    help="Processing batch size",
    default=1_000_000,
    type=int,
)
@click.option(
    "--max-memory",
    help="Maximum memory limit (e.g. '16GB')",
)
@click.option(
    "--max-cpus",
    help="Maximum number of threads",
    type=int,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD020192)",
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@_dataset_metadata_options
def convert_fragpipe_cmd(
    psm_file: Optional[Path],
    ion_file: Optional[Path],
    pg_file: Optional[Path],
    sdrf_file: Optional[Path],
    experiment_annotation_file: Optional[Path],
    output_folder: Path,
    output_prefix: Optional[str],
    batch_size: int,
    max_memory: Optional[str],
    max_cpus: Optional[int],
    project_accession: Optional[str],
    compression: str,
    verbose: bool,
    metadata: _DatasetMetadata,
):
    """Convert FragPipe output to QPX format.

    Reads FragPipe result files and converts them into QPX Parquet format.
    Supports psm.tsv, combined_ion.tsv, and combined_protein.tsv.

    \b
    Examples:
        # Convert PSMs only
        qpxc convert fragpipe \\
            --psm-file psm.tsv \\
            --output-folder ./qpx_output

        # Convert features + protein groups
        qpxc convert fragpipe \\
            --ion-file combined_ion.tsv \\
            --pg-file combined_protein.tsv \\
            --sdrf-file metadata.sdrf.tsv \\
            --experiment-annotation-file experiment_annotation.tsv \\
            --output-folder ./qpx_output
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not any([psm_file, ion_file, pg_file]):
        raise click.ClickException("No input files provided. Supply at least one of --psm-file, --ion-file, or --pg-file.")

    output_folder.mkdir(parents=True, exist_ok=True)

    from qpx.converters.fragpipe import FragPipeConverter

    converter = FragPipeConverter(
        output_directory=output_folder, max_memory=max_memory, max_cpus=max_cpus, compression=compression
    )
    converter.convert(
        psm_file=psm_file,
        ion_file=ion_file,
        pg_file=pg_file,
        sdrf_file=sdrf_file,
        experiment_annotation_file=experiment_annotation_file,
        output_prefix=output_prefix,
        batch_size=batch_size,
        project_accession=project_accession,
    )

    prefix = output_prefix or "fragpipe"
    _finish_conversion(output_folder, prefix, project_accession, metadata)

    click.echo(f"FragPipe conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# mzIdentML
# ---------------------------------------------------------------------------


@convert.command("mzidentml")
@click.option(
    "--mzid-path",
    help="Input mzIdentML (.mzid) file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="mzidentml",
)
@click.option(
    "--mgf-path",
    help="Optional MGF file for spectra attachment",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--include-spectra",
    help="Attach mz_array and intensity_array from MGF to PSM records",
    is_flag=True,
    default=False,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD054720)",
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@_dataset_metadata_options
def convert_mzidentml_cmd(
    mzid_path: Path,
    output_folder: Path,
    output_prefix: str,
    mgf_path: Optional[Path],
    include_spectra: bool,
    project_accession: Optional[str],
    compression: str,
    verbose: bool,
    metadata: _DatasetMetadata,
):
    """Convert mzIdentML file to QPX dataset.

    Supports both standard mzIdentML (1.1/1.2) and mzIdentML 1.3 with
    cross-linking extensions (inter-peptide, looplinks, noncovalent).
    Produces a full QPX dataset including PSM, pepmap, provenance,
    ontology, and dataset metadata.

    \b
    Examples:
        # Convert a standard mzIdentML file
        qpxc convert mzidentml \\
            --mzid-path results.mzid \\
            --output-folder ./qpx_output

        # Convert an XL-MS mzIdentML 1.3 file
        qpxc convert mzidentml \\
            --mzid-path crosslinks.mzid \\
            --output-folder ./qpx_output \\
            --output-prefix xl_experiment

        # Convert with spectra from MGF
        qpxc convert mzidentml \\
            --mzid-path results.mzid \\
            --mgf-path spectra.mgf \\
            --include-spectra \\
            --output-folder ./qpx_output

        # Convert with project accession
        qpxc convert mzidentml \\
            --mzid-path results.mzid \\
            --output-folder ./qpx_output \\
            --project-accession PXD054720
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    from qpx.converters.mzidentml.converter import MzIdentMLConverter

    converter = MzIdentMLConverter(compression=compression)
    converter.convert(
        mzid_path=mzid_path,
        output_folder=output_folder,
        output_prefix=output_prefix,
        mgf_path=mgf_path,
        include_spectra=include_spectra,
        project_accession=project_accession,
    )

    _finish_conversion(output_folder, output_prefix, project_accession, metadata)

    click.echo(f"mzIdentML conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# OpenMS
# ---------------------------------------------------------------------------


@convert.command("openms", deprecated=True)
@click.option(
    "--qpx-dir",
    help="Directory containing OpenMS -out_qpx parquet files (*.psm.parquet, *.feature.parquet, *.pg.parquet)",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path (for sample/run generation)",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--consensusxml",
    "consensusxml_path",
    help="OpenMS .consensusXML (the -out_cxml companion of -out_qpx); its "
    "ColumnHeaders give the authoritative channel count/order for relabeling. "
    "Peptide identifications, when present, recover PSM runs using exact spectrum matches.",
    required=False,
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for the full QPX dataset",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="openms",
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD001819)",
)
@click.option(
    "--enrich-pride",
    help="Fetch project metadata from PRIDE API after conversion",
    is_flag=True,
    default=False,
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_openms_cmd(**kwargs):
    r"""[DEPRECATED] Enrich OpenMS ProteomicsLFQ -out_qpx output into a full QPX dataset.

    DEPRECATED: Some OpenMS exporters assign PSMs to the first run or emit
    conflicting identities (OpenMS#9872, OpenMS#9871). Prefer
    ``qpxc convert openms-consensus`` to read the original consensusXML directly.

    With --consensusxml, uniquely matching spectrum evidence restores PSM runs.
    Missing or ambiguous spectrum matches fail before existing core outputs are
    replaced. Duplicate identities are retained with warnings; use
    ``qpxc validate`` to report them as errors. Rows are not discarded
    or assigned arbitrary IDs to resolve conflicts.

    Validates the existing psm/feature/pg parquet files, copies them to the
    output folder, and generates the missing metadata tables (run, sample,
    ontology, provenance, dataset) from the SDRF file.

    \b
    Examples:
        # Enrich OpenMS QPX output
        qpxc convert openms \\
            --qpx-dir ./openms_qpx_output \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_full

        # With project accession
        qpxc convert openms \\
            --qpx-dir ./openms_qpx_output \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_full \\
            --project-accession PXD001819
    """
    qpx_dir = kwargs["qpx_dir"]
    sdrf_file = kwargs["sdrf_file"]
    output_folder = kwargs["output_folder"]
    output_prefix = kwargs["output_prefix"]
    project_accession = kwargs["project_accession"]
    enrich_pride = kwargs["enrich_pride"]
    compression = kwargs["compression"]
    verbose = kwargs["verbose"]

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    converter = OpenMSConverter(
        qpx_dir=qpx_dir,
        sdrf_path=sdrf_file,
        consensusxml_path=kwargs.get("consensusxml_path"),
        compression=compression,
    )
    try:
        converter.convert(
            output_folder=output_folder,
            output_prefix=output_prefix,
            project_accession=project_accession,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    _maybe_enrich_pride(output_folder, project_accession, enrich_pride, output_prefix)

    click.echo(f"OpenMS QPX enrichment complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# OpenMS consensusXML (interim path: consensusXML + SDRF -> QPX)
# ---------------------------------------------------------------------------


def _consensusxml_list(ctx, param, value) -> tuple[Path, ...]:
    """Split a comma-separated list of consensusXML paths and check that each file exists."""
    paths = tuple(Path(item.strip()) for item in value.split(",") if item.strip())
    if not paths:
        raise click.BadParameter("no consensusXML file given", ctx=ctx, param=param)
    for path in paths:
        if not path.is_file():
            raise click.BadParameter(f"file '{path}' does not exist", ctx=ctx, param=param)
    return paths


@convert.command("openms-consensus")
@click.option(
    "--consensusxml",
    "consensusxml_paths",
    required=True,
    callback=_consensusxml_list,
    help=(
        "OpenMS .consensusXML file (peptide feature intensities + IDs + protein inference), or a "
        "comma-separated list of files (e.g. one per sample group) to convert into one dataset; "
        "no run may appear in more than one file."
    ),
)
@click.option(
    "--sdrf-file",
    "sdrf_path",
    required=False,
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="SDRF metadata (run/sample views + grouped_runs fraction grouping). Optional.",
)
@click.option(
    "--output-folder",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Output directory for the QPX views.",
)
@click.option("--output-prefix", default="openms", show_default=True)
@click.option(
    "--structures",
    default="feature,psm,pg,run,sample",
    show_default=True,
    help="Comma-separated views to write (run/sample require --sdrf-file).",
)
@click.option(
    "--pg-top",
    default=0,
    show_default=True,
    type=int,
    help=(
        "Peptides used for the interim pg intensity (unnormalized unique-peptide sum): "
        "0 = all; 3 mirrors the quantms ProteomicsLFQ/IsobaricWorkflow default."
    ),
)
@click.option(
    "--streaming/--no-streaming",
    default=None,
    help=(
        "Force the low-memory streaming reader on/off. Default: auto (streaming for "
        "consensusXML files > 4 GB, which pyopenms would otherwise load whole into "
        "~0.8x-file RAM)."
    ),
)
@click.option("--verbose", is_flag=True, help="Enable verbose logging.")
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD001819)",
)
@click.option(
    "--include-unassigned-psms/--no-include-unassigned-psms",
    default=True,
    show_default=True,
    help=(
        "Keep unassigned peptide identifications by default to preserve identification "
        "evidence. Use --no-include-unassigned-psms to exclude them. Feature links "
        "are only written when both feature and PSM views are exported; a null "
        "feature_id means no recorded link. Protein inference uses every "
        "identification either way."
    ),
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option(
    "--fasta",
    help=(
        "Optional protein FASTA used for the search. When given, fills null pg.sequence_coverage, "
        "pg.molecular_weight and feature.pg_positions for target rows after conversion; proteins "
        "absent from it (e.g. DIA-NN internal decoys) stay null. See 'qpxc transform protein-properties'."
    ),
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--mudata/--no-mudata",
    default=True,
    show_default=True,
    help=(
        "Write the dataset's MuData (.h5mu) view after conversion. Best-effort: a view that "
        "cannot be built is reported and the conversion still succeeds. Use --no-mudata to skip "
        "the build on very large datasets."
    ),
)
@_dataset_metadata_options
def convert_openms_consensus_cmd(
    consensusxml_paths,
    sdrf_path,
    output_folder,
    output_prefix,
    structures,
    pg_top,
    streaming,
    verbose,
    project_accession,
    include_unassigned_psms,
    compression,
    fasta,
    mudata,
    metadata: _DatasetMetadata,
):
    """Convert one or more OpenMS consensusXML files (+ SDRF) to one QPX dataset.

    Interim quantms path while OpenMS -out_qpx is pre-1.1. The feature view carries
    the per-run/channel peptide intensities from the consensusXML; the pg view
    carries an interim unnormalized unique-peptide-sum protein intensity (stamped
    with a quantification_method cv_param) until OpenMS provides the authoritative
    protein quant. Several comma-separated --consensusxml inputs are written into the same views.
    """
    if verbose:
        logging.basicConfig(level=logging.INFO)
    from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter

    structs = tuple(s.strip() for s in structures.split(",") if s.strip())
    try:
        written = OpenMSConsensusConverter().convert(
            consensusxml_path=[str(path) for path in consensusxml_paths],
            output_folder=str(output_folder),
            output_prefix=output_prefix,
            include_unassigned_psms=include_unassigned_psms,
            sdrf_path=str(sdrf_path) if sdrf_path else None,
            structures=structs,
            pg_top=pg_top,
            streaming=streaming,
            project_accession=project_accession,
            compression=compression,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    _finish_conversion(
        output_folder,
        output_prefix,
        project_accession,
        metadata,
        fasta=fasta,
        mudata=mudata,
        written=bool(written),
    )

    click.echo(f"consensusXML conversion complete. Wrote: {sorted(written)}")


# ---------------------------------------------------------------------------
# SDRF
# ---------------------------------------------------------------------------


@convert.command("sdrf")
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="sdrf",
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_sdrf_cmd(
    sdrf_file: Path,
    output_folder: Path,
    output_prefix: str,
    compression: str,
    verbose: bool,
):
    """Convert SDRF metadata to QPX sample.parquet and run.parquet.

    Reads a Sample and Data Relationship Format (SDRF) file and produces
    the QPX sample and run data structures as Parquet files.

    \b
    Example:
        qpxc convert sdrf \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_output
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)

    from qpx.converters.sdrf import SdrfConverter

    with SdrfConverter(compression=compression) as converter:
        converter.convert(
            sdrf_path=str(sdrf_file),
            sample_output=str(output_folder / f"{output_prefix}.sample.parquet"),
            run_output=str(output_folder / f"{output_prefix}.run.parquet"),
        )

    click.echo(f"SDRF conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# QuantMS MSstats
# ---------------------------------------------------------------------------


@convert.command("quantms-msstats")
@click.option(
    "--msstats-file",
    help="QuantMS-generated *_msstats_in.csv file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option("--output-prefix", help="Prefix for output file names")
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD007683)",
)
@click.option(
    "--max-memory",
    help="Maximum DuckDB memory limit",
    default="16GB",
    show_default=True,
)
@click.option(
    "--max-cpus",
    help="Maximum number of DuckDB threads",
    default=4,
    show_default=True,
    type=int,
)
@click.option(
    "--batch-size",
    help="Number of Feature records written per batch",
    default=50_000,
    show_default=True,
    type=click.IntRange(min=1),
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@_dataset_metadata_options
def convert_quantms_msstats_cmd(
    msstats_file: Path,
    sdrf_file: Path,
    output_folder: Path,
    output_prefix: Optional[str],
    project_accession: Optional[str],
    max_memory: str,
    max_cpus: int,
    batch_size: int,
    compression: str,
    verbose: bool,
    metadata: _DatasetMetadata,
):
    """Convert QuantMS MSstats and SDRF inputs to QPX.

    Produces Feature, Sample, Run, Dataset, Provenance, and Ontology Parquet
    files. MSstats does not contain the evidence required to produce PSM or
    protein-group views, so those views are intentionally not generated.

    \b
    Example:
        qpxc convert quantms-msstats \
            --msstats-file PXD007683.sdrf_openms_design_msstats_in.csv \
            --sdrf-file PXD007683.sdrf.tsv \
            --output-folder ./qpx_output
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    from qpx.converters.quantms_msstats import QuantmsMsstatsConverter

    converter = QuantmsMsstatsConverter(
        max_memory=max_memory,
        max_cpus=max_cpus,
        compression=compression,
    )
    prefix = converter.convert(
        msstats_file=msstats_file,
        sdrf_file=sdrf_file,
        output_folder=output_folder,
        output_prefix=output_prefix,
        project_accession=project_accession,
        batch_size=batch_size,
    )

    _finish_conversion(output_folder, prefix, project_accession, metadata)

    click.echo(f"QuantMS MSstats conversion complete. Output prefix: {prefix}")


# ---------------------------------------------------------------------------
# Spectronaut
# ---------------------------------------------------------------------------


@convert.command("spectronaut")
@click.option(
    "--report-path",
    help="Spectronaut report file path (TSV or Parquet)",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--qvalue-threshold",
    help=(
        "Optional q-value threshold. Unset (the default) converts the report as "
        "reported — no filtering. When given, the feature view filters on the "
        "precursor q-value."
    ),
    default=None,
    type=float,
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option("--output-prefix", help="Prefix for output file names")
@click.option(
    "--max-memory",
    help="Maximum memory limit (e.g. '16GB')",
)
@click.option(
    "--max-cpus",
    help="Maximum number of threads",
    type=int,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD020192)",
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@_dataset_metadata_options
def convert_spectronaut_cmd(
    report_path: Path,
    sdrf_file: Optional[Path],
    qvalue_threshold: Optional[float],
    output_folder: Path,
    output_prefix: Optional[str],
    max_memory: Optional[str],
    max_cpus: Optional[int],
    project_accession: Optional[str],
    compression: str,
    verbose: bool,
    metadata: _DatasetMetadata,
):
    r"""Convert Spectronaut report to QPX format.

    Reads a Spectronaut report TSV file and converts feature-level and
    protein-group quantification data into QPX Parquet format.

    \b
    Examples:
        # Feature + protein groups
        qpxc convert spectronaut \
            --report-path report.tsv \
            --output-folder ./qpx_output

        # With SDRF metadata
        qpxc convert spectronaut \
            --report-path report.tsv \
            --sdrf-file data.sdrf.tsv \
            --output-folder ./qpx_output
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)
    prefix = output_prefix or "spectronaut"

    from qpx.converters.spectronaut import SpectronautConverter

    converter = SpectronautConverter(
        report_path=report_path,
        sdrf_path=sdrf_file,
        max_memory=max_memory,
        max_cpus=max_cpus,
        compression=compression,
    )
    converter.convert_features(
        qvalue_threshold=qvalue_threshold,
        output_folder=output_folder,
        output_prefix=prefix,
    )
    converter.convert_pg(
        output_folder=output_folder,
        output_prefix=prefix,
    )

    if sdrf_file:
        converter.convert_sdrf(output_folder=output_folder, prefix=prefix)
    converter.write_ontology(output_folder, prefix=prefix)
    converter.write_provenance(output_folder, prefix=prefix)
    converter.write_dataset(output_folder, prefix=prefix, project_accession=project_accession)

    _finish_conversion(output_folder, prefix, project_accession, metadata)

    click.echo(f"Spectronaut conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# CDAP
# ---------------------------------------------------------------------------


@convert.command("cdap")
@click.option(
    "--psm-dir",
    help="Directory containing CDAP *.psm files for one study",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="cdap",
    show_default=True,
)
@click.option(
    "--structures",
    help="Comma-separated list of structures to produce (psm, feature, pg). Default: all.",
    default=None,
)
@click.option(
    "--batch-size",
    help="Processing batch size",
    default=200_000,
    type=int,
    show_default=True,
)
@click.option(
    "--project-accession",
    help="PDC / ProteomeXchange accession (e.g. PDC000227)",
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option(
    "--max-memory",
    help="Maximum memory limit (e.g. '16GB')",
    default="16GB",
    show_default=True,
)
@click.option(
    "--max-cpus",
    help="Maximum number of threads",
    default=4,
    type=int,
    show_default=True,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@_dataset_metadata_options
def convert_cdap_cmd(
    psm_dir: Path,
    output_folder: Path,
    output_prefix: str,
    structures: Optional[str],
    batch_size: int,
    project_accession: Optional[str],
    compression: str,
    max_memory: str,
    max_cpus: int,
    verbose: bool,
    metadata: _DatasetMetadata,
):
    """Convert CPTAC CDAP .psm files to QPX format.

    Reads all *.psm tab-separated files in a study directory (e.g. PDC000440/)
    and writes psm.parquet, feature.parquet, and pg.parquet.

    \b
    Examples:
        # Convert a single CPTAC study
        qpxc convert cdap \\
            --psm-dir /data/CPTAC/PDC000440 \\
            --output-folder ./qpx_output/PDC000440 \\
            --project-accession PDC000440

        # Only produce PSM and Feature (skip PG)
        qpxc convert cdap \\
            --psm-dir /data/CPTAC/PDC000440 \\
            --output-folder ./qpx_output \\
            --structures psm,feature
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)

    requested = None
    if structures:
        requested = [s.strip() for s in structures.split(",")]

    from qpx.converters.cdap import CdapConverter

    converter = CdapConverter(
        max_memory=max_memory,
        max_cpus=max_cpus,
        compression=compression,
    )
    converter.convert(
        psm_dir=psm_dir,
        output_folder=output_folder,
        output_prefix=output_prefix,
        structures=requested,
        batch_size=batch_size,
        project_accession=project_accession,
    )

    _finish_conversion(output_folder, output_prefix, project_accession, metadata)

    click.echo(f"CDAP conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# mz (full spectra)
# ---------------------------------------------------------------------------


@convert.command("mz")
@click.option(
    "--mzml-dir",
    help="Directory containing mzML / mzML.gz spectra files (one study)",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output",
    "-o",
    help="Output .mz.parquet file path",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--ms-levels",
    default=None,
    help="Comma-separated MS levels to include (e.g. '2' or '1,2'). Default: all levels.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_mz_cmd(
    mzml_dir: Path,
    output: Path,
    ms_levels: Optional[str],
    verbose: bool,
):
    """Convert a directory of mzML spectra to a QPX mz.parquet (full spectra).

    Reads every .mzML / .mzML.gz file in --mzml-dir and writes scan-level
    spectral data (m/z and intensity arrays, precursor info) to a single
    mz.parquet. Each spectrum carries run_file_name + scan so it can be linked
    back to PSM / feature records.

    \b
    Examples:
        # All MS levels (full spectra; required for quantms precursor LFQ)
        qpxc convert mz \\
            --mzml-dir /data/CPTAC/PDC000109/mzml \\
            --output ./qpx_output/PDC000109.mz.parquet

        # Only MS2 spectra
        qpxc convert mz \\
            --mzml-dir /data/CPTAC/PDC000109/mzml \\
            --output ./PDC000109.mz.parquet \\
            --ms-levels 2
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    levels = [int(x.strip()) for x in ms_levels.split(",")] if ms_levels else None

    output.parent.mkdir(parents=True, exist_ok=True)

    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    mapper = SpectraMappingTransform(mzml_directory=mzml_dir)
    mapper.write_mz_parquet_from_dir(output, ms_levels=levels)

    click.echo(f"mz conversion complete. Output: {output}")
