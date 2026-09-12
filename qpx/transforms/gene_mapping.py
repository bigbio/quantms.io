"""Gene mapping transform — FASTA protein-to-gene annotation.

This transform maps protein accessions to gene names using the ``GN=`` field of
UniProt FASTA headers. The FASTA supplied by the caller is the only source: no
network service is queried, and no optional dependency is required.

The FASTA header format expected (UniProt):
    >sp|P12345|PROT_HUMAN Some description GN=BRCA1 PE=1 SV=2

The 'GN=' field is extracted as the gene name (symbol).

Usage:
    mapping = GeneMappingTransform(fasta_path="proteins.fasta")

    # Annotate a Feature or PG structure
    annotated_df = mapping.annotate_dataframe(feature_df)

    # Annotate via Dataset (returns new DataFrame, does not modify originals)
    annotated_df = mapping.annotate_dataset_features(dataset)

    # Write annotated Parquet
    mapping.write_annotated_features(dataset, "output.feature.parquet")
"""

from __future__ import annotations

import gzip
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)

_GENE_NAME_RE = re.compile(r"\bGN=(\S+)")


def _parse_fasta_header(header: str) -> tuple[str, str, Optional[str]]:
    """Split one FASTA header line into (accession, entry name, gene name).

    Handles the UniProt ``>db|ACCESSION|NAME description`` form, the two-field
    ``>db|ACCESSION`` form, and a bare ``>IDENTIFIER``. The gene name is the
    ``GN=`` field when present, else None.
    """
    header = header[1:] if header.startswith(">") else header
    header = header.strip()
    identifier = header.split(None, 1)[0] if header else ""
    parts = identifier.split("|")

    if len(parts) >= 3:
        accession, name = parts[1], parts[2]
    elif len(parts) == 2:
        accession = name = parts[1]
    else:
        accession = name = identifier

    match = _GENE_NAME_RE.search(header)
    return accession, name, match.group(1) if match else None


def _parse_gene_names_from_fasta(
    fasta_path: str,
    map_by: str = "accession",
) -> dict[str, set[str]]:
    """
    Parse FASTA file and build a protein -> gene name mapping.

    Extracts gene names from the 'GN=' field in UniProt FASTA headers. Only
    header lines are read (sequences are skipped), and ``.gz`` files are opened
    transparently.

    Args:
        fasta_path: Path to the FASTA file (optionally gzip-compressed).
        map_by: How to key the mapping:
            - "accession": use UniProt accession (e.g., P12345)
            - "name": use UniProt entry name (e.g., PROT_HUMAN)

    Returns:
        Dict mapping protein identifier to set of gene names.
    """
    gene_map: dict[str, set[str]] = defaultdict(set)

    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    with opener(fasta_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            accession, name, gene_name = _parse_fasta_header(line)

            if map_by == "accession":
                gene_map[accession].add(gene_name)
            elif map_by == "name":
                gene_map[name].add(gene_name)
            else:
                # Map by both accession and name for maximum matching
                gene_map[accession].add(gene_name)
                gene_map[name].add(gene_name)

    logger.info("Parsed gene names for %d protein identifiers from %s", len(gene_map), fasta_path)
    return gene_map


def _normalize_protein_list(accessions) -> Optional[list]:
    """Normalize a protein accessions value into a list suitable for _resolve_gene_names.

    - None / float NaN  → None
    - str               → [str]
    - dict / struct     → [item]   (list(dict) would give keys, not the struct)
    - iterable          → list(item)
    - other             → [item]
    """
    if accessions is None or (isinstance(accessions, float) and accessions != accessions):
        return None
    if isinstance(accessions, str):
        return [accessions]
    if isinstance(accessions, dict) or (
        hasattr(accessions, "dtype") and hasattr(accessions.dtype, "names") and accessions.dtype.names
    ):
        return [accessions]
    if hasattr(accessions, "__iter__"):
        return list(accessions)
    return [accessions]


def _extract_accession_from_raw(raw) -> Optional[str]:
    """Extract a bare accession string from a raw protein entry (str, dict, or struct).

    Returns None for None, NaN, or unsupported types.
    """
    if raw is None or (isinstance(raw, float) and raw != raw):
        return None
    if isinstance(raw, str):
        return raw
    if hasattr(raw, "__getitem__"):
        try:
            return raw["accession"]
        except (KeyError, TypeError):
            return None
    return None


def _resolve_gene_names(
    protein_accessions: list[str],
    gene_map: dict[str, set[str]],
) -> Optional[list[str]]:
    """
    Resolve gene names for a list of protein accessions.

    Args:
        protein_accessions: List of protein accessions from pg_accessions.
        gene_map: Protein-to-gene mapping from _parse_gene_names_from_fasta.

    Returns:
        List of gene names, or None if no mappings found.
    """
    if protein_accessions is None:
        return None

    gene_names = []
    for raw in protein_accessions:
        accession = _extract_accession_from_raw(raw)
        if accession is None:
            continue
        # Try full accession first, then extract short UniProt ID (sp|P12345|NAME → P12345)
        if accession not in gene_map:
            parts = accession.split("|")
            if len(parts) >= 2:
                accession = parts[1]
        if accession in gene_map:
            names = gene_map[accession]
            for name in names:
                if name is not None and name not in gene_names:
                    gene_names.append(name)

    return gene_names if gene_names else None


class GeneMappingTransform:
    """
    Map protein accessions to gene names from a FASTA database.

    This transform reads UniProt FASTA headers to extract gene symbols (GN= field)
    and annotates QPX Feature or PG data structures with a gg_names column. The
    FASTA is the only source of truth, so annotation is offline and reproducible.

    Usage:
        mapping = GeneMappingTransform(fasta_path="proteins.fasta")

        # Get the raw gene map
        gene_map = mapping.gene_map

        # Annotate a DataFrame (adds a gg_names column)
        annotated_df = mapping.annotate_dataframe(df, protein_col="pg_accessions")

        # Annotate a Dataset's features
        feature_df = mapping.annotate_dataset_features(dataset)

        # Write annotated features to a new Parquet file
        mapping.write_annotated_features(dataset, "output.feature.parquet")
    """

    def __init__(
        self,
        fasta_path: Union[str, Path],
        map_by: str = "accession",
    ):
        """
        Initialize the gene mapping transform.

        Args:
            fasta_path: Path to the UniProt FASTA file (optionally gzip-compressed).
            map_by: Mapping strategy ("accession" or "name").
        """
        self._fasta_path = Path(fasta_path)
        if not self._fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        self._map_by = map_by

        # Lazily computed
        self._gene_map: Optional[dict[str, set[str]]] = None

    @property
    def gene_map(self) -> dict[str, set[str]]:
        """Protein-to-gene name mapping (lazy-loaded from FASTA)."""
        if self._gene_map is None:
            self._gene_map = _parse_gene_names_from_fasta(
                str(self._fasta_path),
                map_by=self._map_by,
            )
        return self._gene_map

    def annotate_dataframe(
        self,
        df: pd.DataFrame,
        protein_col: str = "pg_accessions",
    ) -> pd.DataFrame:
        """
        Add a gg_names column to a DataFrame.

        The protein_col should contain either:
        - A list of protein accessions (e.g., from QPX pg_accessions column)
        - A single protein accession string

        A FASTA carries no genomic accessions, so ``gg_accessions`` is left
        exactly as the converter wrote it (and created as NULL when absent).

        Args:
            df: Input DataFrame with protein identifiers.
            protein_col: Column name containing protein accessions.

        Returns:
            DataFrame with an added gg_names column.
        """
        if protein_col not in df.columns:
            raise ValueError(f"Column '{protein_col}' not found in DataFrame.")

        result = df.copy()
        gene_map = self.gene_map

        # Map protein accessions to gene names
        result["gg_names"] = result[protein_col].apply(
            lambda accessions: _resolve_gene_names(
                _normalize_protein_list(accessions),
                gene_map,
            )
        )

        if "gg_accessions" not in result.columns:
            result["gg_accessions"] = None

        n_mapped = result["gg_names"].notna().sum()
        logger.info(
            "Mapped gene names for %d/%d rows (%.1f%%)",
            n_mapped,
            len(result),
            n_mapped / len(result) * 100,
        )
        return result

    def annotate_dataset_features(
        self,
        dataset,
    ) -> pd.DataFrame:
        """
        Annotate a Dataset's Feature data with gene names.

        Materializes the Feature data as a DataFrame, adds gene annotations,
        and returns the annotated DataFrame.

        Args:
            dataset: A qpx.Dataset with feature data.

        Returns:
            Annotated Feature DataFrame with gg_names.
        """
        if dataset.feature is None:
            raise ValueError("Dataset does not contain feature data.")

        feature_df = dataset.feature.to_df()
        return self.annotate_dataframe(
            feature_df,
            protein_col="pg_accessions",
        )

    def write_annotated_features(
        self,
        dataset,
        output_path: Union[str, Path],
    ) -> Path:
        """
        Write gene-annotated Feature data to a new Parquet file.

        Uses the FeatureWriter to produce a schema-validated output file.

        Args:
            dataset: A qpx.Dataset with feature data.
            output_path: Path for the output .feature.parquet file.

        Returns:
            Path to the written file.
        """
        from qpx.writers.feature import FeatureWriter

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        annotated_df = self.annotate_dataset_features(dataset)

        source_composite = dataset.feature.file_metadata.get("identity_composite")
        identity_composite = tuple(source_composite.split(",")) if source_composite else None
        with FeatureWriter(
            output_path,
            override_provided_ids=False,
            identity_composite=identity_composite,
        ) as writer:
            writer.write_dataframe(annotated_df)

        logger.info("Wrote gene-annotated features to %s", output_path)
        return output_path

    def annotate_dataset_pg(
        self,
        dataset,
    ) -> pd.DataFrame:
        """
        Annotate a Dataset's protein-group data with gene names.

        Args:
            dataset: A qpx.Dataset with pg data.

        Returns:
            Annotated PG DataFrame with gg_names.
        """
        if dataset.pg is None:
            raise ValueError("Dataset does not contain protein group data.")

        pg_df = dataset.pg.to_df()
        return self.annotate_dataframe(
            pg_df,
            protein_col="pg_accessions",
        )

    def write_annotated_pg(
        self,
        dataset,
        output_path: Union[str, Path],
    ) -> Path:
        """
        Write gene-annotated protein-group data to a new Parquet file.

        Uses the PgWriter to produce a schema-validated output file, preserving
        the source file's identity recipe so pg_id values do not change.

        Args:
            dataset: A qpx.Dataset with pg data.
            output_path: Path for the output .pg.parquet file.

        Returns:
            Path to the written file.
        """
        from qpx.writers.pg import PgWriter

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        annotated_df = self.annotate_dataset_pg(dataset)

        source_composite = dataset.pg.file_metadata.get("identity_composite")
        identity_composite = tuple(source_composite.split(",")) if source_composite else None
        with PgWriter(
            output_path,
            override_provided_ids=False,
            identity_composite=identity_composite,
        ) as writer:
            writer.write_dataframe(annotated_df)

        logger.info("Wrote gene-annotated protein groups to %s", output_path)
        return output_path
