"""PublicOntology — bionty-style CV term registry backed by Parquet + DuckDB.

Three-layer resolution:
    1. Repo checkout: ``data/ontology/`` at the repo root (developers)
    2. Auto-updated cache: ``~/.qpx/ontology/`` (pip users, 24h version check)
    3. Bundled fallback: ``qpx/core/ontology/data/`` (offline / first use)

Usage::

    from qpx.core.ontology import PublicOntology

    ms = PublicOntology("psi_ms")
    ms.search("percolator")           # fuzzy search
    ms.lookup("MS:1001492")            # exact by accession
    ms.lookup_name("percolator:score") # case-insensitive name
    ms.lookup_normalized("percolator_score")
    ms.validate(["percolator_score", "unknown"])
    ms.standardize("percolator:score") # raw → accession
    ms.scores()                        # all score terms
    ms.df()                            # full DataFrame
    ms.version                         # "4.1.235"
"""

from __future__ import annotations

import json
import logging
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import duckdb
import pandas as pd
import yaml

from qpx.core.obo import read_parquet_ontology_metadata
from qpx.core.sql import escape_path as _escape_path
from qpx.core.sql import sql_build, validate_table

# Allowed column names for search/validate to prevent SQL injection
_ALLOWED_FIELDS = frozenset(
    {
        "name",
        "normalized_name",
        "accession",
        "source",
        "definition",
    }
)

# Regex for safe identifiers (view names, source names)
_SAFE_IDENTIFIER = re.compile(r"^[a-zA-Z_][a-zA-Z0-9_]*$")


logger = logging.getLogger(__name__)

# Module-level singleton cache
_instances: dict[str, PublicOntology] = {}

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

_THIS_DIR = Path(__file__).parent
_SOURCES_PATH = _THIS_DIR / "sources.yaml"


def _load_sources_config() -> dict:
    """Load sources.yaml shipped with the package."""
    with open(_SOURCES_PATH, encoding="utf-8") as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Three-layer resolution helpers
# ---------------------------------------------------------------------------


def _find_repo_checkout(source_name: str) -> Optional[Path]:
    """Walk up from this file looking for data/ontology/{source}.parquet at repo root."""
    current = _THIS_DIR
    for _ in range(10):  # limit depth
        candidate = current / "data" / "ontology" / f"{source_name}.parquet"
        if candidate.is_file():
            return candidate
        parent = current.parent
        if parent == current:
            break
        current = parent
    return None


def _cache_dir(config: dict) -> Path:
    """Resolve the cache directory from config."""
    raw = config.get("cache_dir", "~/.qpx/ontology")
    return Path(raw).expanduser()


def _cache_meta_path(cache: Path) -> Path:
    return cache / "cache_meta.json"


def _load_cache_meta(cache: Path) -> dict:
    path = _cache_meta_path(cache)
    if path.is_file():
        try:
            return json.loads(path.read_text("utf-8"))
        except (json.JSONDecodeError, OSError):
            pass
    return {"last_version_check": None, "sources": {}}


def _save_cache_meta(cache: Path, meta: dict) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    path = _cache_meta_path(cache)
    path.write_text(json.dumps(meta, indent=2), "utf-8")


def _version_check_due(meta: dict, interval_hours: int) -> bool:
    last_check = meta.get("last_version_check")
    if last_check is None:
        return True
    try:
        last_dt = datetime.fromisoformat(last_check)
        elapsed = (datetime.now(timezone.utc) - last_dt).total_seconds()
        return elapsed > interval_hours * 3600
    except (ValueError, TypeError):
        return True


def _fetch_versions_yaml(base_url: str, timeout: int = 5) -> Optional[dict]:
    """Fetch versions.yaml from the repo (lightweight, ~200 bytes)."""
    from urllib.error import URLError

    from qpx.core.http import safe_urlopen

    url = f"{base_url}/versions.yaml"
    try:
        with safe_urlopen(url, timeout=timeout) as resp:
            text = resp.read().decode("utf-8")
        return yaml.safe_load(text)
    except (URLError, OSError, TimeoutError) as exc:
        logger.debug("Could not fetch versions.yaml: %s", exc)
        return None


def _download_file(url: str, dest: Path, timeout: int = 30) -> bool:
    """Download a file from URL to dest. Returns True on success."""
    from urllib.error import URLError

    from qpx.core.http import safe_urlopen

    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        with safe_urlopen(url, timeout=timeout) as resp:
            dest.write_bytes(resp.read())
        return True
    except (URLError, OSError, TimeoutError) as exc:
        logger.warning("Download failed %s: %s", url, exc)
        return False


# ---------------------------------------------------------------------------
# PublicOntology
# ---------------------------------------------------------------------------


class PublicOntology:
    """Bionty-style public ontology backed by Parquet + DuckDB.

    Provides search, lookup, validation, and standardization of CV terms
    from proteomics ontologies (PSI-MS, PRIDE CV, and future sources).

    Args:
        source_name: Ontology source key (e.g., ``"psi_ms"``, ``"pride_cv"``).
        auto_update: If ``True``, check for newer versions every 24h.
            Set to ``False`` for CI, offline, or test use.
    """

    def __init__(self, source_name: str, auto_update: bool = True):
        if not _SAFE_IDENTIFIER.match(source_name):
            raise ValueError(f"Invalid source_name '{source_name}': must be alphanumeric/underscore")
        self._source_name = source_name
        self._config = _load_sources_config()
        self._conn = duckdb.connect(":memory:")
        self._view_name = f"ontology_{source_name}"

        # Resolve and register Parquet
        parquet_path = self._resolve_parquet(auto_update)
        self._parquet_path = parquet_path
        self._conn.execute(
            sql_build(
                "CREATE OR REPLACE VIEW $view AS SELECT * FROM read_parquet('$path')",
                view=validate_table(self._view_name),
                path=_escape_path(str(parquet_path)),
            )
        )

        # Read version from Parquet metadata
        meta = read_parquet_ontology_metadata(parquet_path)
        self._version = meta.get("qpx_ontology_version", "unknown")

        # Populate singleton cache
        _instances[source_name] = self

    def _resolve_parquet(self, auto_update: bool) -> Path:
        """Three-layer resolution: repo checkout > cache > bundled."""
        source = self._source_name

        # Layer 1: Repo checkout
        repo_path = _find_repo_checkout(source)
        if repo_path is not None:
            logger.debug("Using repo checkout: %s", repo_path)
            return repo_path

        # Layer 2: Auto-updated cache
        cache = _cache_dir(self._config)
        cached_path = cache / f"{source}.parquet"

        if auto_update:
            meta = _load_cache_meta(cache)
            interval = self._config.get("check_interval_hours", 24)
            timeout = self._config.get("version_check_timeout_seconds", 5)
            base_url = self._config["repo_base_url"]

            if cached_path.is_file() and not _version_check_due(meta, interval):
                logger.debug("Using cached (not due for check): %s", cached_path)
                return cached_path

            # Version check is due (or no cache yet)
            remote_versions = _fetch_versions_yaml(base_url, timeout=timeout)

            if remote_versions and source in remote_versions:
                remote_ver = remote_versions[source].get("version", "")
                cached_ver = meta.get("sources", {}).get(source, {}).get("version", "")

                if remote_ver != cached_ver or not cached_path.is_file():
                    # Download new version
                    url = f"{base_url}/{source}.parquet"
                    if _download_file(url, cached_path):
                        meta["sources"][source] = {
                            "version": remote_ver,
                            "cached_at": datetime.now(timezone.utc).isoformat(),
                        }

                # Update check timestamp
                meta["last_version_check"] = datetime.now(timezone.utc).isoformat()
                _save_cache_meta(cache, meta)

            elif cached_path.is_file():
                # Network failed but we have a cached copy
                logger.warning("Could not check for updates; using cached %s", source)
                meta["last_version_check"] = datetime.now(timezone.utc).isoformat()
                _save_cache_meta(cache, meta)
            else:
                # No cache, no network — try to download anyway
                url = f"{base_url}/{source}.parquet"
                _download_file(url, cached_path)

        if cached_path.is_file():
            logger.debug("Using cached: %s", cached_path)
            return cached_path

        # Layer 3: Bundled fallback
        bundled = _THIS_DIR / "data" / f"{source}.parquet"
        if bundled.is_file():
            logger.debug("Using bundled fallback: %s", bundled)
            return bundled

        raise FileNotFoundError(
            f"No Parquet file found for ontology '{source}'. "
            f"Checked: repo checkout, cache ({cache}), bundled ({bundled}). "
            f"Run 'python -m qpx.core.ontology.build --source {source}' to generate."
        )

    # --- Properties ---

    @property
    def version(self) -> str:
        """Ontology version from Parquet metadata."""
        return self._version

    @property
    def source_name(self) -> str:
        return self._source_name

    # --- Query methods ---

    def search(self, query: str, field: str = "name", top_k: int = 10) -> pd.DataFrame:
        """Search for terms matching *query* using DuckDB ILIKE.

        Args:
            query: Search string (supports SQL wildcards).
            field: Column to search (``"name"``, ``"normalized_name"``, ``"accession"``).
            top_k: Maximum results to return.
        """
        if field not in _ALLOWED_FIELDS:
            raise ValueError(f"Invalid field '{field}': must be one of {sorted(_ALLOWED_FIELDS)}")
        stmt = sql_build(
            "SELECT * FROM $view WHERE $field ILIKE '%' || ? || '%' AND is_obsolete = false LIMIT $k",
            view=validate_table(self._view_name),
            field=field,
            k=str(int(top_k)),
        )
        return self._conn.execute(stmt, [query]).fetchdf()

    def lookup(self, accession: str) -> Optional[dict]:
        """Exact lookup by CV accession (e.g., ``"MS:1001492"``)."""
        stmt = sql_build("SELECT * FROM $view WHERE accession = ?", view=validate_table(self._view_name))
        return self._fetchone_dict(stmt, [accession])

    def lookup_name(self, name: str) -> Optional[dict]:
        """Case-insensitive lookup by original CV term name."""
        stmt = sql_build(
            "SELECT * FROM $view WHERE LOWER(name) = LOWER(?) AND is_obsolete = false",
            view=validate_table(self._view_name),
        )
        return self._fetchone_dict(stmt, [name])

    def lookup_normalized(self, name: str) -> Optional[dict]:
        """Exact lookup by normalized snake_case name."""
        stmt = sql_build(
            "SELECT * FROM $view WHERE normalized_name = ? AND is_obsolete = false",
            view=validate_table(self._view_name),
        )
        return self._fetchone_dict(stmt, [name])

    def _fetchone_dict(self, sql: str, params: list) -> Optional[dict]:
        """Execute a query and return the first row as a dict, or None."""
        result = self._conn.execute(sql, params)
        columns = [desc[0] for desc in result.description]
        rows = result.fetchall()
        if not rows:
            return None
        return dict(zip(columns, rows[0]))

    def validate(self, names: list[str], field: str = "normalized_name") -> pd.Series:
        """Check which names exist in the ontology.

        Returns a boolean Series indexed by the input names.
        """
        if field not in _ALLOWED_FIELDS:
            raise ValueError(f"Invalid field '{field}': must be one of {sorted(_ALLOWED_FIELDS)}")
        if not names:
            return pd.Series([], dtype=bool)
        placeholders = ", ".join(["?"] * len(names))
        stmt = sql_build(
            "SELECT DISTINCT $field FROM $view WHERE $field IN ($ph) AND is_obsolete = false",
            field=field,
            view=validate_table(self._view_name),
            ph=placeholders,
        )
        found = set(self._conn.execute(stmt, names).fetchdf()[field].tolist())
        return pd.Series([n in found for n in names], index=names, dtype=bool)

    def standardize(self, raw_name: str) -> Optional[str]:
        """Normalize a raw name and look up its CV accession.

        Returns the accession string, or ``None`` if not found.
        """
        from qpx.core.obo import _normalize_name

        normalized = _normalize_name(raw_name)
        term = self.lookup_normalized(normalized)
        if term is not None:
            return term["accession"]
        # Also try case-insensitive name lookup
        term = self.lookup_name(raw_name)
        if term is not None:
            return term["accession"]
        return None

    def scores(self) -> pd.DataFrame:
        """Return all terms flagged as scores."""
        stmt = sql_build(
            "SELECT * FROM $view WHERE is_score = true AND is_obsolete = false",
            view=validate_table(self._view_name),
        )
        return self._conn.execute(stmt).fetchdf()

    def df(self, include_obsolete: bool = False) -> pd.DataFrame:
        """Return the full ontology as a DataFrame."""
        vn = validate_table(self._view_name)
        if include_obsolete:
            stmt = sql_build("SELECT * FROM $view", view=vn)
        else:
            stmt = sql_build("SELECT * FROM $view WHERE is_obsolete = false", view=vn)
        return self._conn.execute(stmt).fetchdf()

    # --- Update methods ---

    def check_update(self) -> Optional[str]:
        """Check if a newer version is available in the repo.

        Returns the new version string, or ``None`` if already up to date.
        """
        base_url = self._config["repo_base_url"]
        timeout = self._config.get("version_check_timeout_seconds", 5)
        remote = _fetch_versions_yaml(base_url, timeout=timeout)
        if remote and self._source_name in remote:
            remote_ver = remote[self._source_name].get("version", "")
            if remote_ver and remote_ver != self._version:
                return remote_ver
        return None

    def update(self) -> PublicOntology:
        """Force download the latest version from the repo.

        Returns a new ``PublicOntology`` instance with the updated data.
        """
        cache = _cache_dir(self._config)
        base_url = self._config["repo_base_url"]
        url = f"{base_url}/{self._source_name}.parquet"
        dest = cache / f"{self._source_name}.parquet"

        if _download_file(url, dest):
            # Update cache meta
            meta = _load_cache_meta(cache)
            remote = _fetch_versions_yaml(base_url)
            ver = "unknown"
            if remote and self._source_name in remote:
                ver = remote[self._source_name].get("version", "unknown")
            meta["sources"][self._source_name] = {
                "version": ver,
                "cached_at": datetime.now(timezone.utc).isoformat(),
            }
            meta["last_version_check"] = datetime.now(timezone.utc).isoformat()
            _save_cache_meta(cache, meta)

        # Return fresh instance (bypass singleton)
        new_instance = PublicOntology(self._source_name, auto_update=False)
        _instances[self._source_name] = new_instance
        return new_instance

    # --- Class methods ---

    @classmethod
    def from_parquet(cls, path: str | Path) -> PublicOntology:
        """Load a PublicOntology from an arbitrary Parquet file."""
        path = Path(path)
        meta = read_parquet_ontology_metadata(path)
        raw_source = meta.get("qpx_ontology_source", path.stem)
        # Sanitize source for use as identifier
        source = re.sub(r"[^a-zA-Z0-9_]", "_", raw_source)

        instance = cls.__new__(cls)
        instance._source_name = source
        instance._config = _load_sources_config()
        instance._conn = duckdb.connect(":memory:")
        instance._view_name = f"ontology_{source}"
        instance._parquet_path = path
        instance._version = meta.get("qpx_ontology_version", "unknown")
        instance._conn.execute(
            sql_build(
                "CREATE OR REPLACE VIEW $view AS SELECT * FROM read_parquet('$path')",
                view=validate_table(instance._view_name),
                path=_escape_path(str(path)),
            )
        )
        return instance

    @classmethod
    def from_obo(cls, path: str | Path, source: str = "MS") -> PublicOntology:
        """Build a PublicOntology from a local OBO file.

        The temporary Parquet file is cleaned up when the instance is closed.
        """
        import tempfile

        from qpx.core.obo import (
            extract_obo_version,
            parse_obo,
            write_terms_parquet,
        )

        text = Path(path).read_text("utf-8")
        version = extract_obo_version(text)
        terms = parse_obo(text, source=source)

        # Write to temp dir (tracked for cleanup)
        tmpdir = tempfile.mkdtemp()
        tmp = Path(tmpdir) / f"{source}.parquet"
        write_terms_parquet(terms, tmp, source=source, version=version)

        instance = cls.from_parquet(tmp)
        instance._tmpdir = tmpdir  # track for cleanup
        return instance

    # --- Lifecycle ---

    def close(self):
        if self._conn:
            try:
                self._conn.close()
            except Exception:
                logger.debug("Failed to close DuckDB connection", exc_info=True)
            self._conn = None
        # Clean up temp dir if created by from_obo()
        tmpdir = getattr(self, "_tmpdir", None)
        if tmpdir:
            import shutil

            shutil.rmtree(tmpdir, ignore_errors=True)
            self._tmpdir = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __repr__(self):
        return f"PublicOntology('{self._source_name}', version='{self._version}', path='{self._parquet_path}')"

    def __len__(self):
        stmt = sql_build(
            "SELECT COUNT(*) FROM $view WHERE is_obsolete = false",
            view=validate_table(self._view_name),
        )
        return self._conn.execute(stmt).fetchall()[0][0]
