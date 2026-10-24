"""Shared read-back + assertions for the feature<->pg softlink round-trip checks.

The feature->pg association is a *computed softlink* (bigbio/qpx#269): converters
no longer materialize ``feature.pg_ids``; ``Dataset.link_feature_pg()`` derives it
on read via a label-aware join of the registered ``feature`` and ``pg`` views.
These helpers open a converted dataset, run the softlink, and assert every
``(feature_id, pg_id, label)`` edge is a real pg row whose canonical
``pg_accessions`` membership equals the feature's, whose ``grouped_runs`` contains
the feature's ``run_file_name``, and whose ``label`` is one of the labels the
feature actually carries (so a feature never over-links to a same-membership pg
row for a channel it lacks).
"""

from __future__ import annotations

import qpx


def open_converted(output_dir, prefix=None):
    """Open the converted dataset (feature + pg views) for softlink checks."""
    return qpx.Dataset(str(output_dir), file_prefix=prefix, structures=["feature", "pg"])


def feature_rows(ds):
    """``(feature_id, sequence, run_file_name, canonical membership, [labels])``."""
    return ds.sql(
        "SELECT feature_id, sequence, run_file_name, "
        "list_sort(list_transform(pg_accessions, lambda x: x.accession)) AS membership, "
        "list_transform(intensities, lambda i: i.label) AS labels "
        "FROM feature ORDER BY feature_id"
    ).fetchall()


def pg_rows(ds):
    """``(pg_id, canonical membership, [grouped_runs], label)``."""
    return ds.sql("SELECT pg_id, list_sort(pg_accessions) AS membership, grouped_runs, label FROM pg").fetchall()


def link_rows(ds):
    """``(feature_id, pg_id, label)`` from the computed softlink."""
    return ds.link_feature_pg().fetchall()


def assert_softlink_valid(ds):
    """Assert every softlink edge is a real, label-consistent pg row.

    For each ``(feature_id, pg_id, label)`` the softlink produces: the pg row
    exists, its canonical ``pg_accessions`` membership equals the feature's, the
    feature's run is in the pg row's ``grouped_runs``, and — the #269 guarantee —
    the pg row's ``label`` is one the feature actually carries in its
    ``intensities``. Returns ``(feat, pg, link)`` for further per-test checks.
    """
    feat = feature_rows(ds)
    pg = pg_rows(ds)
    link = link_rows(ds)

    feat_by_id = {fid: (set(memb or []), run, set(labels or [])) for fid, _seq, run, memb, labels in feat}
    pg_by_id = {pid: (set(memb or []), list(runs or []), label) for pid, memb, runs, label in pg}

    for fid, pid, label in link:
        assert fid in feat_by_id, f"softlink references unknown feature_id {fid}"
        assert pid in pg_by_id, f"softlink references unknown pg_id {pid}"
        f_memb, f_run, f_labels = feat_by_id[fid]
        p_memb, p_runs, p_label = pg_by_id[pid]
        assert f_memb == p_memb, f"pg {pid} membership {sorted(p_memb)} != feature {fid} {sorted(f_memb)}"
        assert f_run in p_runs, f"feature run {f_run!r} not in pg {pid} grouped_runs {p_runs}"
        assert p_label == label, f"softlink label {label!r} != pg row label {p_label!r}"
        assert label in f_labels, f"feature {fid} lacks linked label {label!r} (carries {sorted(f_labels)})"
    return feat, pg, link


def pg_ids_by_sequence(feat, link):
    """``sequence -> sorted list of linked pg_ids`` (unlinked features -> ``[]``).

    Assumes each sequence maps to a single feature (true for the synthetic
    fixtures); a feature with no softlink edge yields an empty list.
    """
    seq_by_id = {fid: seq for fid, seq, *_ in feat}
    result: dict[str, set] = {seq: set() for _fid, seq, *_ in feat}
    for fid, pid, _label in link:
        result[seq_by_id[fid]].add(pid)
    return {seq: sorted(ids) for seq, ids in result.items()}
