"""SQL safety tests for MuData label detection."""

from unittest.mock import Mock

import pytest

from qpx.mudata import _is_multiplexed


@pytest.mark.parametrize(
    ("table", "label_field"),
    [
        ("feature", "channel"),
        ("feature", "label"),
        ("pg", "channel"),
        ("pg", "label"),
    ],
)
def test_is_multiplexed_uses_static_query_whitelist(table, label_field):
    """Supported identifier pairs select a fixed SQL statement."""
    engine = Mock()
    engine.execute.return_value.fetchall.return_value = [(2,)]

    assert _is_multiplexed(engine, table, label_field)
    engine.execute.assert_called_once()
    query = engine.execute.call_args.args[0]
    assert f"FROM {table}" in query
    if table == "pg":
        # pg is flattened since 1.1: a scalar ``label`` column (no intensities
        # list to UNNEST), so both label_field variants read ``label`` directly.
        assert "label" in query
    else:
        assert f"i.{label_field}" in query


@pytest.mark.parametrize(
    ("table", "label_field"),
    [
        ("feature; DROP TABLE pg", "label"),
        ("feature", "label) FROM pg; --"),
    ],
)
def test_is_multiplexed_rejects_unlisted_identifiers(table, label_field):
    """Unlisted identifiers are rejected before SQL execution."""
    engine = Mock()

    assert not _is_multiplexed(engine, table, label_field)
    engine.execute.assert_not_called()
