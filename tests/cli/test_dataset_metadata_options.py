"""Shared dataset-metadata options of the convert commands."""

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from qpx.cli.convert import _apply_software, convert
from qpx.dataset import Dataset
from tests.converters.test_openms_consensus import _TMT_CONSENSUSXML

_METADATA_OPTIONS = (
    "--software-name",
    "--software-version",
    "--software-config",
    "--enrich-pride",
    "--compute-integrity",
)
_METADATA_COMMANDS = (
    "diann",
    "maxquant",
    "fragpipe",
    "mzidentml",
    "openms-consensus",
    "quantms-msstats",
    "spectronaut",
    "cdap",
)
_DIANN_SMALL = Path(__file__).resolve().parents[1] / "examples" / "diann" / "small"


def _run(base, *extra):
    base.mkdir(parents=True, exist_ok=True)
    cx = base / "input.consensusXML"
    cx.write_text(_TMT_CONSENSUSXML)
    folder = base / "out"
    result = CliRunner().invoke(
        convert,
        [
            "openms-consensus",
            "--consensusxml",
            str(cx),
            "--output-folder",
            str(folder),
            "--output-prefix",
            "t",
            "--structures",
            "feature,pg,psm",
            *extra,
        ],
    )
    return result, folder


def test_software_options_record_workflow_in_dataset_and_provenance(tmp_path):
    params = tmp_path / "params.json"
    params.write_text(json.dumps({"b": 1, "a": [1, 2]}))
    baseline, baseline_folder = _run(tmp_path / "baseline")
    assert baseline.exit_code == 0, baseline.output
    before = pq.read_table(baseline_folder / "t.provenance.parquet").to_pylist()

    result, folder = _run(
        tmp_path,
        "--software-name",
        "nf-core/mhcquant",
        "--software-version",
        "3.1.0",
        "--software-config",
        str(params),
    )
    assert result.exit_code == 0, result.output

    dataset = pq.read_table(folder / "t.dataset.parquet").to_pylist()[0]
    assert dataset["software_name"] == "nf-core/mhcquant"
    assert dataset["software_version"] == "3.1.0"

    steps = pq.read_table(folder / "t.provenance.parquet").to_pylist()
    assert steps[0]["step_category"] == "workflow"
    assert steps[0]["step_name"] == "pipeline_execution"
    assert steps[0]["tool_name"] == "nf-core/mhcquant"
    assert steps[0]["tool_version"] == "3.1.0"
    assert steps[0]["config"] == '{"a":[1,2],"b":1}'
    assert [step["step_order"] for step in steps] == list(range(1, len(steps) + 1))
    assert steps[1]["step_name"] == before[0]["step_name"]
    assert steps[1]["tool_name"] == before[0]["tool_name"]


def test_software_version_without_name_is_usage_error(tmp_path):
    result, folder = _run(tmp_path, "--software-version", "1.0")

    assert result.exit_code == 2
    assert "require --software-name" in result.output
    assert not list(folder.glob("*.parquet"))


def test_software_config_must_be_valid_json(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")

    result, _ = _run(tmp_path, "--software-name", "x", "--software-config", str(bad))

    assert result.exit_code == 2
    assert "not a readable JSON file" in result.output


def test_compute_integrity_records_files_and_verifies(tmp_path):
    result, folder = _run(tmp_path, "--compute-integrity")

    assert result.exit_code == 0, result.output
    record = pq.read_table(folder / "t.dataset.parquet").to_pylist()[0]
    assert "t.feature.parquet" in dict(record["file_checksums"])
    assert "t.feature.parquet" in dict(record["file_row_counts"])
    assert "t.feature.parquet" in dict(record["file_sizes_bytes"])
    assert record["total_structures"] == len(dict(record["file_checksums"]))
    assert record["packaged_at"]
    assert Dataset(str(folder)).verify_integrity()["errors"] == []


def test_enrich_pride_sets_project_title(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "qpx.core.pride.fetch_pride_metadata",
        lambda acc: {"project_title": "T", "project_description": "D", "pubmed_id": "1"},
    )

    result, folder = _run(tmp_path, "--enrich-pride", "--project-accession", "PXD000001")

    assert result.exit_code == 0, result.output
    assert pq.read_table(folder / "t.dataset.parquet").to_pylist()[0]["project_title"] == "T"


def test_enrich_pride_failure_is_a_warning(tmp_path, monkeypatch):
    def _unavailable(acc):
        raise RuntimeError("PRIDE unavailable")

    monkeypatch.setattr("qpx.core.pride.fetch_pride_metadata", _unavailable)

    result, _ = _run(tmp_path, "--enrich-pride", "--project-accession", "PXD000001")

    assert result.exit_code == 0, result.output
    assert "Could not enrich from PRIDE" in result.output


def test_apply_software_without_metadata_files_warns(tmp_path, capsys):
    _apply_software(tmp_path, "t", "nf-core/mhcquant", "3.1.0", None)

    assert "WARNING: no dataset/provenance metadata found; --software-name skipped" in capsys.readouterr().err


@pytest.mark.parametrize("command", _METADATA_COMMANDS)
def test_convert_command_lists_metadata_options(command):
    result = CliRunner().invoke(convert, [command, "--help"])

    assert result.exit_code == 0, result.output
    for option in _METADATA_OPTIONS:
        assert option in result.output


@pytest.mark.parametrize("command", ["sdrf", "mz"])
def test_convert_command_without_metadata_options(command):
    result = CliRunner().invoke(convert, [command, "--help"])

    assert result.exit_code == 0, result.output
    assert "--software-name" not in result.output


def test_convert_diann_records_workflow_and_integrity(tmp_path):
    folder = tmp_path / "out"
    result = CliRunner().invoke(
        convert,
        [
            "diann",
            "--report-path",
            str(_DIANN_SMALL / "diann_report.tsv"),
            "--sdrf-file",
            str(_DIANN_SMALL / "PXD019909-DIA.sdrf.tsv"),
            "--output-folder",
            str(folder),
            "--software-name",
            "quantmsdiann",
            "--software-version",
            "2.0.0",
            "--compute-integrity",
        ],
    )

    assert result.exit_code == 0, result.output
    dataset = pq.read_table(folder / "diann.dataset.parquet").to_pylist()[0]
    assert dataset["software_name"] == "quantmsdiann"
    steps = pq.read_table(folder / "diann.provenance.parquet").to_pylist()
    assert steps[0]["step_category"] == "workflow"
    assert steps[0]["tool_name"] == "quantmsdiann"
    assert "DIA-NN" in [step["tool_name"] for step in steps]
    assert Dataset(str(folder)).verify_integrity()["errors"] == []
    assert (folder / "diann.h5mu").exists() == ("MuData view: diann.h5mu" in result.output)
