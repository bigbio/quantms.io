"""consensusXML header tools become the provenance steps and the dataset software."""

import pyarrow.parquet as pq
import pytest

pytest.importorskip("pyopenms")

import qpx  # noqa: E402
from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter  # noqa: E402
from tests.converters.test_openms_consensus import _TMT_CONSENSUSXML  # noqa: E402

_FULL_HEADER_CONSENSUSXML = (
    _TMT_CONSENSUSXML.replace(
        "  <IdentificationRun ",
        """  <dataProcessing completion_time="2026-10-08T08:01:12">
    <software name="FeatureLinkerUnlabeledKD" version="3.5.0-pre-exported-20251212" />
    <processingAction name="Feature grouping" />
    <UserParam type="string" name="parameter: out" value="merged.consensusXML"/>
  </dataProcessing>
  <dataProcessing completion_time="2026-10-08T08:01:55">
    <software name="IDConflictResolver" version="3.5.0-pre-exported-20251212" />
    <processingAction name="Data filtering" />
  </dataProcessing>
  <IdentificationRun """,
    )
    .replace(
        'search_engine="" search_engine_version=""',
        'search_engine="Percolator" search_engine_version="3.07"',
    )
    .replace(
        "    </SearchParameters>",
        """      <UserParam type="string" name="CometAdapter:1:in" value="run_01.mzML"/>
      <UserParam type="string" name="SE:Comet" value="# comet_version 2024.01 rev. 1"/>
    </SearchParameters>""",
    )
)
_LABEL_FREE_CONSENSUSXML = _TMT_CONSENSUSXML.replace('label="tmt6plex_126"', 'label="label-free"').replace(
    'label="tmt6plex_127"', 'label="label-free"'
)

_STRUCTURES = ("feature", "psm")
_FEATURE_PSM = ["feature", "psm"]
_QPX_STEP = ("format_conversion", "openms_consensus_qpx_conversion", "qpx", qpx.__version__, ["ontology"])


def _without_line(xml: str, fragment: str) -> str:
    """Return the XML text with every line that contains ``fragment`` removed."""
    return "".join(line for line in xml.splitlines(keepends=True) if fragment not in line)


def _summary(steps: list[dict]) -> list[tuple]:
    """Return (category, step name, tool, version, output views) for each provenance step."""
    return [
        (step["step_category"], step["step_name"], step["tool_name"], step["tool_version"], step["output_views"])
        for step in steps
    ]


def _convert(tmp_path, paths, streaming: bool) -> tuple[list[dict], dict]:
    """Convert the consensusXML path(s) and return the provenance steps and dataset record."""
    out = tmp_path / "out"
    OpenMSConsensusConverter().convert(paths, str(out), structures=_STRUCTURES, streaming=streaming)
    provenance = pq.read_table(out / "openms.provenance.parquet").to_pylist()
    dataset = pq.read_table(out / "openms.dataset.parquet").to_pylist()[0]
    return provenance, dataset


@pytest.mark.parametrize("streaming", [False, True], ids=["memory", "streaming"])
def test_header_tools_become_provenance_steps(tmp_path, streaming):
    """Comet, Percolator and the two header processing tools are recorded in file order."""
    source = tmp_path / "input.consensusXML"
    source.write_text(_FULL_HEADER_CONSENSUSXML, encoding="utf-8")

    provenance, dataset = _convert(tmp_path, str(source), streaming)

    assert _summary(provenance) == [
        ("database_search", "sequence_search", "Comet", "# comet_version 2024.01 rev. 1", ["psm"]),
        ("rescoring", "psm_rescoring", "Percolator", "3.07", ["psm"]),
        ("quantification", "feature_grouping", "FeatureLinkerUnlabeledKD", "3.5.0-pre-exported-20251212", _FEATURE_PSM),
        ("quantification", "data_filtering", "IDConflictResolver", "3.5.0-pre-exported-20251212", _FEATURE_PSM),
        _QPX_STEP,
    ]
    assert [step["step_order"] for step in provenance] == [1, 2, 3, 4, 5]
    assert dataset["software_name"] == "OpenMS"
    assert dataset["software_version"] == "3.5.0-pre-exported-20251212"


def test_tools_of_every_input_are_listed_once_per_version(tmp_path):
    """Shared tools appear once; a tool with two versions appears once per version."""
    first = tmp_path / "first.consensusXML"
    second = tmp_path / "second.consensusXML"
    first.write_text(_FULL_HEADER_CONSENSUSXML, encoding="utf-8")
    second_xml = _FULL_HEADER_CONSENSUSXML.replace("run_01.mzML", "run_02.mzML").replace(
        'name="IDConflictResolver" version="3.5.0-pre-exported-20251212"', 'name="IDConflictResolver" version="3.6.0"'
    )
    second.write_text(second_xml, encoding="utf-8")

    provenance, dataset = _convert(tmp_path, [str(first), str(second)], streaming=False)

    names = [step["tool_name"] for step in provenance]
    assert names.count("FeatureLinkerUnlabeledKD") == 1
    assert names.count("IDConflictResolver") == 2
    assert names.count("Comet") == 1
    assert names.count("Percolator") == 1
    assert {step["tool_version"] for step in provenance if step["tool_name"] == "IDConflictResolver"} == {
        "3.5.0-pre-exported-20251212",
        "3.6.0",
    }
    assert dataset["software_name"] == "OpenMS"
    assert dataset["software_version"] is None


def test_percolator_without_original_engine_is_the_database_search(tmp_path):
    """Without an SE:<engine> parameter the identification engine is the database search."""
    source = tmp_path / "no_original.consensusXML"
    source.write_text(_without_line(_FULL_HEADER_CONSENSUSXML, "SE:Comet"), encoding="utf-8")

    provenance, _ = _convert(tmp_path, str(source), streaming=False)

    assert _summary(provenance)[0] == ("database_search", "sequence_search", "Percolator", "3.07", ["psm"])


@pytest.mark.parametrize(
    ("xml", "step_name", "software_name"),
    [
        (_TMT_CONSENSUSXML, "isobaric_quantification", "OpenMS/IsobaricWorkflow"),
        (_LABEL_FREE_CONSENSUSXML, "label_free_quantification", "OpenMS/ProteomicsLFQ"),
    ],
    ids=["isobaric", "label-free"],
)
def test_header_without_tools_falls_back_to_the_workflow(tmp_path, xml, step_name, software_name):
    """No dataProcessing and no search engine leaves only the workflow fallback step."""
    source = tmp_path / "no_tools.consensusXML"
    source.write_text(xml, encoding="utf-8")

    provenance, dataset = _convert(tmp_path, str(source), streaming=False)

    assert _summary(provenance) == [
        ("quantification", step_name, software_name, None, _FEATURE_PSM),
        _QPX_STEP,
    ]
    assert dataset["software_name"] == software_name
    assert dataset["software_version"] is None
