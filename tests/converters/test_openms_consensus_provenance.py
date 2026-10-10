"""consensusXML header tools become the provenance steps and the dataset software."""

import re

import pyarrow.parquet as pq
import pytest

pytest.importorskip("pyopenms")

import qpx  # noqa: E402
from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter  # noqa: E402

_CONSENSUSXML = """<?xml version="1.0" encoding="ISO-8859-1"?>
<consensusXML version="1.7" experiment_type="label-free"
  xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dataProcessing completion_time="2026-10-08T08:01:12">
    <software name="FeatureLinkerUnlabeledKD" version="3.5.0-pre-exported-20251212" />
    <processingAction name="Feature grouping" />
    <UserParam type="string" name="parameter: out" value="merged.consensusXML"/>
  </dataProcessing>
  <dataProcessing completion_time="2026-10-08T08:01:55">
    <software name="IDConflictResolver" version="3.5.0-pre-exported-20251212" />
    <processingAction name="Data filtering" />
  </dataProcessing>
  <IdentificationRun id="PI_0" date="0000-00-00T00:00:00" search_engine="Percolator" search_engine_version="3.07">
    <SearchParameters db="db.fasta" db_version="" taxonomy="" mass_type="monoisotopic" charges="2:3"
      enzyme="unspecific cleavage" missed_cleavages="0" precursor_peak_tolerance="0"
      precursor_peak_tolerance_ppm="false" peak_mass_tolerance="0" peak_mass_tolerance_ppm="false">
      <UserParam type="string" name="CometAdapter:1:in" value="run_01.mzML"/>
      <UserParam type="string" name="SE:Comet" value="# comet_version 2024.01 rev. 1"/>
    </SearchParameters>
    <ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0">
      <ProteinHit id="PH_0" accession="P12345" score="0" sequence="">
      </ProteinHit>
    </ProteinIdentification>
  </IdentificationRun>
  <mapList count="1">
    <map id="0" name="run_01.mzML" unique_id="1" label="" size="1">
    </map>
  </mapList>
  <consensusElementList>
    <consensusElement id="e_0" quality="0.0" charge="2">
      <centroid rt="100.123456" mz="450.251234" it="0.0"/>
      <groupedElementList>
        <element map="0" id="0" rt="100.123456" mz="450.251234" it="1000.0"/>
      </groupedElementList>
      <PeptideIdentification identification_run_ref="PI_0" score_type=""
        higher_score_better="true" significance_threshold="0" MZ="450.26" RT="100"
        spectrum_reference="controllerType=0 controllerNumber=1 scan=42">
        <PeptideHit score="0" sequence="PEPTIDEK" charge="2" protein_refs="PH_0">
          <UserParam type="string" name="target_decoy" value="target"/>
          <UserParam type="float" name="Posterior Error Probability_score" value="1.0e-03"/>
        </PeptideHit>
      </PeptideIdentification>
    </consensusElement>
  </consensusElementList>
</consensusXML>
"""

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
    source.write_text(_CONSENSUSXML, encoding="utf-8")

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
    first.write_text(_CONSENSUSXML, encoding="utf-8")
    second_xml = _CONSENSUSXML.replace("run_01.mzML", "run_02.mzML").replace(
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
    source.write_text(_without_line(_CONSENSUSXML, "SE:Comet"), encoding="utf-8")

    provenance, _ = _convert(tmp_path, str(source), streaming=False)

    assert _summary(provenance)[0] == ("database_search", "sequence_search", "Percolator", "3.07", ["psm"])


def test_header_without_tools_falls_back_to_the_label_free_workflow(tmp_path):
    """No dataProcessing and no search engine leaves only the label-free fallback step."""
    xml = re.sub(r"  <dataProcessing .*?</dataProcessing>\n", "", _CONSENSUSXML, flags=re.DOTALL)
    xml = xml.replace('search_engine="Percolator" search_engine_version="3.07"', 'search_engine="" search_engine_version=""')
    source = tmp_path / "no_tools.consensusXML"
    source.write_text(_without_line(xml, "SE:Comet"), encoding="utf-8")

    provenance, dataset = _convert(tmp_path, str(source), streaming=False)

    assert _summary(provenance) == [
        ("quantification", "label_free_quantification", "OpenMS/ProteomicsLFQ", None, _FEATURE_PSM),
        _QPX_STEP,
    ]
    assert dataset["software_name"] == "OpenMS/ProteomicsLFQ"
    assert dataset["software_version"] is None
