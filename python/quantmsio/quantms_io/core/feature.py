"""
The feature file handle and manage the feature table in column format. The main serialization format is Apache Parquet.
The feature file is defined in the docs folder of this repository.
(https://github.com/bigbio/quantms.io/blob/main/docs/FEATURE.md). Among other information, the feature file contains:
    - Peptide sequence
    - Peptide modifications
    - Peptide charge
    - Peptide retention time
    - Peptide intensity
    - Sample accession
The feature file is a column format that defines the peptide quantification/identification and its relation with each
sample in the experiment.
"""
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd

from quantms_io.core.parquet_handler import ParquetHandler

class FeatureHandler(ParquetHandler):
    """
    This class handle protein tables in column format. The main serialization format is Apache Parquet.
    """

    FEATURE_FIELDS = [pa.field("sequence", pa.string(),
                 metadata={"description": "Peptide sequence of the feature"}),
        pa.field("protein_accessions", pa.list_(pa.string()),
                 metadata={"description": "accessions of associated proteins"}),
        pa.field("protein_start_positions", pa.list_(pa.int32()),
                 metadata={"description": "start positions in the associated proteins"}),
        pa.field("protein_end_positions", pa.list_(pa.int32()),
                 metadata={"description": "end positions in the associated proteins"}),
        pa.field("unique", pa.int32(),
                 metadata={"description": "if the peptide is unique to a particular protein"}),
        pa.field("modifications", pa.list_(pa.string()),
                 metadata={"description": "peptide modifications"}),
        pa.field("retention_time", pa.float64(),
                 metadata={"description": "retention time"}),
        pa.field("charge", pa.int32(),
                 metadata={"description": "charge state of the feature"}),
        pa.field("exp_mass_to_charge", pa.float64(),
                 metadata={"description": "experimentally measured mass-to-charge ratio"}),
        pa.field("calc_mass_to_charge", pa.float64(),
                 metadata={"description": "calculated mass-to-charge ratio"}),
        pa.field("peptidoform", pa.string(),
                 metadata={"description": "peptidoform in proforma notation"}),
        pa.field("posterior_error_probability", pa.float64(),
                 metadata={"description": "posterior error probability"}),
        pa.field("global_qvalue", pa.float64(),
                 metadata={"description": "global q-value"}),
        pa.field("is_decoy", pa.int32(),
                 metadata={"description": "flag indicating if the feature is a decoy (1 is decoy, 0 is not decoy)"}),
        pa.field("best_id_score", pa.string(),
                 metadata={"description": "best identification score as key value pair"}),
        pa.field("intensity", pa.float64(),
                 metadata={"description": "intensity value"}),
        pa.field("spectral_count", pa.int32(),
                 metadata={"description": "number of spectral counts"}),
        pa.field("sample_accession", pa.string(),
                 metadata={"description": "accession of the associated sample"}),
        pa.field("condition", pa.string(),
                 metadata={"description": "experimental condition, value of the experimental factor"}),
        pa.field("fraction", pa.string(),
                 metadata={"description": "fraction information"}),
        pa.field("biological_replicate", pa.string(),
                 metadata={"description": "biological replicate information"}),
        pa.field("fragment_ion", pa.string(),
                 metadata={"description": "fragment ion information"}),
        pa.field("isotope_label_type", pa.string(),
                 metadata={"description": "type of isotope label"}),
        pa.field("run", pa.string(),
                 metadata={"description": "experimental run information"}),
        pa.field("channel", pa.string(),
                 metadata={"description": "experimental channel information"}),
        pa.field("id_scores", pa.list_(pa.string()),
                 metadata={"description": "identification scores as key value pairs"}),
        pa.field("consensus_support", pa.float64(),
                 metadata={"description": "consensus support value"}),
        pa.field("reference_file_name", pa.string(),
                 metadata={"description": "file name of the reference file"}),
        pa.field("scan_number", pa.string(),
                 metadata={"description": "scan number of the best PSM"}),
        pa.field("mz", pa.list_(pa.float64()),
                 metadata={"description": "mass-to-charge ratio values"}),
        pa.field("intensity_array", pa.list_(pa.float64()),
                 metadata={"description": "intensity array values"}),
        pa.field("num_peaks", pa.int32(),
                 metadata={"description": "number of peaks"}),
        pa.field("gene_accessions", pa.list_(pa.string()),
                 metadata={"description": "accessions of associated genes"}),
        pa.field("gene_names", pa.list_(pa.string()),
                 metadata={"description": "names of associated genes"}),
                      ]

    def __init__(self, parquet_path: str = None):
        self.schema = self._create_schema()
        self.parquet_path = parquet_path
        self.dataset = None

    def _create_schema(self):
        """
        Create the schema for the protein file. The schema is defined in the docs folder of this repository.
        (https://github.com/bigbio/quantms.io/blob/main/docs/FEATURE.md)
        """
        return pa.schema(FeatureHandler.FEATURE_FIELDS, metadata={"description": "Feature file in quantms.io format"})

    def read_feature_table(self) -> pa.Table:
        table = pq.ParquetDataset(self.parquet_path, use_legacy_dataset=False, schema=self.schema).read() # type: pa.Table
        return table

    def create_feature_table(self, feature_list: list):
        return pa.Table.from_pandas(pd.DataFrame(feature_list), schema=self.schema)

    def describe_schema(self):
        schema_description = []
        for field in self.schema:
            field_description = {
                "name": field.name,
                "type": str(field.type),
                "description": field.metadata.get("description", "")
            }
            schema_description.append(field_description)
        return schema_description

