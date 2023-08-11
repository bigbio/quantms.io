"""
The protein file is a file that handle and manage the protein information in quantms.io format. The protein file is defined
in the docs folder of this repository. (https://github.com/bigbio/quantms.io/blob/main/docs/PROTEIN.md). The Protein
information is a column format that defines the protein quantification/identification and its relation with each
specific sample in the experiment. Among other information, the protein file contains:
 - Protein accession
 - Protein description
 - Sample accession
"""
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd


class ProteinHandler:
    """
    This class handle protein tables in column format. The main serialization format is Apache Parquet.
    """
    def __init__(self, parquet_path: str = None):
        self.schema = self._create_schema()
        self.parquet_path = parquet_path

    def _create_schema(self):
        """
        Create the schema for the protein file. The schema is defined in the docs folder of this repository.
        (https://github.com/bigbio/quantms.io/blob/main/docs/PROTEIN.md)
        """
        fields = [
            pa.field("protein_accessions", pa.list_(pa.string()),
                     metadata={"description": "accessions of the protein"}),
            pa.field("sample_accession", pa.string(),
                     metadata={"description": "accession of the sample in the SDRF file"}),
            pa.field("abundance", pa.float64(),
                     metadata={"description": "protein abudance value selected by the workflow (e.g. MaxLFQ, iBAQ, etc)"}),
            pa.field("global_qvalue", pa.float64(),
                     metadata={"description": "global q-value"}),
            pa.field("is_decoy", pa.int32(),
                     metadata={"description": "flag indicating if the protein is a decoy (1 is decoy, 0 is not decoy)"}),
            pa.field("best_id_score", pa.string(),
                     metadata={"description": "best identification score as key value pair"}),
            pa.field("gene_accessions", pa.list_(pa.string()),
                     metadata={"description": "accessions of associated genes"}),
            pa.field("gene_names", pa.list_(pa.string()),
                     metadata={"description": "names of associated genes"}),
            pa.field("number_of_peptides", pa.int32(),
                     metadata={"description": "number of peptides associated with the protein in the given sample"}),
            pa.field("number_of_psms", pa.int32(),
                     metadata={"description": "number of peptide spectrum matches in the given sample"}),
            pa.field("number_of_unique_peptides", pa.int32(),
                     metadata={"description": "number of unique peptides associated with the protein"}),
            pa.field("protein_descriptions", pa.list_(pa.string()),
                     metadata={"description": "protein descriptions"}),
            pa.field("ibaq", pa.float64(),
                     metadata={"description": "intensity-based absolute quantification value"}),
            pa.field("ribaq", pa.float64(),
                     metadata={"description": "normalized intensity-based absolute quantification value"}),
            pa.field("intensity", pa.float64(),
                     metadata={"description": "sum of all peptide intensity value"}),
        ]
        return pa.schema(fields)

    def read_proteins(self):
        table = pq.read_table(self.parquet_path)
        return table.to_pandas()

    def create_protein(self, protein_data):
        table = pa.Table.from_pandas(pd.DataFrame([protein_data]), schema=self.schema)
        writer = pq.ParquetWriter(self.parquet_path, self.schema)
        writer.write_table(table)
        writer.close()

    def append_protein(self, protein_data):
        table = pa.Table.from_pandas(pd.DataFrame([protein_data]))
        writer = pq.ParquetWriter(self.parquet_path, self.schema, append=True)
        writer.write_table(table)
        writer.close()

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