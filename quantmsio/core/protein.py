"""
The protein file is a file that handle and manage the protein information in quantms.io format. The protein file is defined
in the docs folder of this repository. (https://github.com/bigbio/quantms.io/blob/main/docs/PROTEIN.md). The Protein
information is a column format that defines the protein quantification/identification and its relation with each
 sample in the experiment. Among other information, the protein file contains:
 - Protein accession
 - Protein description
 - Sample accession
"""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from quantmsio.core.parquet_handler import ParquetHandler


class ProteinHandler(ParquetHandler):
    """
    This class handle protein tables in column format. The main serialization format is Apache Parquet.
    """

    PROTEIN_FIELDS = [
        pa.field(
            "protein_accessions",
            pa.list_(pa.string()),
            metadata={"description": "accessions of the protein"},
        ),
        pa.field(
            "sample_accession",
            pa.string(),
            metadata={"description": "accession of the sample in the SDRF file"},
        ),
        pa.field(
            "abundance",
            pa.float64(),
            metadata={"description": "protein abundance value selected by the workflow (e.g. MaxLFQ, iBAQ, etc)"},
        ),
        pa.field("global_qvalue", pa.float64(), metadata={"description": "global q-value"}),
        pa.field(
            "is_decoy",
            pa.int32(),
            metadata={"description": "flag indicating if the protein is a decoy (1 is decoy, 0 is not decoy)"},
        ),
        pa.field(
            "best_id_score",
            pa.string(),
            metadata={"description": "best identification score as key value pair"},
        ),
        pa.field(
            "gene_accessions",
            pa.list_(pa.string()),
            metadata={"description": "accessions of associated genes"},
        ),
        pa.field(
            "gene_names",
            pa.list_(pa.string()),
            metadata={"description": "names of associated genes"},
        ),
        pa.field(
            "number_of_peptides",
            pa.int32(),
            metadata={"description": "number of peptides associated with the protein in the given sample"},
        ),
        pa.field(
            "number_of_psms",
            pa.int32(),
            metadata={"description": "number of peptide spectrum matches in the given sample"},
        ),
        pa.field(
            "number_of_unique_peptides",
            pa.int32(),
            metadata={"description": "number of unique peptides associated with the protein"},
        ),
        pa.field(
            "protein_descriptions",
            pa.list_(pa.string()),
            metadata={"description": "protein descriptions"},
        ),
        pa.field(
            "ibaq",
            pa.float64(),
            metadata={"description": "intensity-based absolute quantification value"},
        ),
        pa.field(
            "ribaq",
            pa.float64(),
            metadata={"description": "normalized intensity-based absolute quantification value"},
        ),
        pa.field(
            "intensity",
            pa.float64(),
            metadata={"description": "sum of all peptide intensity value"},
        ),
    ]

    def __init__(self, parquet_path: str = None):
        super().__init__(parquet_path)
        self.schema = self._create_schema()
        self.parquet_path = parquet_path
        self.dataset = None

    def _create_schema(self):
        """
        Create the schema for the protein file. The schema is defined in the docs folder of this repository.
        (https://github.com/bigbio/quantms.io/blob/main/docs/PROTEIN.md)
        """
        return pa.schema(
            ProteinHandler.PROTEIN_FIELDS,
            metadata={"description": "Protein file in quantms.io format"},
        )

    def read_protein_dataset(self) -> pa.Table:
        table = pq.ParquetDataset(
            self.parquet_path, use_legacy_dataset=False, schema=self.schema
        ).read()  # type: pa.Table
        return table

    def create_proteins_table(self, protein_list: list):
        return pa.Table.from_pandas(pd.DataFrame(protein_list), schema=self.schema)

    def describe_schema(self):
        schema_description = []
        for field in self.schema:
            field_description = {
                "name": field.name,
                "type": str(field.type),
                "description": field.metadata.get("description", ""),
            }
            schema_description.append(field_description)
        return schema_description
