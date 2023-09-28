import json
import random

import fastavro
import fastparquet
import pandas as pd
from faker import Faker
from io import BytesIO

import pyarrow.parquet as pq
import pyarrow as pa

# Load Avro schema from a file
#avro_schema = json.read_json("docs/psm.avsc")

# Number of records you want to generate
num_records = 100

# Create a Faker object
faker = Faker()

# Generate mock Avro data
avro_data = []
for _ in range(num_records):
    accessions = [faker.text(max_nb_chars=6), faker.text(max_nb_chars=6)]
    modifications = ["1|5-Oxidation", "6(localization score:0.98)-Phospho"]

    avro_data.append({
            "sequence": faker.text(max_nb_chars=30),
            "psm_id": faker.random_int(),
            "accessions": accessions,
            "unique": faker.random_int(min=0, max=1),
            "modifications": modifications,
            "retention_time": random.uniform(20.0, 60.0),
            "charge": faker.random_int(min=1, max=3),
            "exp_mass_to_charge": random.uniform(500.0, 600.0),
            "calc_mass_to_charge": random.uniform(490.0, 590.0),
            "psm_usi": faker.uuid4(),
            "posterior_error_probability": random.uniform(0.0, 0.1),
            "global_qvalue": random.uniform(0.0, 0.1),
            "isdecoy": faker.random_int(0, 1),
            "consensus_support": random.uniform(0.5, 1.0),
            "id_score_MS-GF:RawScore": random.uniform(0.8, 1.0),
            "id_score_MS-GF:DeNovoScore": random.uniform(0.8, 1.0),
            "mz": [random.uniform(100.0, 1000.0) for _ in range(3)],
            "intensity": [random.uniform(50.0, 100.0) for _ in range(3)],
            "num_peaks": faker.random_int(min=1, max=3),
    })

# Convert Avro data to Parquet format
# parquet_schema = fastavro.parse_schema(avro_schema)
# parquet_data = BytesIO()
# fastparquet.write(parquet_data, avro_data, parquet_schema)
#
# # Save the Parquet data to a file (replace 'output.parquet' with your desired filename)
# with open('output.parquet', 'wb') as f:
#     f.write(parquet_data.getvalue())

# df = pd.json_normalize(avro_data, sep="_")
# table = pa.Table.from_pandas(df)
# pq.write_table(table, "dump_psm.parquet")  # save json/table as parquet
#

df = pd.DataFrame(avro_data)

# Convert the DataFrame to a PyArrow Table
table = pa.Table.from_pandas(df)

# Write the PyArrow Table to a Parquet file
parquet_file_path = "dump_psm.parquet"
pq.write_table(table, parquet_file_path)
df = pd.read_parquet('dump_psm.parquet', engine='pyarrow')
print(df)