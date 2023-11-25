import duckdb


def get_skip_rows(path):
    skip_rows = 0
    with open(path) as f:
        line = f.readline()
        while line.startswith("#"):
            skip_rows += 1
            line = f.readline()
        return skip_rows


class Database:

    def __init__(self):

        self.feature_db = []
        self.ibaq_db = []

    def load_db(self, path: str):

        if path.endswith(".parquet"):
            self.feature_db.append(duckdb.read_parquet(path))
        elif path.endswith(".absolute.tsv"):
            skip_rows = get_skip_rows(path)
            self.ibaq_db.append(duckdb.read_csv(path, sep='\t', skiprows=skip_rows))

    def get_unique_peptides(self, n):

        feature_db = self.feature_db[n]
        unique_peps = duckdb.sql(f"SELECT DISTINCT sequence FROM feature_db").df()

        return unique_peps['sequence'].tolist()

    def get_unique_proteins(self, n):

        ibaq_db = self.ibaq_db[n]
        unique_prts = duckdb.sql(f"SELECT DISTINCT protein FROM ibaq_db").df()

        return unique_prts['protein'].tolist()

    def query_peptide(self, peptide: str, n):

        feature_db = self.feature_db[n]
        return duckdb.sql(f"SELECT * FROM feature_db WHERE sequence ='{peptide}'").df()

    def query_protein(self, protein: str, n):

        ibaq_db = self.ibaq_db[n]

        return duckdb.sql(f"SELECT * FROM ibaq_db WHERE protein ='{protein}'").df()
