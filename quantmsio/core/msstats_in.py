from quantmsio.core.duckdb import DuckDB
from quantmsio.core.sdrf import SDRFHandler


class MsstatsIN(DuckDB):
    def __init__(self, report_path, sdrf_path, duckdb_max_memory="16GB", duckdb_threads=4):
        super(MsstatsIN, self).__init__(report_path, duckdb_max_memory, duckdb_threads)
        self._sdrf = SDRFHandler(sdrf_path)
        
    def iter_runs(self, file_num=10):
        references = self._sdrf.get_runs()
        ref_list = [references[i : i + file_num] for i in range(0, len(references), file_num)]
        for refs in ref_list:
            batch_df = self.query_field("Reference",refs)
            yield batch_df
    