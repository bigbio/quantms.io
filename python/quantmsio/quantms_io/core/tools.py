'''
tools.py provide some function for generate optional feature:
mz_array: 
intensity_array:
num_peaks:
gene_names:
gene_accessions:
'''
from quantms_io.core.openms import OpenMSHandler
from quantms_io.core.feature import FeatureHandler
import pandas as pd
import pyarrow.parquet as pq
import pyarrow as pa

# optional about spectrum
def map_spectrum_mz(mz_path: str, scan: str, Mzml: OpenMSHandler, mzml_directory: str):
    """
    mz_path: mzML file path
    scan: scan number
    mzml: OpenMSHandler object
    """
    if mzml_directory.endswith('/'):
        mz_path = mzml_directory + mz_path + '.mzML'
    else:
        mz_path = mzml_directory + '/' + mz_path + '.mzML'
    mz_array, array_intensity = Mzml.get_spectrum_from_scan(mz_path, int(scan))
    return mz_array, array_intensity, 0 
 

def generate_features_of_spectrum(parquet_path: str, mzml_directory: str):
    table =  pd.read_parquet(parquet_path)
    Mzml = OpenMSHandler()

    table[['mz', 'array_intensity', 'num_peaks']] = table[
    ['best_psm_reference_file_name', 'best_psm_scan_number']].apply(
                        lambda x: map_spectrum_mz(x['best_psm_reference_file_name'], x['best_psm_scan_number'],
                                                        Mzml,mzml_directory), axis=1, result_type="expand")
    Feature = FeatureHandler()
    parquet_table = pa.Table.from_pandas(table, schema=Feature.schema)
    pq.write_table(parquet_table, parquet_path, compression="snappy", store_schema=True)

#option about gene