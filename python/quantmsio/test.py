from quantms_io.core.feature import FeatureHandler
#from quantms_io.core.convert import FeatureConvertor
import os
import pyarrow.parquet as pq
os.chdir(r'D:\converter')

Hander = FeatureHandler('res_lfq2_no_cache.parquet')
Hander.convert_mztab_msstats_to_feature('lfq2\PXD002854-serum.sdrf_openms_design_msstats_in.csv',"lfq2\PXD002854-serum.sdrf.tsv","lfq2\PXD002854-serum.sdrf_openms_design_openms.mzTab",use_cache= False)

#Hander.convert_mztab_msstats_to_feature("MSV000079033-Blood-Plasma-TMT6.sdrf_openms_design_msstats_in.csv","MSV000079033-Blood-Plasma-TMT6.sdrf.tsv","MSV000079033-Blood-Plasma-TMT6.sdrf_openms_design_openms.mzTab",None,use_cache=True)
#Hander.convert_mztab_msstats_to_feature("lfq3\MSV000087095-DIA.sdrf_openms_design_msstats_in.csv","lfq3\MSV000087095-DIA.sdrf.tsv","lfq3\MSV000087095-DIA.sdrf_openms_design_out.mzTab",None,use_cache=True)
#Hander.convert_mztab_msstats_to_feature("ITRAQ\MSV000079033-IgY14-ITRAQ.sdrf_openms_design_msstats_in.csv","ITRAQ\MSV000079033-IgY14-ITRAQ.sdrf.tsv","ITRAQ\MSV000079033-IgY14-ITRAQ.sdrf_openms_design_openms.mzTab",None,use_cache=True)
