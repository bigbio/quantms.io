import streamlit as st

from quantms_io.core.database import Database
import os
import random
from st_aggrid import AgGrid,ColumnsAutoSizeMode
from st_aggrid.grid_options_builder import GridOptionsBuilder

st.set_page_config(page_title="quantms_analyze", layout="wide")
workdir = st.sidebar.text_input('workdir',value='')

def load_unique_list(file_type,View):
    if file_type.endswith('.parquet'):
        return View.get_unique_peptides(index)
    else:
        return View.get_unique_proteins(index)

def refresh_keys(query_list):
    random.shuffle(query_list)

if os.path.exists(workdir):
    os.chdir(workdir)

    View = Database()

    file_list = os.listdir()
    file_list = [file for file in file_list if file.endswith('.parquet') or file.endswith(".absolute.tsv")]

    parquet_file_list = [file for file in file_list if file.endswith('.parquet')]
    ibaq_list = [file for file in file_list if file.endswith(".absolute.tsv")]

    for file in file_list:
        View.load_db(file)


    file_type = st.sidebar.selectbox(
        "Select the file to load",
        file_list
    )


    index = None
    if file_type.endswith('.parquet'):
        index = parquet_file_list.index(file_type)
    else:
        index = ibaq_list.index(file_type)
    
    query_list = load_unique_list(file_type,View)
    query_key = st.sidebar.selectbox(
        "query peptide or protein",
        query_list[:10]
    )
    st.sidebar.button("Refresh", on_click=refresh_keys(query_list))

    if file_type.endswith('.parquet'):
        df = View.query_peptide(query_key, index)
    else:
        df = View.query_protein(query_key, index)

    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_pagination(paginationAutoPageSize=False,paginationPageSize=15)
    gb.configure_side_bar()
    gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, aggFunc="sum", editable=False)
    gridOptions = gb.build()
    if file_type.endswith('.parquet'):
        AgGrid(df,gridOptions,columns_auto_size_mode=ColumnsAutoSizeMode.FIT_CONTENTS,theme='alpine')
    else:
        AgGrid(df,gridOptions,columns_auto_size_mode=ColumnsAutoSizeMode.FIT_ALL_COLUMNS_TO_VIEW,theme='alpine')
    


