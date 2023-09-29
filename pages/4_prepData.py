import os
import streamlit as st
import pickle as pkl
from helper_functions.uploads import fileuploads
from helper_functions.session_state import ss
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from utils import convert_df
# ###################################################### SESSION STATES ##################################################
st.session_state.update(st.session_state)
ss.initialise_state(state_dict = { 'file_type': 'Raw Counts',
                                    'demo': False,
                                    'view_df':False,
                                    'df_in': None,
                                    'expr_dict': None,
                                    'Deseq2_dict': None,
                                    "df_excel": None,
                                    "meta_in": None,
                                    "meta_excel": None,
                                  })
# ########################################################################################################################

st.header("File Uploader")
# ################################################# File Uploader ########################################################
file_opts = st.sidebar.expander("File Upload Options", expanded = True)
use_demo = file_opts.checkbox("Use demo dataset", value=st.session_state['demo'], on_change=ss.binaryswitch, args = ('demo',))
fileType_list = ['Raw Counts', 'Metadata', 'DeSeq2 Results']
file_type = file_opts.radio(label="Select data type for your upload", options = fileType_list,
                             index = fileType_list.index(st.session_state['file_type']))
ss.save_state(dict(file_type = file_type))

df_query = file_opts.file_uploader('Upload your file', type = ['csv', 'txt', 'xlsx'], accept_multiple_files=True, help='Note that excel files take longer to upload')

clear = file_opts.button("Clear cache", on_click=ss.clear_output)
if clear:
    st.experimental_rerun()

if len(df_query) != 0:
    ss.save_state(dict(df_in = df_query))
elif st.session_state['demo']:
    ss.save_state(dict(df_in = None))
else:
    st.session_state['df_in'] = st.session_state['df_in']

# # Complicated conditions
# ## 1. have files uploaded and not using demo
# ## 2. no files uploaded and using demo
# ## 3. have files uploaded and using demo
# ## 4. no files uploaded and no demo
if (st.session_state['df_in'] is not None) and (not st.session_state['demo']):
    with file_opts:
        cleandict = fileuploads.read_xfile(df_query, ss_excel = "df_excel")
    #cleandict = fileuploads.capslock_genes(cleandict)

    if file_type == "DeSeq2 Results":
            ss.save_state({'Deseq2_dict':cleandict})
    else:
        with file_opts:
            metadata = st.file_uploader("Upload metadata", 
                                        type = ['csv', 'txt', 'xlsx'], accept_multiple_files=True)
            if len(metadata) != 0:
                ss.save_state({'expr_dict':cleandict,'meta_in':metadata})
            else:
                ss.save_state({'expr_dict':cleandict, 'meta_in':st.session_state['meta_in']})

            if st.session_state['meta_in'] is not None:
                meta_dict = fileuploads.read_xfile(st.session_state['meta_in'], ss_excel = 'meta_excel')
                ss.save_state({'meta_dict':meta_dict})
elif st.session_state['df_in'] is None and st.session_state['demo'] is True:
    databank = {
                "Raw Counts":["data/test_counts.csv", "data/test_metadata.csv"],
                }
    if file_type == "Raw Counts":
        exprmatrix = pd.read_csv(databank[st.session_state['file_type']][0], index_col=0)
        exprmeta = pd.read_csv(databank[st.session_state['file_type']][1], index_col=0)
        exprmatrix_dict = {'Demo':exprmatrix}
        exprmeta_dict = {'Demo_metadata':exprmeta}
        ss.save_state(dict(expr_dict = exprmatrix_dict,
                           meta_dict = exprmeta_dict,
                           Deseq2_dict = None))
    else:
        pass

elif st.session_state['df_in'] is not None and st.session_state['demo'] is True:
    st.warning("Deselect demo dataset to use uploaded file or clear cache to remove your file!")
else:
    ss.save_state({'Deseq2_dic':None, 'expr_dict':None, 'meta_dict':None})

view_df = file_opts.checkbox("View demo/uploaded gene expression dataset", value = st.session_state['view_df'], on_change=ss.binaryswitch, args=('view_df', ))

if view_df is True:
    main_expr, meta_expr, Deseq2_expr = st.tabs(['Gene Expression Data', 'Metadata', 'DeSeq2 Results'])
    exprdict, metadatadict, Deseq2dict = st.session_state['expr_dict'], st.session_state['meta_dict'], st.session_state['Deseq2_dict']

    if exprdict is not None:
        for k,v in exprdict.items():
            main_expr.subheader(k)
            main_expr.dataframe(v)
    else:
        main_expr.info("No gene expression counts uploaded")

    if metadatadict is not None:
        for k,v in metadatadict.items():
            meta_expr.subheader(k)
            meta_expr.dataframe(v)
    else:
        meta_expr.info("No metadata uploaded")

    if Deseq2dict is not None:
        for k,v in Deseq2dict.items():
            Deseq2_expr.subheader(k)
            Deseq2_expr.dataframe(v)
    else:
        Deseq2_expr.info("No ratios and p-values uploaded")

        






uploaded_counts = st.file_uploader("Choose a raw counts file")
if uploaded_counts is not None:
    counts_df = pd.read_csv(uploaded_counts, index_col=0)
    st.write("raw counts you put in", counts_df)
    counts_df = counts_df.T
    st.write("transposed it", counts_df)

uploaded_meta = st.file_uploader("Choose a metadata file")
if uploaded_meta is not None:
    metadata = pd.read_csv(uploaded_meta, index_col=0)
    st.write("meta data", metadata)
#check files exist
try:
    st.write("check counts", counts_df)
    st.write("check meta", metadata)
    #filter out nans
    samples_to_keep = ~metadata.condition.isna()
    counts_df = counts_df.loc[samples_to_keep]
    metadata = metadata.loc[samples_to_keep]
    #filter out low reads genes
    reads_cutoff = st.sidebar.slider("reads number", 5, 20, 5)
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= reads_cutoff]
    counts_df = counts_df[genes_to_keep]
    st.write("now you have", counts_df)

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="condition",
        refit_cooks=True,
        n_cpus=8,
    )

    dds.deseq2()
    st.write(dds)

    stat_res = DeseqStats(dds, n_cpus=8)
    stat_res.summary()
    res = convert_df(stat_res.results_df)
    st.download_button('Download contrast CSV', res, 'text/csv')

except NameError:
    print('counts or meta does not exist')    



    