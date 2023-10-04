import os
import streamlit as st
from helper_functions.uploads import fileuploads
from helper_functions.session_state import ss
import pandas as pd


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
                                    "meta_dict": None,
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

# #uploads conditions
# ## 1. have files uploaded and not using demo
# ## 2. no files uploaded and using demo
# ## 3. have files uploaded and using demo
# ## 4. no files uploaded and no demo
if (st.session_state['df_in'] is not None) and (not st.session_state['demo']):
    with file_opts:
        cleandict = fileuploads.read_xfile(df_query, ss_excel = "df_excel")
    #cleandict = fileuploads.capslock_genes(cleandict)

    if file_type == "DeSeq2 Results":
            #the assumption is that the user has completed preprocessing(Deseq2 or Limma)
            #and has the results of the contrast of interest(contrast of interest and p-values)
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
            meta_expr = meta_expr.data_editor(v)#testing editing function
            metadatadict[k] = meta_expr
    else:
        meta_expr.info("No metadata uploaded")

    if Deseq2dict is not None:
        for k,v in Deseq2dict.items():
            Deseq2_expr.subheader(k)
            Deseq2_expr.dataframe(v)
    else:
        Deseq2_expr.info("No ratios and p-values uploaded")




    