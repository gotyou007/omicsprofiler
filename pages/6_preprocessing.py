import streamlit as st
from helper_functions.session_state import ss
from helper_functions.PCA import PC_cluster
from helper_functions.DeSeq2 import Deseq2


import matplotlib.pyplot as plt 

st.session_state.update(st.session_state)
ss.initialise_state({"counts_cutoff": 0,
                     "comparisons": None,
                     "base_condition": None,
                     "contrast": None,
                     "norm_counts": None,
                     "res_df": None,
})

try:
    exprdict, metadatadict, Deseq2dict = st.session_state['expr_dict'], st.session_state['meta_dict'], st.session_state['Deseq2_dict']
    adjusted_dfs = {}
    prep_exp = st.sidebar.expander("Pre-processing Options")
    padj_mtds = {"DeSeq2": "DeSeq", "BH": "BH", "Bonferroni": "bonferroni", "Holm": "holm"}
    # 1. process raw counts and meta
    # 2. process DeSeq2 results
    if exprdict and metadatadict:
        metadata = metadatadict[list(metadatadict.keys())[0]]
        counts_df = exprdict[list(exprdict.keys())[0]]
        counts_df = counts_df.T
        st.write("counts", counts_df)
        st.write("metadata", metadata)

        #filter out low reads genes
        reads_cutoff = prep_exp.slider("reads number", 5, 20, 5)
        padj_cutoff = st.select_slider("padj cutoff", value=0.05, options=[0.001, 0.005, 0.01, 0.05])
        genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= reads_cutoff]
        counts_df = counts_df[genes_to_keep]
        st.write("now you have", counts_df)
        #run deseq2
        dds = Deseq2(counts_df, metadata, reads_cutoff, 0.05, st.session_state['contrast'])
        #1. run the deseq2() method to fit dispersions and LFCs and return, 
        #2. run stats and return results_df,normalized counts
        res_df, deseq2_counts = dds.run_deseq2() 
        
        ##PCA
        
        df_PCA = PC_cluster.PC_s(deseq2_counts.T, 100, metadata)
        selected_contrast = PC_cluster.PC_contrast(metadata)
        ss.save_state({'contrast':selected_contrast})
        st.write("selected contrast", st.session_state['contrast'])

       
        st.write("stat_res", res_df)
        pass
    else:
        pass
except KeyError:
    st.error("Perhaps you forgot to upload a dataset or use the demo data?")