import streamlit as st
from helper_functions.session_state import ss
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt 

st.session_state.update(st.session_state)
ss.initialise_state({"counts_cutoff": 0,
                     "comparisons": None,
                     "base_condition": None,
                     "contrast": None,
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

        #filter out nans
        samples_to_keep = ~metadata.condition.isna()
        counts_df = counts_df.loc[samples_to_keep]
        metadata = metadata.loc[samples_to_keep]
        #filter out low reads genes
        reads_cutoff = prep_exp.slider("reads number", 5, 20, 5)
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
        dds.fit_size_factors()
        dds.deseq2()
        st.write(dds.varm["LFC"])
        stat_res = DeseqStats(dds, n_cpus=8)
        stat_res.summary(lfc_null=0.1, alt_hypothesis="greaterAbs")
        fig, ax = plt.subplots()
        fig = stat_res.plot_MA(s=20)
        st.pyplot(fig)
        pass
    else:
        pass
except KeyError:
    st.error("Perhaps you forgot to upload a dataset or use the demo data?")