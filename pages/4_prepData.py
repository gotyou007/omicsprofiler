import os
import streamlit as st
import pickle as pkl

import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from utils import convert_df
# Replace this with the path to directory where you would like results to be saved
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



    