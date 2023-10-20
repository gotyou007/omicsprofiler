import streamlit as st
import pydeseq2
import pandas as pd
from pydeseq2.ds import DeseqStats
class Deseq2():
    '''
    This class will provide output of PCA and the contrast of interest for the user.
    Parameter
        ---------
        df_norm: either from the st.file_uploader output or the output from the DeSeq2 normalization.
        top_n_genes: most variable genes to be used for PCA
        metadata: metadata file from the st.file_uploader output

        Returns
        -------
        PCA plot and PCA table
        Contrast of interest
        
    '''
    def __init__(self, counts_df, metadata, reads_cutoff, padj_cutoff, contrast):
        self.counts_df = counts_df
        self.metadata = metadata
        self.reads_cutoff = reads_cutoff
        self.padj_cutoff = padj_cutoff
        self.contrast = contrast

    def run_norm(self):
        '''
        this function will run deseq2 normalization
        '''
        deseq2_counts, _ = pydeseq2.preprocessing.deseq2_norm(self.counts_df)
        return deseq2_counts


    def run_deseq2(self):
        #filter out nans
        samples_to_keep = ~self.metadata.condition.isna()
        self.counts_df = self.counts_df.loc[samples_to_keep]
        self.metadata = self.metadata.loc[samples_to_keep]
        #filter out low reads genes
        genes_to_keep = self.counts_df.columns[self.counts_df.sum(axis=0) >= self.reads_cutoff]
        self.counts_df = self.counts_df[genes_to_keep]
        #run deseq2
        dds = pydeseq2.dds.DeseqDataSet(
            counts=self.counts_df,
            metadata=self.metadata,
            design_factors=self.contrast[0],
            refit_cooks=True,
            n_cpus=8,
        )
        dds.deseq2()
        ## run stats (requires modifcation)
        stat_res = DeseqStats(dds, n_cpus=8, 
                              alpha=self.padj_cutoff, 
                              contrast=self.contrast)
        stat_res.summary()
        res_df = stat_res.results_df
        
        return res_df

        #########