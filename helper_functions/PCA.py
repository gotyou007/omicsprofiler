import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import numpy as np
import plotly.express as px


class PCana():
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
    @st.cache_data

    def PC_s(_self, df_norm, top_n_genes, metadata):
        #remove duplicates 
        df_PCA = df_norm[~df_norm.index.duplicated(keep='first')]
        #most variable genes
        order = list(df_PCA.var(axis = 1).sort_values(ascending=False).index)
        #reindex with order and get the most variable n genes
        df_PCA = df_PCA.reindex(order).head(top_n_genes)
        # First center and scale the data
        scaled_data = scale(df_PCA.T)
        pca = PCA(n_components=2) # create a PCA object
        pca.fit(scaled_data) # do the math
        pca_data = pca.transform(scaled_data)
        #plot scatter on PC1 and PC2
        per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1)
        labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]
        df_PCA = pd.DataFrame(pca_data, index=df_PCA.columns, columns=labels)
        df_PCA['label'] = metadata.condition
        st.write("PCA table", df_PCA)
        st.write("label", metadata[['condition']])
        condition_number = len(metadata.condition.unique())
        color=list(np.random.choice(range(256), size=condition_number))
        cmap = dict(zip(df_PCA.label, color))
        fig2 = px.scatter(df_PCA, x="PC1", 
                          y="PC2",
                          hover_name=df_PCA.index,
                          color='label', color_discrete_map=cmap, 
                          labels={'PC1': f'PC1 (variance explained {per_var[0]})', 
                                  'PC2': f'PC2 (variance explained {per_var[1]})'},
                          )
        fig2.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            xaxis=dict(showgrid=False),
            yaxis=(dict(showgrid=False))
        )
        st.plotly_chart(fig2, use_container_width=True)
        #PCA table
        return df_PCA
    def PC_contrast(_self, metadata):
        condition_list = metadata.condition.unique()
        selected_baseline = st.selectbox(label="Select baseline condition",
                                           options=list(condition_list))
        selected_trt = st.selectbox(label="Select treated condition",
                                           options=list(condition_list))
        selected_contrast = ["condition", selected_baseline, selected_trt]
        return selected_contrast
    
PC_cluster = PCana()
                                            


    