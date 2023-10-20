import streamlit as st
import numpy as np
import pandas as pd
import seaborn as sns
import os
import plotly.express as px  # pip install plotly-express
import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf
from matplotlib.colors import ListedColormap
from scipy.stats import zscore
from utils import get_data_from_excel
from utils import show_code

def volcano() -> None:

    # Interactive Streamlit elements, like these sliders, return their value.
    # This gives you an extremely simple interaction model.
    sheet = st.text_input("select sheet", 'comb_vs_scr')
    data = get_data_from_excel(sheet)
    st.sidebar.header("Please Filter Here:")
    pvalue = st.sidebar.selectbox('P.adj',(0.001,0.01,0.05))
    lfc = st.sidebar.selectbox('lfc',(0.5, 1, 2))
    label=("up", "not significant", "down")
    data.loc[(data["log2FoldChange"] >= lfc) & (data['padj'] < pvalue), 'label'] = label[0]  # upregulated
    data.loc[(data["log2FoldChange"] <= -lfc) & (data['padj'] < pvalue), 'label'] = label[2]  # downregulated
    data['label'].fillna(label[1], inplace=True)  # intermediate
    #export gene list
    dif_genes = data.loc[data['label'] != label[1]]
    data['-logpadj'] = -np.log10(data['padj'])
    # ---- MAINPAGE ----
    st.title("Volcano Plot")
    st.markdown("####")
    
    color=("red", "grey", "green")
    cmap = dict(zip(label, color))
    fig = px.scatter(data, x="log2FoldChange", y="-logpadj",
                        width=1200, height=800,
                        hover_name='gene_name',
                        color='label', color_discrete_map=cmap)
                        
    fig.update_layout(
        xaxis=dict(tickmode="linear"),
        plot_bgcolor="rgba(0,0,0,0)",
        yaxis=(dict(showgrid=False)),
        xaxis_range=[-4,4],
    )
    @st.cache_data 
    def convert_df(df):
        return df.to_csv().encode('utf-8')
    csv = convert_df(dif_genes)

    st.download_button(
    "Press to Download differential expressed genes",
    csv,
    "file.csv",
    "text/csv",
    key='download-csv'
    )

    st.plotly_chart(fig, use_container_width=True)

    st.markdown("""---""")

st.set_page_config(page_title="RNASeq Volcano", page_icon="📹")
st.markdown("# Volcano plot")
st.sidebar.header("RNASeq Volcano")
st.write(
    """This app shows how you can use Streamlit to build cool animations.
It displays an animated fractal based on the the Julia Set. Use the slider
to tune different parameters."""
)

volcano()