from turtle import title
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
from sklearn import cluster
from sklearn.decomposition import PCA
from sklearn import preprocessing
from streamlit.hello.utils import get_data_from_excel, convert_df
from streamlit.hello.utils import show_code

def animation_demo() -> None:

    # Interactive Streamlit elements, like these sliders, return their value.
    # This gives you an extremely simple interaction model.
    iterations = st.sidebar.slider("Level of detail", 2, 20, 10, 1)
    separation = st.sidebar.slider("Separation", 0.7, 2.0, 0.7885)