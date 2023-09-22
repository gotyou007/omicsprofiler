# Copyright (c) Streamlit Inc. (2018-2022) Snowflake Inc. (2022)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import inspect
import textwrap

import streamlit as st
import numpy as np
import pandas as pd 
from pandas.api.types import is_numeric_dtype
# from pandas_profiling import ProfileReport

def get_data_from_excel(sheet):
        df = pd.read_excel(
            io="pooled_DESeq.xlsx",
            engine="openpyxl",
            sheet_name=sheet
        )
        return df

def convert_df(df):
      return df.to_csv(index=False).encode('utf-8')
def convert_df_noindex(df):
      return df.to_csv(index=False, header=False, sep='\t').encode('utf-8')


def show_code(demo):
    """Showing the code of the demo."""
    show_code = st.sidebar.checkbox("Show code", True)
    if show_code:
        # Showing the code of the demo.
        st.markdown("## Code")
        sourcelines, _ = inspect.getsourcelines(demo)
        st.code(textwrap.dedent("".join(sourcelines[1:])))
