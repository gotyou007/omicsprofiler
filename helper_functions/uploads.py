import streamlit as st
import pandas as pd
from helper_functions.session_state import ss


class FileUploads():
    def read_xfile(self, df_query, ss_excel):
        '''
        Parameter
        ---------
        df_query: from the st.file_uploader output
        ss_excel: potential session state key for any possible excel files

        Returns
        -------

        dict| keys as names of file, values as pd DataFrame
        '''

        df_dict = {}
        for d in df_query:
            head, sep, tail = str(d.name).partition(".")
            if tail == 'csv':
                data = st.cache_data(pd.read_csv)(d, index_col=0)
                df_dict[head] = data

            elif tail == 'txt':
                data = st.cache_data(pd.read_csv)(d, sep='\t', index_col=0)
                df_dict[head] = data

            elif tail == 'xlsx':
                x = st.cache_data(pd.read_excel)(d, index_col=0, sheet_name=None, engine='openpyxl')
                if ss_excel == "df_excel":
                    selected_sheet = st.multiselect(label="Select which sheet to read in", options=list(x.keys()), default = st.session_state[ss_excel])
                    if len(selected_sheet) != 0:
                        for i in selected_sheet:
                            df_dict[f"{head}_{i}"] = x[i]
                        ss.save_state({ss_excel: selected_sheet})
                    else:
                        ss.save_state({ss_excel: st.session_state[ss_excel]})
                else:
                    selected_meta = st.multiselect(label="Select which sheet to read in for metadata", options=list(x.keys()), default = st.session_state[ss_excel])
                    if len(selected_meta) != 0:
                        for i in selected_meta:
                            df_dict[f"{head}_{i}"] = x[i]
                        ss.save_state({ss_excel: selected_meta})
                    else:
                        ss.save_state({ss_excel: st.session_state[ss_excel]})
        return df_dict
    
    
    
fileuploads = FileUploads()