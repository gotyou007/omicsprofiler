import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from statsmodels.stats.multitest import multipletests
from rpy2.robjects import pandas2ri

with localconverter(ro.default_converter + pandas2ri.converter):
    base = importr('base')
    stats = importr('stats')
    limma = importr('limma')
    edgeR = importr('edgeR')

class limmaProcess():
    '''
    this class will output DEG results from limma
    '''
        
    def __init__(self, counts_df, metadata, reads_cutoff, padj_cutoff, contrast):
        self.counts_df = counts_df
        self.metadata = metadata
        self.reads_cutoff = reads_cutoff
        self.padj_cutoff = padj_cutoff
        self.contrast = contrast
        
    def prepare_RO(self):
        with localconverter(ro.default_converter + pandas2ri.converter):
            self.r_data= ro.conversion.py2rpy(self.counts_df)
            self.r_design= ro.conversion.py2rpy(self.metadata.loc[:, [self.contrast[0]]])
            # Use the genes index column from data as a R String Vector
            global genes
            self.genes = ro.StrVector(
                [
                    str(index)
                    #added tovalues to convert to numpy array
                    for index in self.counts_df.index.tolist()
                    #for index in data.index.tolist()
                ]
            )
            #create a model matrix using metadata's group information
            f = base.factor(self.r_design.rx2(self.contrast[0]), 
                            levels=base.unique(self.r_design.rx2(self.contrast[0])))
            form = Formula('~0 + f') #create R formula for fitting
            form.environment['f'] = f #add the factor to the formula
            r_design = stats.model_matrix(form) #create the design matrix
            r_design.colnames = base.levels(f) #set the column names to the factor levels
    
    def run_preprecessing(self):
        d0 = edgeR.DGEList(self.r_data)
        d0 = edgeR.calcNormFactors(d0)
        y = limma.voom(d0, self.r_design)
        pandas2ri.activate()
        y = y[1]
        y = pd.DataFrame(y, index=self.counts_df.index, columns=self.counts_df.columns)
        return y

            


    def run_limma(self):
        data =  self.counts_df
        design = self.metadata.loc[:, [self.contrast[0]]]
        # Convert data and design pandas dataframes to R dataframes
        with localconverter(ro.default_converter + pandas2ri.converter):
            base = importr('base')
            stats = importr('stats')
            limma = importr('limma')
            r_data = ro.conversion.py2rpy(data)
            r_design = ro.conversion.py2rpy(design)
            # Use the genes index column from data as a R String Vector
            genes = ro.StrVector(
                [
                    str(index)
                    #added tovalues to convert to numpy array
                    for index in data.index.tolist()
                    #for index in data.index.tolist()
                ]
            )
        # Create a model matrix using design's Target column using the R formula "~0 + f" to get all the unique factors as columns
        with localconverter(ro.default_converter):
            f = base.factor(r_design.rx2(self.contrast[0]), 
                            levels=base.unique(r_design.rx2(self.contrast[0])))
            form = Formula('~0 + f') #create R formula for fitting
            form.environment['f'] = f #add the factor to the formula
            r_design = stats.model_matrix(form) #create the design matrix
            r_design.colnames = base.levels(f) #set the column names to the factor levels
            # Fit the data to the design using lmFit from limma
            fit = limma.lmFit(r_data, r_design)
            # Make a contrasts matrix with the 1st and the last unique values
            contrast_matrix = limma.makeContrasts(f"{r_design.colnames[0]}-{r_design.colnames[-1]}", levels=r_design)
            # Fit the contrasts matrix to the lmFit data & calculate the bayesian fit
            fit2 = limma.contrasts_fit(fit, contrast_matrix)
            fit2 = limma.eBayes(fit2)
            # topTreat the bayesian fit using the contrasts and add the genelist
            r_output = limma.topTreat(fit2, coef=1, genelist=genes, number=np.Inf)
            r_output = pandas2ri.rpy2py_dataframe(r_output)
        return r_output