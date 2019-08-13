#!/usr/bin/env python
# coding: utf-8

# In[1]:


#packages needed to be imported
import numpy as np 
import math
import scipy.stats
import operator
import pandas as pd
import seaborn


# In[2]:


get_ipython().run_cell_magic('bash', '', 'cp /stor/work/Bio321G_Spring2019/megan_mini/dataMeanFiltered.txt .')


# In[3]:


#this prints out the gene expressed data
data = "dataMeanFiltered.txt"
geneExpressed = pd.read_csv(data)
geneExpressed.head()


# In[4]:


geneExpressed.set_index('Gene_Symbol', inplace=True)


# In[5]:


geneExpressed = geneExpressed.drop(['mean', 'variance'], axis =1)


# In[6]:


geneExpressed = geneExpressed.replace(to_replace=0,value=0.01)
geneExpressed.head()


# In[7]:


annotFile = "anot.txt"
annot = pd.read_table(annotFile, sep="\t", header=0, index_col=0)
annot.head()


# In[8]:


## pull tissue type info from annot DataFrame
tissue = annot["tissue_type"]
## simplify tissue names
tissue = tissue.str.replace("^tissue: ", "")


# In[9]:


## pull tissue type info from annot DataFrame
sample = annot["sample_type"]
## simplify tissue names
sample = sample.str.replace("^sample: ", "")


# In[10]:


#this is responsible for separating the midbrain neurons from the other samples
#will be used later to separate the PD and control samples
midbrain =  (tissue == "midbrain neuron")
midbrain = annot[midbrain]


# In[11]:


#not logged dataframe 
MB_H = midbrain.index
MB_HFrame = geneExpressed[MB_H]
MB_HFrame


# In[12]:


#logged dataframe for the midbrain 
mb = MB_HFrame.apply(np.log10)
mb


# In[13]:


#this is responsible for getting the number of samples of midbrain neurons and sporadic PD
pd_mbSamples =  ((sample == "PD")) & (tissue == "midbrain neuron")
numPD_MB = annot[pd_mbSamples].shape[0]
numPD_MB


# In[14]:


MB_PD =  annot[pd_mbSamples]
MBP = MB_PD.index
MB_LIST = MBP.tolist()
MB_HFrame[MB_LIST]


# In[15]:


#this is responsible for getting the number of samples of midbrain neurons and healthy (control)
h_mbSamples =  ((sample == "control")) & (tissue == "midbrain neuron")
numH_MB = annot[h_mbSamples].shape[0]
numH_MB


# In[16]:


MB_C =  annot[h_mbSamples]
MBC = MB_C.index
MB_C_LIST = MBC.tolist()
MB_HFrame[MB_C_LIST]


# In[17]:


#this gets the columns for midbrain; helps get the condition 
mb.columns


# In[18]:


#this is the list of conditions 1 represents with PD and -1 represents as a control
condition = [1,-1,1,-1,1,-1,-1,1,-1,1,-1,-1,1,1,-1]


# In[19]:


#t test to get corr coeff , t stat and p value
#this is responsble for the t-test for the midbrain 
def mytstatistic(x, y = condition, n = numPD_MB+numH_MB):
    corrcoeff = np.corrcoef(x.values, y)[0,1]
    tstat = ((np.sqrt(n-2))/(np.sqrt(1-np.square(corrcoeff))))*corrcoeff
    p = 2*(1-scipy.stats.t.cdf(np.abs(tstat), n-2))
    return (p)


# In[20]:


p = mytstatistic(mb,condition,numPD_MB+numH_MB)


# In[ ]:





# In[21]:


#the ones that are less than 0.05 are significantly differentially expressed
mb_Ttest = mb.apply(mytstatistic, axis=1)
mb_Ttest


# In[22]:


#this is the dictionary for midbrain 
p_dict_mb = mb_Ttest.to_dict()
p_dict_mb


# In[23]:


#this is the differtially expressed genes list for the midbrain
deg_list_mb = []
for key in p_dict_mb.keys():
    if p_dict_mb[key] < 0.01:
        deg_list_mb.append(key)


# In[24]:


#function that calls the differentialy expressed genes list for the midbrain 
deg_list_mb


# In[25]:


#this is responsible for separating the stem cells from the other samples
stem_cell =  (tissue == "stem cell")
stem_cell = annot[stem_cell]


# In[26]:


#Not logged data
SC_H = stem_cell.index
SC_HFrame = geneExpressed[SC_H]


# In[27]:


#logged data
sc = SC_HFrame.apply(np.log10)
sc


# In[28]:


#this is responsible for getting the number of samples of stem cells and sporadic PD
pd_scSamples =  ((sample == "PD")) & (tissue == "stem cell")
numPD_SC = annot[pd_scSamples].shape[0]
numPD_SC


# In[29]:


S_PD =  annot[pd_scSamples]
PDSC = S_PD.index
SD_LIST = PDSC.tolist()
SC_HFrame[SD_LIST]


# In[30]:


#this is responsible for getting the number of samples of stem cells and healthy (control)
h_scSamples =  ((sample == "control")) & (tissue == "stem cell")
numH_SC = annot[h_scSamples].shape[0]
numH_SC


# In[31]:


S_H =  annot[h_scSamples]
SCC = S_H.index
SCC_LIST = SCC.tolist()
SC_HFrame[SCC_LIST]


# In[32]:


#this is the list of conditions 1 represents with PD and -1 represents as a control
condition = [1, -1, 1, 1, -1, 1, 1, 1, -1, 1, 1,1,-1, -1, 1, -1,1, -1, 1, 1, -1, -1, -1, -1,  -1, -1, -1, -1, 1, 1, 1]


# In[33]:


#t test to get corr coeff , t stat and p value
#this is the t test for the stem cells
def mytstatistic(x, y = condition, n = numH_SC+numPD_SC):
    corrcoeff = np.corrcoef(x.values, y)[0,1]
    tstat = ((np.sqrt(n-2))/(np.sqrt(1-np.square(corrcoeff))))*corrcoeff
    p = 2*(1-scipy.stats.t.cdf(np.abs(tstat), n-2))
    return (p)


# In[34]:


#columns that help get the condition 
sc.columns


# In[35]:


p = mytstatistic(sc,condition,numH_SC+numPD_SC)


# In[36]:


#the stem cells t test
sc_Ttest = sc.apply(mytstatistic, axis=1)
sc_Ttest


# In[37]:


#pvalue dictionary of the stem cells
p_dict_sc = sc_Ttest.to_dict()
p_dict_sc


# In[38]:


#the differentially expressed genes list for stem cells
deg_list_sc = []
for key in p_dict_sc.keys():
    if p_dict_sc[key] < 0.01:
        deg_list_sc.append(key)


# In[39]:


#prints out the list for stem cells
len(deg_list_sc)


# In[40]:


#this is responsible for separating the fibroblasts from the other samples
#will be used later to separate the PD and control samples
fibroblasts =  (tissue == "Fibroblasts")
fibroblasts = annot[fibroblasts]


# In[41]:


#not logged function of fibroblasts cells
F_H = fibroblasts.index
F_HFrame = geneExpressed[F_H]


# In[42]:


#fibroblast cells logged data
F_cells = MB_HFrame.apply(np.log10)
F_cells


# In[43]:


## pull tissue type info from annot DataFrame
sample = annot["sample_type"]
## simplify tissue names
sample = sample.str.replace("^sample: ", "")


# In[44]:


#this is responsible for getting the number of samples of fibroblasts and sporadic PD
pd_fSamples =  ((sample == "PD")) & (tissue == "Fibroblasts")
numPD_F = annot[pd_fSamples].shape[0]
numPD_F


# In[45]:


PD_F =  annot[pd_fSamples]
PDF = PD_F.index
PF_F_LIST = PDF.tolist()
PD_FSET = F_HFrame[PF_F_LIST]


# In[46]:


#this is responsible for getting the number of samples of fibroblasts and healthy (control)
H_FSamples =  ((sample == "control")) & (tissue == "Fibroblasts")
numH_F = annot[H_FSamples].shape[0]
numH_F


# In[47]:


HD_F =  annot[H_FSamples]
HDF = HD_F.index
H_F_LIST = HDF.tolist()
FH_SET = F_HFrame[H_F_LIST]


# In[48]:


#Write a function to calculate log2 fold change using the two conditions provided to you. Use healthy controls as the denominator.
#You can use numpy.mean function to calculate mean expression and math.log function to calculate log2 value.
#You must handle situations where the denominator is 0 to avoid divide by zero errors. You can make all 0 mean expression values to 0.1
#You must handle situations where the fold change is 0 because log2 (0) is undefined. You can again substitute it to 0.1
def getLog2FoldChange(x, HFsamp=H_F_LIST, PDFsamp=PF_F_LIST):
    
    # get the expression values for each condition
    conditionHFsamp = x[HFsamp].values
    conditionPDFsamp = x[PDFsamp].values
    
    # get the means for the expression of each gene for each condition 
    conditionHF = np.mean(conditionHFsamp)
    conditionPDF = np.mean(conditionPDFsamp)

    # if the conditions are 0, then change it to 0.1 to avoid errors
    if (conditionPDF == 0):
        conditionPDF = 0.1
    # get the fold change from the quotient from the means from eahc condition
    foldChange = conditionHF / conditionPDF
     # if the foldChange is 0, then change it to 0.1 to avoid zero errors
    if foldChange == 0:
        foldChange = 0.1
    # return the log 2 of the fold change using the pandas log2 functions
    return np.log2(foldChange)

# apply the log2foldchange function to every value from the subsetted dataframe
log2FoldValues = F_HFrame.apply(getLog2FoldChange, axis=1)


# In[49]:


#this is the list of conditions 1 represents with PD and -1 represents as a control
condition = [1, -1, 1, -1, 1, -1, -1, 1,-1, 1, -1, -1,1, 1, -1]


# In[50]:


log2FoldValues 


# In[51]:


#t test to get corr coeff , t stat and p value
def mytstatistic(x, y = condition, n = numPD_F+numH_F):
    corrcoeff = np.corrcoef(x.values, y)[0,1]
    tstat = ((np.sqrt(n-2))/(np.sqrt(1-np.square(corrcoeff))))*corrcoeff
    p = 2*(1-scipy.stats.t.cdf(np.abs(tstat), n-2))
    return (p)


# In[52]:


#these are the columns that help get the condition 
F_cells.columns


# In[53]:


p = mytstatistic(F_cells,condition,numPD_F+numH_F)


# In[54]:


#this gets the pvalues of the fibro blast cells
F_cells_Stats = F_cells.apply(mytstatistic, axis=1)
F_cells_Stats


# In[82]:


#Find differentially expressed genes uing the criteria: |log2FC|>=2 and p-value<0.05
#Make a subset of the expression dataframe for just differentially expressed genes.
def getDifferentiallyExpressedGenes(log2FoldVals, geneStats=F_cells_Stats):
    # create a dataframe from the log2 fold change calculated above and the p-value
    df = pd.concat([log2FoldVals, geneStats], axis=1)
    # add columns for the new dataframe
    df.columns = ['log', 'p']
    
    # subsets the dataframe using the criteria before turning it into a list and returning it 
    return df.loc[(np.abs((df['log'])) >= 1.5) & (df['p'] < 0.01)].index.tolist()
    #return df.loc[(df['log'] >= 2) & (df['p'] < 0.05)]
        

#a list of the genes that are differentially expressed
diffExpressedIndexes = getDifferentiallyExpressedGenes(log2FoldValues, F_cells_Stats)

# gets the transposed version of the subsetted dataframe from part a because we need to filter based on genes, which are currently 
# the rows not columns 
transLog = F_cells.transpose()
# filter the dataframe based on the genes that are differentially expressed 
diffExpressed = transLog[diffExpressedIndexes]  #[list(diffExpressedIndexes)] 
#transposes the diffExpressed dataframe again to return it to the original state with genes as the rows
diffExpressed = diffExpressed.transpose()
seaborn.clustermap(diffExpressed)


# In[56]:


# this gets the fibroblasts cells into a dictionary with a key and values 
p_dict_F_cells = F_cells_Stats.to_dict()
p_dict_F_cells


# In[57]:


deg_list_F_cells = []
for key in p_dict_F_cells.keys():
    if p_dict_F_cells[key] < 0.01:
        deg_list_F_cells.append(key)


# In[58]:


#this prints out all the significant DEG that are within the fibroblasts cells
deg_list_F_cells


# In[59]:


FCset = F_cells.filter(items = deg_list_F_cells, axis=0)


# In[ ]:





# In[60]:


seaborn.clustermap(FCset)


# In[61]:


MBset = mb.filter(items = deg_list_mb, axis=0)


# In[62]:


seaborn.clustermap(MBset)


# In[63]:


SCset = sc.filter(items = deg_list_sc, axis=0)


# In[64]:


seaborn.clustermap(SCset)


# In[65]:


#mid brain and stem cells
def intersection (deg_list_mb, deg_list_sc):
    list_sc_mb = [value for value in deg_list_mb if value in deg_list_sc]
    return list_sc_mb


# In[66]:


#midbrain and stem cells
print (intersection(deg_list_mb, deg_list_sc))


# In[93]:


#midbrain and fibro
def intersection (deg_list_F_cells, deg_list_mb):
    list_mbbb= [value for value in deg_list_F_cells if value in deg_list_mb]
    return list_mbbb


# In[95]:


#midbrain and fibro
print (len(intersection(deg_list_F_cells, deg_list_mb)))


# In[78]:


def intersection (deg_list_F_cells, deg_list_sc):
    list_f_sc = [value for value in deg_list_F_cells if value in deg_list_sc]
    return list_f_sc


# In[79]:


#fibro and stem cells 
print(intersection(deg_list_F_cells, deg_list_sc))


# In[88]:


def intersection (deg_list_F_cells, deg_list_sc, deg_list_mb):
    list_f_sc_mb = [value for value in deg_list_F_cells if value in deg_list_sc if value in deg_list_mb]
    return list_f_sc_mb


# In[89]:


print(intersection(deg_list_F_cells, deg_list_sc, deg_list_mb))


# In[ ]:




