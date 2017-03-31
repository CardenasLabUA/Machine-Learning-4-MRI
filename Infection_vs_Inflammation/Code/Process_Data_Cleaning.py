# Import Modules as needed
import numpy as np
#import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt
from mylocal_functions import *

# ======== T2 MSME============= #
# Make list of all T2.txt files
T2_list = get_ipython().getoutput('ls ../Study_03_CBA/*T2.txt')

# Allocate variables needed for analysis
T2DF = pd.DataFrame()
TR = np.linspace(.012,.012*12,12)
# Fit T2 and construct dataframe
for names in T2_list:
    #Convert txt file to array
    YDataMatrix = txt_2_array(names)
    #Estimate T2
    T2time = fitT2(TR,YDataMatrix)
    #convert to data frame
    df_T2 = pd.DataFrame(T2time.T,columns = ["Infected","Healthy_R","St_Inf","Healthy_L"])
    #df_T2=pd.DataFrame(T2time.T,columns=["ROI-1","ROI-2","ROI-3","ROI-4"])
    df_info = name_2_df(names)
    df_final = pd.concat([df_T2,df_info], axis = 1)
    T2DF = T2DF.append(df_final,ignore_index=True)

# Plot T2 Density ROIs 1 and 2
#T2DF[T2DF.Slice==1].iloc[:,:4].plot.density(); title("Slice 01"); xlim((0.025,.15))
#T2DF[T2DF.Slice==2].iloc[:,:4].plot.density(); title("Slice 02"); xlim((0.025,.15))
#T2DF[T2DF.Slice==3].iloc[:,:4].plot.density(); title("Slice 03"); xlim((0.025,.15))
#T2DF[T2DF.Slice==4].iloc[:,:4].plot.density(); title("Slice 04"); xlim((0.025,.15))
#T2DF[T2DF.Slice==5].iloc[:,:4].plot.density(); title("Slice 05"); xlim((0.025,.15))


# ======== CEST============= #
# Make list of all T2.txt files
CEST_list = get_ipython().getoutput('ls ../Study_03_CBA/*CEST.txt')
CEST_DF = pd.DataFrame()
Z = np.zeros((4,110))

def normalize_data(DataMatrix):
    rows,cols = DataMatrix.shape
    newData = np.zeros_like(DataMatrix)
    for row in range(rows):
        newData[row,:] = DataMatrix[row,:]/DataMatrix[row,8]
    return newData

for names in CEST_list:
    #Convert txt file to array
    D = txt_2_array(names);
    Zn = normalize_data(D.T)
    Z = np.concatenate((Z,Zn))

Z = Z[4::,9::]

# define offsets in ppm
a1 = np.linspace(-55,-50,9)
ppm = np.linspace(-8,8,101)
full_ppm = np.concatenate((a1, ppm))


# fit CEST data.
y = Z[12,:]
p = fit_L2_scale(ppm,y)
Yhat = Lscale(ppm,p[0],p[1],p[2],p[3],p[4],p[5],p[6]);

plt.figure(figsize = (10,6))
plt.plot(ppm,y,'o',label='Signal');
plt.plot(ppm,1-Yhat,'-',label='Fit');
plt.legend()

## ======   BUILD CEST Predictors ======== #####
CEST_predictors = np.zeros_like(Z)
rows,cols = CEST_predictors.shape
Tissue_Class = np.zeros((4,rows))

for i in range(rows):
    p = fit_L2_scale(ppm,Z[i,:])
    CEST_predictors[i,:] = Lscale(ppm,p[0],p[1],p[2],p[3],p[4],p[5],p[6]);

Tissue_Class = np.zeros((64,1))

for i in range(4):
    Tissue_Class[i::4]=i

CEST_Dataframe=pd.DataFrame(CEST_predictors)
CEST_Dataframe["Tissue_Class"]=Tissue_Class

pd.DataFrame.to_csv(CEST_Dataframe,"CEST_infections.csv",header=True,index=False)
