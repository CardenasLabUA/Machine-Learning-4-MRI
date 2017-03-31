# Import Modules as needed
import numpy as np
#import seaborn as sn
import pandas as pd
from mylocal_functions import *
import matplotlib.pyplot as plt

# ======== CEST============= #
# Make list of all CEST.txt files
CEST_list=get_ipython().getoutput('ls ../Study_03_CBA/*CEST.txt')

# create function to normalize data
def normalize_data(DataMatrix):
    rows,cols = DataMatrix.shape
    newData = np.zeros_like(DataMatrix)
    for row in range(rows):
        newData[row,:]=DataMatrix[row,:]/DataMatrix[row,8]
    return newData

# extract CEST data as a 4 X 110 matrix
Z=np.zeros((4,110))
for names in CEST_list:
    D=txt_2_array(names);       #Convert txt file to array
    Zn=normalize_data(D.T)
    Z=np.concatenate((Z,Zn))
Z=Z[4::,9::]

# define offsets in ppm
a1=np.linspace(-55,-50,9);  ppm=np.linspace(-8,8,101);  full_ppm = np.concatenate((a1, ppm))


# Fit data to center it
CEST_centered=np.zeros_like(Z); rows,cols = CEST_centered.shape;
CEST_integral=np.zeros((rows,1))

for i in range(rows):
    p=fit_L2_scale(ppm,Z[i,:])
    L=Lscale(ppm,p[0],p[1],p[2],p[3],p[4],p[5],p[6]);
    CEST_centered[i,:]=L
    CEST_integral[i,0]=np.sum(L)
# create tissue label    
TissueLabel=np.zeros((rows,1))
for i in range(1,4):
    TissueLabel[i::4]=i

CEST_integral_df= pd.DataFrame(data=CEST_integral,columns=["CEST_integral"]); 
CEST_integral_df["Tissue"]=TissueLabel

Y=np.zeros((16,4))
for i in range(4):
    df=CEST_integral_df[CEST_integral_df["Tissue"]==i]
    Y[:,i]=df.CEST_integral.values;
    



# ======== T2 MSME============= #
# Make list of all T2.txt files
T2_list = get_ipython().getoutput('ls ../Study_03_CBA/*S3*T2.txt')
T2_matrix=np.zeros( (len(T2_list),4) )

TR=np.linspace(.012,.012*12,12)

# Fit T2 
for i in range(len(T2_list)):
    YDataMatrix=txt_2_array(T2_list[i])
    #Estimate T2
    T2time=fitT2(TR,YDataMatrix)
    T2_matrix[i,:]=T2time.T
    
#======== create violing plots ============= #
Tissues=["Infected","Healthy R","Sterile Infl.","Healthy K"]
fig = plt.figure(); ax = fig.add_subplot(111)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(Tissues)
plt.violinplot(Y); plt.ylabel("CEST Integral")

# create violing plot
fig = plt.figure(); ax = fig.add_subplot(111)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(Tissues)
plt.violinplot(T2_matrix);  plt.ylabel("T2 time")

