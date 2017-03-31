# Import Modules as needed
import numpy as np
from mylocal_functions import *
import matplotlib.pyplot as plt

# ======== CEST============= #
CEST_list=get_ipython().getoutput('ls ../Study_03_CBA/*CEST.txt')
CEST_Int_matrix=np.zeros((len(CEST_list),4))
NOE=np.zeros((1,4)); NOEmatrix=np.zeros_like(CEST_Int_matrix)
NOE_w=np.zeros_like(NOE);  NOEmatrix_w=np.zeros_like(NOEmatrix)
ppm=np.linspace(-8,8,101);  


for i in range( len(CEST_list) ):
    D=txt_2_array(CEST_list[i]);       #Convert txt file to array
    Zn=normalize_data(D.T,8);          Zn=Zn[:,9::]
    M=np.zeros([1,4])
    
    for j in range(4):
        p=fit_L3_scale(ppm,Zn[j,:])
        L=Lscale3(ppm,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9]);
        #CEST_centered[i,:]=L
        #CEST_integral[i,0]=np.sum(L)
        M[0,j]=np.sum(L)
        NOE[0,j]=p[3]
        NOE_w[0,j]=p[4]
        
    CEST_Int_matrix[i,:]=M
    NOEmatrix[i,:]=NOE
    NOEmatrix_w[i,:]=NOE_w

#======== create violing plots ============= #
Tissues=["Infected","Healthy R","Sterile Infl.","Healthy L"]

# Set dimensions of plot
fig = plt.figure(1,figsize=(10,10));
# CEST
ax = fig.add_subplot(3,1,1);  ax.set_xticks([1, 2, 3, 4]);  ax.set_xticklabels(Tissues)
plt.violinplot(CEST_Int_matrix, showextrema=True,showmedians=True);
plt.ylabel("CEST Integral")
#T2
ax = fig.add_subplot(3,1,2);  ax.set_xticks([1, 2, 3, 4]);  ax.set_xticklabels(Tissues)
plt.violinplot(NOEmatrix_w,showextrema=True,showmedians=True);
plt.violinplot(NOEmatrix_w,showextrema=True,showmedians=True);
plt.violinplot(NOEmatrix_w,showextrema=True,showmedians=True);
plt.ylabel("NOEmatrix_w")
#T2ex
ax = fig.add_subplot(3,1,3);  ax.set_xticks([1, 2, 3, 4]);  ax.set_xticklabels(Tissues)
plt.violinplot(NOEmatrix,showextrema=True,showmedians=True);
plt.violinplot(NOEmatrix,showextrema=True,showmedians=True);
plt.ylabel("NOE Amp")




