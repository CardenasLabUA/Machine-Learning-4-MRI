# Import Modules as needed
import numpy as np
from mylocal_functions import *
import matplotlib.pyplot as plt

# ======== CEST============= #
CEST_list=get_ipython().getoutput('ls ../Study_03_CBA/*CEST.txt')
CEST_Int_matrix=np.zeros((len(CEST_list),4))

ppm=np.linspace(-8,8,101);  


for i in range( len(CEST_list) ):
    D=txt_2_array(CEST_list[i]);       #Convert txt file to array
    Zn=normalize_data(D.T,8);          Zn=Zn[:,9::]
    M=np.zeros([1,4])
    
    for j in range(4):
        p=fit_L2_scale(ppm,Zn[j,:])
        L=Lscale(ppm,p[0],p[1],p[2],p[3],p[4],p[5],p[6]);
        #CEST_centered[i,:]=L
        #CEST_integral[i,0]=np.sum(L)
        M[0,j]=np.sum(L)
        
    CEST_Int_matrix[i,:]=M

# ======== T2 MSME============= #a
# Make list of all T2.txt files
T2_list = get_ipython().getoutput('ls ../Study_03_CBA/*T2.txt')
T2_matrix=np.zeros( (len(T2_list),4) )

TR=np.linspace(.012,.012*12,12)

# Fit T2 
for i in range(len(T2_list)):
    YDataMatrix=txt_2_array(T2_list[i])
    #Estimate T2
    T2time=fitT2(TR,YDataMatrix)
    T2_matrix[i,:]=T2time.T

# ======== T2ex DCE============= #
# Make list of all T2.txt files
T2ex_list = get_ipython().getoutput('ls ../Study_03_CBA/*T2exDCE.txt')
T2ex_Int_matrix=np.zeros( (len(T2ex_list),4) )

# T2ex integral
for i in range( len(T2ex_list) ):
    D=txt_2_array(T2ex_list[i]);       #Convert txt file to array
    Zn=normalize_data(D.T,0);          Zn=Zn[:,9::]
    T2ex_Int_matrix[i,:]=np.sum(Zn-1,axis=1)


#======== create violing plots ============= #
Tissues=["Infected","Healthy R","Sterile Infl.","Healthy K"]

# Set dimensions of plot
fig = plt.figure(1,figsize=(10,10));
# CEST
ax = fig.add_subplot(3,1,1);  ax.set_xticks([1, 2, 3, 4]);  ax.set_xticklabels(Tissues)
plt.violinplot(CEST_Int_matrix, showextrema=True,showmedians=True);
plt.ylabel("CEST Integral")
#T2
ax = fig.add_subplot(3,1,2);  ax.set_xticks([1, 2, 3, 4]);  ax.set_xticklabels(Tissues)
plt.violinplot(T2_matrix,showextrema=True,showmedians=True);
plt.violinplot(T2_matrix,showextrema=True,showmedians=True);
plt.violinplot(T2_matrix,showextrema=True,showmedians=True);
plt.ylabel("T2 time")
#T2ex
ax = fig.add_subplot(3,1,3);  ax.set_xticks([1, 2, 3, 4]);  ax.set_xticklabels(Tissues)
plt.violinplot(T2ex_Int_matrix,showextrema=True,showmedians=True);
plt.violinplot(T2ex_Int_matrix,showextrema=True,showmedians=True);
plt.ylabel("T2ex Integral")

# plot non_neg only

# ======== T2ex DCE ALL============= #

# Make list of all T2.txt files

def plotvio(slice_num):
    p1='ls ../Study_03_CBA/*S'
    p2=str(slice_num)
    p3='*T2exDCE.txt'
    file_names=p1+p2+p3
    T2ex_list = get_ipython().getoutput(file_names)
    T2ex_Int_matrix=np.zeros( (len(T2ex_list),4) )
# T2ex integral
    for i in range( len(T2ex_list) ):
        D=txt_2_array(T2ex_list[i]);       #Convert txt file to array
        Zn=normalize_data(D.T,0);          #Zn=Zn[:,9::]
        T2ex_Int_matrix[i,:]=np.sum(Zn-1,axis=1)
    
    plt.violinplot(T2ex_Int_matrix,showextrema=True,showmedians=True);
        
for i in range(5):
    n=i+1
    plt.figure(99,figsize=(15,15));
    plt.subplot(5,1,n); plt.title("Slice_0"+str(n))
    plotvio(n)


