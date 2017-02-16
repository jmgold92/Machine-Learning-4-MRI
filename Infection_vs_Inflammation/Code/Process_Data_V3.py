# Import Modules as needed
import numpy as np
from mylocal_functions import *
import matplotlib.pyplot as plt

# ======== CEST============= #
CEST_list=get_ipython().getoutput('ls ../Study_03_CBA/*CEST.txt')
CEST_Int_matrix=np.zeros((len(CEST_list),4))

for i in range( len(CEST_list) ):
    D=txt_2_array(CEST_list[i]);       #Convert txt file to array
    Zn=normalize_data(D.T,8);          Zn=Zn[:,9::]
    CEST_Int_matrix[i,:]=np.sum(Zn,axis=1) 

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

# ======== T2ex DCE============= #
# Make list of all T2.txt files
T2ex_list = get_ipython().getoutput('ls ../Study_03_CBA/*S3*T2exDCE.txt')
T2ex_Int_matrix=np.zeros( (len(T2ex_list),4) )

# T2ex integral
for i in range( len(T2ex_list) ):
    D=txt_2_array(T2ex_list[i]);       #Convert txt file to array
    Zn=normalize_data(D.T,0);          #Zn=Zn[:,9::]
    T2ex_Int_matrix[i,:]=np.sum(Zn-1,axis=1)


#======== create violing plots ============= #
Tissues=["Infected","Healthy R","Sterile Infl.","Healthy K"]

# Set dimensions of plot
fig = plt.figure(10,figsize=(10,10));
# CEST
ax = fig.add_subplot(3,1,1);  ax.set_xticks([1, 2, 3, 4]);  ax.set_xticklabels(Tissues)
plt.violinplot(CEST_Int_matrix, showextrema=True,showmedians=True);
plt.ylabel("CEST Integral")
#T2
ax = fig.add_subplot(3,1,2);  ax.set_xticks([1, 2, 3, 4]);  ax.set_xticklabels(Tissues)
plt.violinplot(T2_matrix,showextrema=True,showmedians=True);
plt.ylabel("T2 time")
#T2ex
ax = fig.add_subplot(3,1,3);  ax.set_xticks([1, 2, 3, 4]);  ax.set_xticklabels(Tissues)
plt.violinplot(T2ex_Int_matrix,showextrema=True,showmedians=True);
plt.ylabel("T2ex Integral")

p=plt.violinplot(T2ex_Int_matrix,showextrema=True,showmedians=True);
plt.ylabel("T2ex Integral")
