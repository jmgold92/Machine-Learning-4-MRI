
# coding: utf-8

# In[1]:

import numpy as np
import pandas as pd
from pylab import *
get_ipython().magic('matplotlib inline')
from mylocal_functions import *    


# In[2]:

# Make list of all DCE MRI files
dce_list = get_ipython().getoutput('ls ../Study_03_CBA/*T2*.txt')
df=Txt2DataFrame(dce_list[0])
dce_list[0].split("_")[-1].split(".")[0]


# In[ ]:


    for i in range(num_rows):
        l=all_lines[i].split()[1::3]
        np_array[i,:]=array(l)
    # Define Column Names
    column_names=list()
    for i in range(num_cols):
        column_names.append( "ROI_"+str(i+1))
    df=pd.DataFrame(np_array,columns=column_names)
    #df.plot()
    return df


# In[ ]:



