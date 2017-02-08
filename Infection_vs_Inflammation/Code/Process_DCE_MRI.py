# import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
#%matplotlib inline
# Define Function(s)
def Txt2DataFrame(file_name):
    f1 = open(file_name,'r')

    # Read all lines
    all_lines = f1.readlines()
    f1.close()

    # Select from Fourht line to the end
    first_row_with_data = 3
    all_lines=all_lines[first_row_with_data:]
    # Calculate number of rows
    num_rows=len(all_lines)
    f1.close()
    num_cols = len( all_lines[0].split()[1::3] )
    np_array=np.zeros( (num_rows,num_cols) )

    # allocate line by line
    for i in range(num_rows):
        l=all_lines[i].split()[1::3]
        np_array[i,:]=np.array(l)

    # Define Column Names
    column_names=list()
    for i in range(num_cols):
        column_names.append( "ROI_"+str(i+1))

    # Normalize
    for i in range(num_cols):
        np_array[:,i]= np_array[:,i]/np_array[0,i]

    df=pd.DataFrame(np_array,columns=column_names)
    
    Mouse=file_name.split("_")[3][1]; Mouse=int(Mouse)
    Slice=file_name.split("_")[4][1]; Slice=int(Slice)
    Year=file_name.split("_")[5][1];  Year=int(Year)
    Month=file_name.split("_")[6];    Month=int(Month)
    Day=file_name.split("_")[7];      Day=int(Day)
    
    df["Mouse"]=Mouse; df["Slice"]=Slice;
    df["Year"]=Year;   df["Day"]=Day
    return df

# Make list of all DCE MRI files
# This is equivalent to dce_list=!ls ../Study_03_CBA/*T2*.txt
# in the Jupyter Notebook

file_list = get_ipython().getoutput('ls ../Study_03_CBA/*T2.txt')
s0=Txt2DataFrame(file_list[0])
s1=Txt2DataFrame(file_list[1])
#plt.plot(s0,'--');plt.plot(s1,'o');


s0["Mouse"]=int(file_list[0].split("_")[3][1])





