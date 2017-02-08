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
    Year=file_name.split("_")[5]  
    Month=file_name.split("_")[6];    Month=int(Month)
    Day=file_name.split("_")[7];      Day=int(Day)
    
    df["Mouse"]=Mouse; df["Slice"]=Slice;
    df["Day"]=Day; df["Month"]=Month; df["Year"]=Year;   
    
    return df
#  ========= DCE MRI ================
# List of all files
dce_file_list = get_ipython().getoutput('ls ../Study_03_CBA/*T2.txt')

# Create Massive DataFrame
DF=Txt2DataFrame(dce_file_list[0])
TE=np.linspace(.12,.144,12)
DF["TE"]=TE

for i in range(1,len(dce_file_list)):
    print(dce_file_list[i])
    df=Txt2DataFrame(dce_file_list[i]); df["TE"]=TE
    DF=DF.append(df)
    
DF=DF.set_index("TE")


#  ========= CEST MRI ================
# List of all files
cest_file_list = get_ipython().getoutput('ls ../Study_03_CBA/*CEST.txt')

# Create Massive DataFrame
DF_CEST=Txt2DataFrame(cest_file_list[0])
for i in range(1,len(cest_file_list)):
    df=Txt2DataFrame(cest_file_list[i]); 
    DF_CEST=DF_CEST.append(df)



plt.plot(DF_CEST.iloc[:,[0,1,2,3]].as_matrix())


