# Define Function(s)
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def txt_2_array(file_name):
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
    data_array=np.zeros( (num_rows,num_cols) )

    # allocate line by line
    for i in range(num_rows):
        l=all_lines[i].split()[1::3]
        data_array[i,:]=np.array(l)

    return data_array

def name_2_dict(file_name):
    # create dictionary with information
    file_info = {'Mouse': file_name.split("_")[3], \
                'Slice':  int(file_name.split("_")[4][1]),\
                'Year':   int(file_name.split("_")[5]),\
                'Month':    int(file_name.split("_")[6]),\
                'Day':    int(file_name.split("_")[7])}
    return file_info

def T2decay(xdata, T1,Mz):
    R1=1.0/T1
    return Mz * np.exp(-xdata*R1)

def fitT2(xdata,Ydata):
    num_rois=Ydata.shape[1]
    initial_guess=[1,.1]
    T2_estimates=np.zeros((num_rois,1))

    for i in range(num_rois):
        ydata=Ydata[:,i] / Ydata[0,i]
        pars_hat, cov = curve_fit(T2decay, xdata, ydata,p0=initial_guess,bounds=(0, 1.0))
        T2_estimates[i,0]=pars_hat[0]
    return T2_estimates

def name_2_df(file_name):
    # create dataframe with information from file name
    my_dict=name_2_dict(file_name)
    my_df=pd.DataFrame()
    my_df=my_df.from_dict(my_dict,orient='index').transpose()
    return my_df


def Lorentzian(sat_offset,Amp,Width,Center):
    Width = Width**2; Width=Width/4
    xdata = (sat_offset-Center)**2
    return (Amp*Width) / (Width +xdata )

def Lorentzian2(sat_offset,a1,w1,c1,a2,w2,c2):
    return Lorentzian(sat_offset,a1,w1,c1) + Lorentzian(sat_offset,a2,w2,c2)


def Lscale(sat_offset,a1,w1,c1,a2,w2,c2,scale):
    return Lorentzian2(sat_offset,a1,w1,c1,a2,w2,c2) + scale

def fit_L2_scale(offsets,ydata):
    Signal=1-ydata
    # fix xdata
    xdata=offsets-offsets[Signal.argmax()]
    # allocate fitting based on this
    A10, W10, C10 =  0.90, 1, 0
    A20, W20, C20 =  .1, 1, -4

    A1L, W1L, C1L =  0.5, .1, -.1
    A2L, W2L, C2L =  0, .1,   -6

    A1U, W1U, C1U =  1.0, 5, +.1
    A2U, W2U, C2U =  1.0, 5, -1.0

    scale0, scaleL, scaleU = 0, -1, +1

    initial_guess = [A10, W10, C10,      A20, W20, C20,    scale0]
    lb            = [A1L, W1L, C1L,      A2L, W2L, C2L,    scaleL]
    ub            = [A1U, W1U, C1U,      A2U, W2U, C2U,    scaleU]

    p, cov = curve_fit(Lscale, xdata, Signal,p0=initial_guess,bounds=(lb, ub))

    return p;
