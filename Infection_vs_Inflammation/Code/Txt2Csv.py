def Txt2Csv(main_dir,file_name):
    f1 = open(main_dir + file_name,'r')
    # Read all lines
    all_lines = f1.readlines()
    # Select from Fourht line to the end
    first_row_with_data = 3
    all_lines=all_lines[first_row_with_data:]
    # Calculate number of rows
    num_rows=len(all_lines)
    # Use first line to estimate number of columns
    num_columns = len(all_lines[0].split()[1::3])
    # Pre allocate Array
    Data_Array = np.zeros( (num_rows,num_columns) )
    Data_Array.shape
    # Allocate Array
    for i in range(num_rows):
        Data_Array[i-first_row_with_data]=np.array(all_lines[i].split()[1::3])

    # Create a DataFrame and save it as CSV
    csv_name=file_name[0:-4] + ".csv"
    # Create Heading
    df=pd.DataFrame()
    for i in range(Data_Array.shape[1]):
        df["ROI_" + str(i+1)]=Data_Array[:,i]
        # Write File
        df.to_csv(csv_name)
     return df
