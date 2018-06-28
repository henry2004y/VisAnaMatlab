# Code for 'Script'
# We will use NumPy to read the csv file.
# Refer to NumPy documentation for genfromtxt() for details on
# customizing the CSV file parsing.
import numpy as np
# assuming data.csv is a CSV file with the 1st row being the names names for
# the columns
data = np.genfromtxt("test.csv", dtype=None, names=True, delimiter=',', autostrip=True)

for name in data.dtype.names:
    array = data[name]
    
    # You can directly pass a NumPy array to the pipeline.
    # Since ParaView expects all arrays to be named, you
    # need to assign it a name in the 'append' call.
    output.RowData.append(array , name)

