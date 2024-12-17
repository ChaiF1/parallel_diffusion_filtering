# We will load in data from output.out using the pandas dataframe and make some plots
# This output.out has 5 columns "it, error, time, np, n" and a header
# The first column is the iteration number, the second column is the error, the third column is the time taken, the fourth column is the number of processors and the fifth column is the number of grid points


import pandas as pd

# Load in data
data = pd.read_csv("output.out", sep=",", header=0)

# Plot the time taken versus the number of pixels, and add for the amount of processors
import matplotlib.pyplot as plt

serial_data_np = [512,1024,2048,5096]
serial_data_time = [0.43, 2, 8.7, 39]

plt.figure()
for i in data["np"].unique():
	plt.plot(data[data["np"]==i]["n"], data[data["np"]==i]["time"], '.', label="np="+str(i))

# Plot the serial data with a red line
plt.plot(serial_data_np, serial_data_time, 'rx', label="Serial")
plt.xlabel("Number of pixels")
plt.ylabel("Time taken (s)")
plt.legend()
plt.savefig("time_vs_np.png")
