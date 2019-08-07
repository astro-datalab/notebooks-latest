

# Load modules
import numpy as np
import astropy.units as u


# Define a function that reads in the data file
# USAGE :   time, total, data = Read("filename")
def Read(filename):

	# open the file 
	file = open(filename,'r')
    
	#read header info line by line (line will be a string)
	# read first two lines FIRST and store as variable
    
	# read in first line and store time
	line1 = file.readline()
	label, value = line1.split()
	time = float(value)*u.Myr

	# read in 2nd line and store total number of particles
	line2 = file.readline()
	label, value = line2.split()
	total = float(value)
    
	# close file
	file.close()
	
	# read the remainder of the file, 
	# "dtype=None" means line is split using white spaces
	# "skip_header=3"  skipping the first 3 lines 
	# the flag "names=True" creates arrays to store the date
	# with the column headers given in line 4 like "m", "x"
	data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
	# return the time of the snapshot, 
	# total number of particles 
	#and an array that stores the remainder of the data. 
	return time, total, data






