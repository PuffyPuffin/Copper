# Import packages
import numpy as np
from netCDF4 import Dataset as nc
import math
import statistics as stat
import pandas as pd
import xarray as xr
import csv

#Import files from computer
## file 1 = no sinking (ns)
## file 2 = sinking (s)
file1 = 'data/release_2018.05.01.nc'
file2 = 'data/release_2018.05.01.nc'
# Read netCDF4 files for sticking
ns = nc(file1)
s = nc(file2)

# Read netCDF4 files for pandas
pd_ns = xr.open_dataset(file1, decode_times=False)
pd_s = xr.open_dataset(file2, decode_times=False)

# Convert relative depth float64 to data files 
ns_reldepth = ns["cs"][:].data
s_reldepth = s["cs"][:].data

# Filter data by relative depth <= -.0.99
ns_depthB = ns_reldepth <= -0.99
# Find indices relating to first instance of reaching <= -0.99
ns_ind = np.argmax(ns_depthB, axis = 0)
# For particles that never reach <= -0.99, set index to final value
ns_ind[ns_ind==0] = 720
# Convert to tuple for faster iteration
ns_ind = tuple(ns_ind)

# Repeat for sinking
s_depthB = s_reldepth <= -0.99
s_ind = np.argmax(s_depthB, axis = 0)
s_ind[s_ind==0] = 720
s_ind = tuple(s_ind)

# Covert read files to pandas dataframe
pd_ns = pd_ns.to_dataframe()
pd_s = pd_s.to_dataframe()

#Select latitude and longitude from variables
pd_ns = pd_ns.take([0,1], axis=1)
pd_s = pd_s.take([0,1], axis=1)
# Create variables for # of coulmns and rows of interest
n, m = 1000, 720
# Create empty arrays to fill in for loops
## For distance between each time point
ns_dist, s_dist = np.zeros((n,m)), np.zeros((n,m))
## For total distance before sticking
ns_sum_dist, s_sum_dist = np.zeros((n,1)), np.zeros((n,1))
# Create nested for loop to calculate distance (km) traveled between each time point 
## no sink
for i in range(n): # columns
    for j in range(m): # rows
        start_latt, end_latt = pd_ns.loc[(i,j),"lat"], pd_ns.loc[(i,j+1),"lat"] # set the two latitudes
        start_long, end_long = pd_ns.loc[(i,j),"lon"], pd_ns.loc[(i,j+1),"lon"] # set two longitudes
        d_long, d_latt = (end_long - start_long)*math.pi/180, (end_latt - start_latt)*math.pi/180 # find difference and convert to radians
        a = math.sin(d_latt/2)**2 + math.cos(start_latt) * math.cos(end_latt) * math.sin(d_long/2)**2 # part 1 cosine-Haversine formula
        c = 2 * math.asin(math.sqrt(a)) # part 2 cosine-Haversine formula
        ns_dist[i,j] = 6371 * c # part 3 cosine-haversine + save calculated distance to array
# Create for loop that sums distances before first instance of relative depth <= -.99       
for i in range(n):
    ns_sum_dist[i] = ns_dist[i,0:ns_ind[i]].sum()
# Create nested for loop to calculate distance (km) traveled between each time point 
## sink
for i in range(n): # columns
    for j in range(m): # rows
        start_latt, end_latt = pd_s.loc[(i,j),"lat"], pd_s.loc[(i,j+1),"lat"] # set the two latitudes
        start_long, end_long = pd_s.loc[(i,j),"lon"], pd_s.loc[(i,j+1),"lon"] # set two longitudes
        d_long, d_latt = (end_long - start_long)*math.pi/180, (end_latt - start_latt)*math.pi/180 # find difference and convert to radians
        a = math.sin(d_latt/2)**2 + math.cos(start_latt) * math.cos(end_latt) * math.sin(d_long/2)**2 # part 1 cosine-Haversine formula
        c = 2 * math.asin(math.sqrt(a)) # part 2 cosine-Haversine formula
        s_dist[i,j] = 6371 * c # part 3 cosine-haversine + save calculated distance to array
# Create for loop that sums distances before first instance of relative depth <= -.99 
for i in range(n):
    s_sum_dist[i] = s_dist[i,0:s_ind[i]].sum()
#Create column names
cols = ["No Sinking", "Sinking"]
# Convert summmed distance to pandas dataframe and get rid of brackets
distance = pd.DataFrame(list(zip(ns_sum_dist.flatten(), s_sum_dist.flatten())))
# Save data to csv in Data/ folder with following naming scheme:
## MonthYearLocation_dist.csv
### Month = 3 letter acronymn, ex. January = Jan
### Year = last two digits unless pre-2000, ex. 2018 = 18
### Location = shortened name, ex. Puget Sound Naval Shipyard = Navy
distance.to_csv("data/May18Navy_dist.csv", header = cols)

