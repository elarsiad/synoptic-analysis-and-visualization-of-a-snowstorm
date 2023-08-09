#!/usr/bin/env python
# coding: utf-8

# In[84]:


import xarray as xr
import numpy as np
import metpy
import metpy.calc as mpcalc
import metpy.units as units
from datetime import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches

from metpy.units import units
#from siphon.catalog import TDSCatalog
#from siphon.ncss import NCSS

# In[2]:


# Case Study Date
year = 2022
month = 2
day = 24
hour = 0

dt = datetime(year, month, day, hour)

# In[26]:


# Open data with xarray, and parse it with MetPy
# ds = xr.open_dataset(xr.backends.NetCDF4DataStore(data)).metpy.parse_cf()

era5_pres_fname = ('C:/Users/ejhas/OneDrive/Documents/BIMTEK/BMKG_OFS/era5.data_on_pressure_levels.2022022400_2022022521.nc')
ds = xr.open_dataset(era5_pres_fname).metpy.parse_cf()

ds

# In[27]:


vtime = ds['time'].metpy.time[0]
vtime

# In[31]:


# This is the time we're using
vtime = ds['time'].metpy.time[0]

# Grab lat/lon values from file as unit arrays
lats = ds['latitude'].metpy.unit_array
lons = ds['longitude'].metpy.unit_array

# Calculate distance between grid points
# will need for computations later
dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

# Grabbing data for specific variable contained in file (as a unit array)
# 700 hPa Geopotential Heights
# hght_400 = ds['z'].metpy.sel(level=300 * units.hPa, time=vtime)


# 400 hPa Temperature
tmpk_400 = ds['t'].metpy.sel(level=400 * units.hPa, time=vtime)

# 400 hPa u-component_of_wind
uwnd_400 = ds['u'].metpy.sel(level=400 * units.hPa, time=vtime)

# 400 hPa v-component_of_wind
vwnd_400 = ds['v'].metpy.sel(level=400 * units.hPa, time=vtime)

# 400 hPa Geopotential Heights
hght_500 = ds['z'].metpy.sel(level=500 * units.hPa, time=vtime)

# 500 hPa u-component_of_wind
uwnd_500 = ds['u'].metpy.sel(level=500 * units.hPa, time=vtime)

# 500 hPa v-component_of_wind
vwnd_500 = ds['v'].metpy.sel(level=500 * units.hPa, time=vtime)

# 600 hPa u-component_of_wind
uwnd_600 = ds['u'].metpy.sel(level=600 * units.hPa, time=vtime)

# 600 hPa v-component_of_wind
vwnd_600 = ds['v'].metpy.sel(level=600 * units.hPa, time=vtime)

# In[6]:


# Set constant values that will be needed in computations

# Set default static stability value
sigma = 2.0e-6 * units('m^2 Pa^-2 s^-2')

# Set f-plane at typical synoptic f0 value
f0 = 1e-4 * units('s^-1')

# Use dry gas constant from MetPy constants
Rd = mpconstants.Rd

# In[7]:


# Smooth Heights
# For calculation purposes we want to smooth our variables
# a little to get to the "synoptic values" from higher
# resolution datasets

# Number of repetitions of smoothing function
n_reps = 50

# Apply the 9-point smoother
hght_700s = mpcalc.smooth_n_point(hght_700, 9, n_reps)  # .metpy.unit_array
hght_500s = mpcalc.smooth_n_point(hght_500, 9, n_reps)  # .metpy.unit_array

tmpk_700s = mpcalc.smooth_n_point(tmpk_700, 9, n_reps)  # .metpy.unit_array
tmpc_700s = tmpk_700s.metpy.convert_units('degC')

uwnd_700s = mpcalc.smooth_n_point(uwnd_700, 9, n_reps)  # .metpy.unit_array
vwnd_700s = mpcalc.smooth_n_point(vwnd_700, 9, n_reps)  # .metpy.unit_array

uwnd_500s = mpcalc.smooth_n_point(uwnd_500, 9, n_reps)  # .metpy.unit_array
vwnd_500s = mpcalc.smooth_n_point(vwnd_500, 9, n_reps)  # .metpy.unit_array

uwnd_900s = mpcalc.smooth_n_point(uwnd_900, 9, n_reps)  # .metpy.unit_array
vwnd_900s = mpcalc.smooth_n_point(vwnd_900, 9, n_reps)  # .metpy.unit_array

# In[11]:


# Absolute Vorticity Calculation
avor_900 = mpcalc.absolute_vorticity(uwnd_900s, vwnd_900s, dx, dy, lats)
avor_500 = mpcalc.absolute_vorticity(uwnd_500s, vwnd_500s, dx, dy, lats)

# Advection of Absolute Vorticity
vortadv_900 = mpcalc.advection(avor_900, (uwnd_900s, vwnd_900s), (dx, dy))
vortadv_500 = mpcalc.advection(avor_500, (uwnd_500s, vwnd_500s), (dx, dy))

# Differential Vorticity Advection between two levels
diff_avor = ((vortadv_900 - vortadv_500) / (400 * units.hPa))

# Calculation of final differential vorticity advection term
term_A = (-f0 / sigma * diff_avor)
print(term_A.units)

# In[ ]:


# In[ ]:


# In[ ]:




