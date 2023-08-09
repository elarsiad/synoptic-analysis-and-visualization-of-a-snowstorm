import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
#import metview as mv
from metpy.units import units
import numpy as np
import xarray as xr
import os
os.environ['proj_lib'] = 'C:/Users/ejhas/anaconda3/pkgs/proj-9.0.0-h1cfcee9_1/Library/share/proj'

######################################################################
# Use Xarray to access GFS data from THREDDS resource and uses
# metpy accessor to parse file to make it easy to pull data using
# common coordinate names (e.g., vertical) and attach units.
era5_pres_fname = ('C:/Users/ejhas/OneDrive/Documents/BIMTEK/BMKG_OFS/era5.data_on_pressure_levels.2022022400_2022022521.nc')
ds = xr.open_dataset(era5_pres_fname)
ds_subset = ds.sel(time=np.datetime64('2022-02-25T12:00:00'), level=850).metpy.parse_cf()

# Subset data based on latitude and longitude values, calculate potential
# temperature for frontogenesis calculation.
#

# Set subset slice for the geographic extent of data to limit download
lon_slice = slice(200, 350)
lat_slice = slice(85, 10)

# Grab lat/lon values (GFS will be 1D)
lats = ds.latitude.sel(latitude=lat_slice).values
lons = ds.longitude.sel(longitude=lon_slice).values

#lats = ds.latitude.data
#lons = ds.longitude.data


level = 850 * units.hPa
#z850 = mpcalc.geopotential_to_height(ds_subset['z'])
#u850 = ds_subset['u']
#v850 = ds_subset['v']
#t850 = ds_subset['t']

zz850 = ds.Geopotential_height_isobaric.metpy.sel(
    vertical=level, latitude=lat_slice, longitude=lon_slice).metpy.unit_array.squeeze()
u850 = ds.Temperature_isobaric.metpy.sel(
    vertical=level, latitude=lat_slice, longitude=lon_slice).metpy.unit_array.squeeze()
v850 = ds['u-component_of_wind_isobaric'].metpy.sel(
    vertical=level, latitude=lat_slice, longitude=lon_slice).metpy.unit_array.squeeze()
t850 = ds['v-component_of_wind_isobaric'].metpy.sel(
    vertical=level, latitude=lat_slice, longitude=lon_slice).metpy.unit_array.squeeze()

# Convert temperatures to degree Celsius for plotting purposes
#tmpc_850 = t850.to('degC')

# Calculate potential temperature for frontogenesis calculation
#thta_850 = mpcalc.potential_temperature(level, t850)
thta_850 =mpcalc.potential_temperature(level, t850)

# Get a sensible datetime format
#vtime = ds.time.data[0].astype('datetime64[ms]').astype('O')
# Calculate frontogenesis
# -----------------------
#
# Frontogenesis calculation in MetPy requires temperature, wind
# components, and grid spacings. First compute the grid deltas using MetPy
# functionality, then put it all together in the frontogenesis function.
#
# Note: MetPy will give the output with SI units, but typically
# frontogenesis (read: GEMPAK) output this variable with units of K per
# 100 km per 3 h; a conversion factor is included here to use at plot time
# to reflect those units.
#
dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

fronto_850 = mpcalc.frontogenesis(thta_850, u850, v850, dx, dy, dim_order='yx')

# A conversion factor to get frontogensis units of K per 100 km per 3 h
convert_to_per_100km_3h = 1000*100*3600*3

######################################################################
# Plotting Frontogenesis
# ----------------------
#
# Using a Lambert Conformal projection from Cartopy to plot 850-hPa
# variables including frontogenesis.
#
mapcrs = ccrs.LambertConformal(central_longitude=-110.0, central_latitude=35.0)
# Set projection of the data (GFS is lat/lon)
datacrs = ccrs.PlateCarree()


# Start figure and limit the graphical area extent
fig = plt.figure(1, figsize=(14, 12))
ax = plt.subplot(111, projection=mapcrs)
ax.set_extent([-180, -30, 10, 70], ccrs.PlateCarree())

# Add map features of Coastlines and States
ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))

# Plot 850-hPa Frontogenesis
clevs_tmpc = np.arange(-40, 41, 2)
cf = ax.contourf('longitude', 'latitude', fronto_850*convert_to_per_100km_3h, np.arange(-8, 8.5, 0.5),
                 cmap=plt.cm.bwr, extend='both', transform=datacrs)
cb = plt.colorbar(cf, orientation='horizontal', pad=0, aspect=50, extendrect=True)
cb.set_label('Frontogenesis K / 100 km / 3 h')

# Plot 850-hPa Geopotential Heights
clevs_850_hght = np.arange(0, 8000, 30)
cs = ax.contour('longitude', 'latitude', z850, clevs_850_hght, colors='black', transform=datacrs)
plt.clabel(cs, fmt='%d')

# Plot 850-hPa Wind Barbs only plotting every fifth barb
wind_slice = (slice(None, None, 5), slice(None, None, 5))
ax.barbs('longitude'[wind_slice[0]], 'latitude'[wind_slice[1]],
         u850[wind_slice].to('kt').m, v850[wind_slice].to('kt').m,
         color='black', transform=datacrs)

# Plot some titles
plt.title('ERA5 850-hPa Geopotential Heights (m), Temp (C), and Winds', loc='left')
plt.title('Valid Time: {}'.format(vtime), loc='right')


plt.show()