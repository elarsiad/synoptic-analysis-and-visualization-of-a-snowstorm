#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:

##Plotting 850 geopotential height, temperature (F) and wind (kt)
era5_pres_fname = ('C:/Users/ejhas/OneDrive/Documents/BIMTEK/BMKG_OFS/era5.data_on_pressure_levels.2022022400_2022022521.nc')
ds1 = xr.open_dataset(era5_pres_fname).metpy.parse_cf()


# In[3]:


ds1


# In[4]:


proj = ccrs.LambertConformal(central_longitude=-110.0, central_latitude=35.0)
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection=proj)
ax.stock_img()
ax.add_patch(
mpatches.Rectangle(xy=[-170, 20], width=120, height=45, facecolor='none',
edgecolor='red', transform=ccrs.PlateCarree()))
ax.set_extent([-180, -30, 10, 70], crs=ccrs.PlateCarree())
ax.gridlines()
ax.coastlines()


# In[5]:


ds1_subset = ds1.sel(time=np.datetime64('2022-02-25T15:00:00'), level=850).metpy.parse_cf()


# In[6]:


ds1_subset


# In[7]:


z850 = mpcalc.geopotential_to_height(ds1_subset['z'])
u850 = ds1_subset['u'].metpy.quantify()
v850 = ds1_subset['v'].metpy.quantify()
t850 = ds1_subset['t'].metpy.quantify()
a850 = mpcalc.absolute_vorticity(u850, v850)


# In[8]:


fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111, projection=ccrs.LambertConformal())
ax.set_extent([-130, -65, 20, 55])
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=1)
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor='black')
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=1, edgecolor='black')
plt1 = ax.contourf(t850['longitude'], t850['latitude'],t850.data.to('degF'),
levels=np.linspace(-15, 75, 25),cmap='nipy_spectral', transform=ccrs.PlateCarree())


# In[9]:


plt2 = ax.contour(z850['longitude'], z850['latitude'], z850.data.to('decameter'),
levels=np.linspace(130, 160, 11),colors='black', transform=ccrs.PlateCarree())
ax.clabel(plt2, np.linspace(130, 160, 6), inline=True, fmt='%d', fontsize=14)

cb = fig.colorbar(plt1, ax=ax, orientation='horizontal', pad=0.05, aspect=30)
cb.set_ticks(np.linspace(5, 50, 10))


# In[28]:


plt3 = ax.barbs(u850['longitude'], u850['latitude'], u850.data.to('knot'),v850.data.to('knot'), flagcolor='r', barbcolor='k', flip_barb=False, transform=ccrs.PlateCarree(),barb_increments=dict(half=10, full=20, flag=100),length=5, sizes={'spacing':1},pivot='middle')
#, transform=ccrs.PlateCarree())


#, flagcolor='r', barbcolor='k', flip_barb=False, transform=ccrs.PlateCarree(),
#barb_increments=dict(half=10, full=20, flag=100),length=5, sizes={'spacing':1},pivot='middle')


# In[29]:


gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, dms=True, x_inline=False,
y_inline=False, linewidth=1, color='black', linestyle=':')
gl.xlocator = mticker.FixedLocator([-120, -110, -100, -90, -80, -70])
gl.ylocator = mticker.FixedLocator([25, 30, 35, 40, 45])
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 16, 'rotation': 20}
gl.ylabel_style = {'size': 16}


# In[30]:


plt.show()


# In[ ]:




