# Upstate New York Snowstorm Case Study

This repository contains code and analysis examining the meteorological factors leading to a snowstorm event across upstate New York on February 25, 2022. 

## Overview

On February 25, 2022, a widespread 6-12 inch snowfall occurred across upstate New York driven by the development of a coastal cyclone off the East Coast. This repository examines the synoptic patterns as well as mesoscale processes like frontogenesis that contributed to the heavy snowfall rates observed.

## Contents

- `report.pdf`: A full report analyzing the evolution of the surface cyclone, fronts, and vertical atmospheric structure associated with the snowstorm. Includes a discussion of frontogenesis and its role in snowband development.

- `analysis.py`: Python script using MetPy and Cartopy to visualize conditions at 850 hPa including temperatures, geopotential heights, and winds leading up to and during the snowstorm event.

- `frontogenesis.py`: Python code to calculate frontogenetical forcing and differential vorticity advection using model analysis data. 

## Usage

The Python scripts require the following packages:

- xarray
- MetPy
- Cartopy
- matplotlib
- scipy

The data used in `analysis.py` is read in from a sample NetCDF file. The frontogenesis calculations in `frontogenesis.py` can be applied to model analysis data in NetCDF format.

## Resources

Relevant data sources:

- [ERA5 Reanalysis](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5)
- [NCEP GFS Analysis](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs)

## Acknowledgements

This analysis was completed as as final project for ATMS 502 at University of Washington.