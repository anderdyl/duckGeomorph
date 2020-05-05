

import xarray as xr
import matplotlib.pyplot as plt

inmcm4_rhsmin_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_macav2metdata_rhsmin_inmcm4_r1i1p1_rcp45_2006_2099_CONUS_daily.nc'
inmcm4_rhsmax_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_macav2metdata_rhsmax_inmcm4_r1i1p1_rcp45_2006_2099_CONUS_daily.nc'

inmcm4_rhsmin = xr.open_dataset(inmcm4_rhsmin_url)
inmcm4_rhsmax = xr.open_dataset(inmcm4_rhsmax_url)


test = inmcm4_rhsmax.isel(time=0)
rh = test['relative_humidity']
lat = test['lat']
lon = test['lon']

plt.pcolor(lon,lat,rh)
