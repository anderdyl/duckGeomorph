from PIL import Image
import numpy
from osgeo import gdal
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import elevation
# im = Image.open('USGS_13_n37w076.tiff')
# imarray = numpy.array(im)
# imarray.shape
# im.show()

import rasterio
# from matplotlib import pyplot as plt
# from rasterio.plot import show
# from rasterio.plot import show_hist
#
#
# with rasterio.open('USGS_13_n37w076.tiff') as src:
#     data = src.read(1)
#     print(src.width, src.height)
#     print(src.crs)
#     print(src.transform)
#     print(src.count)
#     print(src.indexes)
#     #print(src.window(**src.window_bounds(((100, 200), (100, 200)))))
# # show((src, 1), cmap='viridis')# type(raster)
#
# # show_hist(src) #, bins=50, lw=0.0, stacked=False, alpha=0.3,histtype='stepfilled', title="Histogram")
#
# plt.imshow(data)
# # plt.show()
#
# filename = "USGS_13_n37w076.tiff"
# gdal_data = gdal.Open(filename)
# gdal_band = gdal_data.GetRasterBand(1)
# nodataval = gdal_band.GetNoDataValue()
#
# # convert to a numpy array
# data_array = gdal_data.ReadAsArray().astype(np.float)
#
#
# # replace missing values if necessary
# if np.any(data_array == nodataval):
#     data_array[data_array == nodataval] = np.nan
#
# #Plot our data with Matplotlib's 'contourf'
# fig = plt.figure(figsize = (12, 8))
# ax = fig.add_subplot(111)
# plt.contourf(data_array, cmap = "viridis",
#             levels = list(range(-10, 20, 1)))
# plt.title("Elevation Contours of Mt. Shasta")
# cbar = plt.colorbar()
# plt.gca().set_aspect('equal', adjustable='box')
# plt.show()

from netCDF4 import Dataset
largeMap = Dataset('crm_vol2.nc')
x = largeMap.variables['x'][:]
y = largeMap.variables['y'][:]
z = largeMap.variables['z'][:]

import cmocean
cmap = cmocean.cm.topo
newcmap = cmocean.tools.crop_by_percent(cmap, 10, which='both', N=None)
fig = plt.figure(figsize = (12, 8))
ax = plt.subplot2grid((4,4),(0,0),rowspan=4,colspan=2)

mappable = ax.contourf(x, y, z, cmap = newcmap, levels = np.arange(-60,61,0.5)) #list(range(-60, 60, 0.5))
ax.set_ylim([34.8,37])
ax.set_xlim([-76.5, -74.8])

# plt.contourf(largeMap.variables['x'][:],largeMap.variables['y'][:],largeMap.variables['z'][:], cmap = "viridis",
#             levels = list(range(-50, 100, 5)))
# plt.title("Elevat")
# cbar = plt.colorbar(mappable, ax=ax)
plt.gca().set_aspect('equal', adjustable='box')
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
axin = inset_axes(ax, width="40%", height="40%", loc=4)
# Global inset map.
inmap = Basemap(projection='ortho', lon_0=-77.5, lat_0=35.5,
                ax=axin, anchor='NE')
inmap.drawcountries(color='white')
inmap.fillcontinents(color='gray')
# bx, by = inmap(m.boundarylons, m.boundarylats)
# xy = list(zip(bx, by))
coord = [[-78, 34], [-74.5, 34], [-74.5, 37], [-78, 37], [-78, 34]]
xs, ys = zip(*coord)
# mapboundary = Polygon(coord, edgecolor='k', linewidth=1, fill=False)
# inmap.ax.plot(coord) #add_patch(mapboundary)
axin.plot(np.array(xs),np.array(ys),'k-') #add_patch(mapboundary)


jp2s = ['m_3607550_se_18_1_20141005_20141118.jp2','m_3607551_sw_18_1_20141005_20141118.jp2','m_3607550_ne_18_1_20141005_20141118.jp2']
jp21 = 'm_3607551_sw_18_1_20141005_20141118.jp2'
jp22 = 'm_3607550_se_18_1_20141005_20141118.jp2'
jp23 = 'm_3607550_ne_18_1_20141005_20141118.jp2'

with rasterio.open(jp21) as f:
    im1r = f.read(1)
    im1g = f.read(2)
    im1b = f.read(3)
with rasterio.open(jp22) as f:
    im2r = f.read(1)
    im2g = f.read(2)
    im2b = f.read(3)
with rasterio.open(jp23) as f:
    im3 = f.read(1)

dataIM1r = np.array(im1r, dtype=im1r[0].dtype)
dataIM1g = np.array(im1g, dtype=im1g[0].dtype)
dataIM1b = np.array(im1b, dtype=im1b[0].dtype)
dataIM2r = np.array(im2r, dtype=im2r[0].dtype)
dataIM2g = np.array(im2g, dtype=im2g[0].dtype)
dataIM2b = np.array(im2b, dtype=im2b[0].dtype)

# dataIM3 = np.array(im3, dtype=im1[0].dtype)

mergedr = np.hstack((dataIM2r,dataIM1r))
mergedg = np.hstack((dataIM2g,dataIM1g))
mergedb = np.hstack((dataIM2b,dataIM1b))
n, m = np.shape(mergedb)
import cv2

rgbArray = np.zeros((n, m, 3), 'uint8')
rgbArray[..., 0] = mergedr  # imR[:, :, 0] * 255
rgbArray[..., 1] = mergedg  # imG[:, :, 0] * 255
rgbArray[..., 2] = mergedb  # imB[:, :, 0] * 255
imGrayscale = cv2.cvtColor(rgbArray, cv2.COLOR_BGR2GRAY)
ax2 = plt.subplot2grid((4,4),(0,2),rowspan=3,colspan=2)
# ax2.imshow(imGrayscale)
ax2.imshow(rgbArray)

ax2.set_xlim([5500,8500])
ax2.set_ylim([2000,0])

# arrs = []
# for jp2 in jp2s:
#     with rasterio.open(jp2) as f:
#         arrs.append(f.read(1))
#
# data = np.array(arrs, dtype=arrs[0].dtype)
# data

plt.show()



