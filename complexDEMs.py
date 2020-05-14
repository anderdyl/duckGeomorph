#   Lets load a couple DEMs and get a work flow running for 2D complex EOFs
#
from getdatatestbed import getDataFRF
import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from glob import glob
import os
from scipy.interpolate import interp1d
import sandBarTool.morphLib as mL

#from sandBarTool import morphLib
#from downloads import pyeemd



geomorphdir = '/media/dylananderson/Elements/FRF_DEMs/'
files = os.listdir(geomorphdir)
files.sort()
files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]

def getBathy(file):
    bathy = Dataset(file)

    xs_bathy = bathy.variables['xFRF'][:]
    ys_bathy = bathy.variables['yFRF'][:]
    zs_bathy = bathy.variables['elevation'][:]
    ts_bathy = bathy.variables['time'][:]
    #pr_bathy = bathy.variables['profileNumber'][:]

    # zs_bathy = np.ma.masked_where((pr_bathy > upper), zs_bathy)
    # ys_bathy = np.ma.masked_where((pr_bathy > upper), ys_bathy)
    # xs_bathy = np.ma.masked_where((pr_bathy > upper), xs_bathy)
    # pr_bathy = np.ma.masked_where((pr_bathy > upper), pr_bathy)
    # ts_bathy = np.ma.masked_where((pr_bathy > upper), ts_bathy)
    #
    # zs_bathy = np.ma.masked_where((pr_bathy < lower), zs_bathy)
    # ys_bathy = np.ma.masked_where((pr_bathy < lower), ys_bathy)
    # xs_bathy = np.ma.masked_where((pr_bathy < lower), xs_bathy)
    # pr_bathy = np.ma.masked_where((pr_bathy < lower), pr_bathy)
    # ts_bathy = np.ma.masked_where((pr_bathy < lower), ts_bathy)

    output = dict()
    output['x'] = xs_bathy
    output['y'] = ys_bathy
    output['z'] = zs_bathy
    # output['pr'] = pr_bathy
    output['t'] = ts_bathy

    return output


getFile = os.path.join(geomorphdir, files[-150])
data = getBathy(getFile)

plt.figure()
plt.pcolor(data['x'],data['y'],data['z'][0,:,:])
