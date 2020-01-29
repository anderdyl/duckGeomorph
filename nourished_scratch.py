#   Lets identify a transect to the north to check on how much nourishment evolution was caught
#
from getdatatestbed import getDataFRF
import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from glob import glob
import os

geomorphdir = '/media/dylananderson/Elements/FRF_Geomorph/'

files = os.listdir(geomorphdir)

files.sort()

files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]

def getBathy(file, lower, upper):
    bathy = Dataset(file)

    xs_bathy = bathy.variables['xFRF'][:]
    ys_bathy = bathy.variables['yFRF'][:]
    zs_bathy = bathy.variables['elevation'][:]
    ts_bathy = bathy.variables['time'][:]
    pr_bathy = bathy.variables['profileNumber'][:]

    zs_bathy = np.ma.masked_where((pr_bathy > upper), zs_bathy)
    ys_bathy = np.ma.masked_where((pr_bathy > upper), ys_bathy)
    xs_bathy = np.ma.masked_where((pr_bathy > upper), xs_bathy)
    pr_bathy = np.ma.masked_where((pr_bathy > upper), pr_bathy)
    ts_bathy = np.ma.masked_where((pr_bathy > upper), ts_bathy)

    zs_bathy = np.ma.masked_where((pr_bathy < lower), zs_bathy)
    ys_bathy = np.ma.masked_where((pr_bathy < lower), ys_bathy)
    xs_bathy = np.ma.masked_where((pr_bathy < lower), xs_bathy)
    pr_bathy = np.ma.masked_where((pr_bathy < lower), pr_bathy)
    ts_bathy = np.ma.masked_where((pr_bathy < lower), ts_bathy)

    output = dict()
    output['x'] = xs_bathy
    output['y'] = ys_bathy
    output['z'] = zs_bathy
    output['pr'] = pr_bathy
    output['t'] = ts_bathy

    return output


time = []
for i in range(len(files)):
    temp = files[i].split('_')
    time = np.append(time, DT.datetime.strptime(temp[1], '%Y%m%d'))


# Lets find just the files after January 2017
recent = np.nonzero((time > DT.datetime(2017,1,1)) & (time < DT.datetime(2018,1,1)))
subset = files[recent[0][0]:recent[0][-1]].copy()


colormap = plt.cm.gist_ncar

fig, ax = plt.subplots(2,1)

ax[0].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subset)))))
labels = []
for i in range(len(subset)):

    data = getBathy(os.path.join(geomorphdir, subset[i]), lower=1070, upper=1100)
    temp = subset[i].split('_')
    ax[0].plot(data['x'], data['z'], label=temp[1])

ax[0].legend(loc='upper right')
ax[0].set_xlim([50, 800])
ax[0].set_ylim([-8, 4])
ax[0].set_title('Profile variability in 2017')

# Lets find just the files after January 2018
recent = np.nonzero((time > DT.datetime(2018, 1, 1)) & (time < DT.datetime(2019, 1, 1)))
subset1 = files[recent[0][0]:recent[0][-1]].copy()


colormap = plt.cm.gist_ncar

ax[1].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subset1)))))
labels = []
for i in range(len(subset1)):

    data = getBathy(os.path.join(geomorphdir, subset1[i]), lower=1070, upper=1100)
    temp = subset1[i].split('_')
    ax[1].plot(data['x'], data['z'], label=temp[1])

ax[1].legend(loc='upper right')
ax[1].set_xlim([50, 800])
ax[1].set_ylim([-8, 4])
ax[1].set_title('Profile variability in 2018')




