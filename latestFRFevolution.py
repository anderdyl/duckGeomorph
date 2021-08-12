#   Lets identify a transect to the north to check on how much nourishment evolution was caught
#
from getdatatestbed import getDataFRF
import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from glob import glob
import os

geomorphdir = '/media/dylananderson/Elements/filteredFRF_Geomorph/'

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
recent = np.nonzero((time > DT.datetime(2017,7,22)) & (time < DT.datetime(2019,4,1)))
# recent = np.nonzero((time > DT.datetime(2018,5,22)) & (time < DT.datetime(2018,9,1)))

subset = files[recent[0][0]:recent[0][-1]].copy()


colormap = plt.cm.gist_ncar
plt.style.use('default')
fig, ax = plt.subplots(2,1)

ax[0].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subset)))))
labels = []
for i in range(len(subset)):
    data = getBathy(os.path.join(geomorphdir, subset[i]), lower=1070, upper=1100)

    # data = getBathy(os.path.join(geomorphdir, subset[i]), lower=700, upper=800)
    temp = subset[i].split('_')
    ax[0].plot(data['x'], data['z'], label=temp[1])

ax[0].legend(loc='upper right')
ax[0].set_xlim([50, 800])
ax[0].set_ylim([-8, 4])
ax[0].set_title('Profile variability in northern edge of property')

# Lets find just the files after January 2018
recent = np.nonzero((time > DT.datetime(2020, 5, 1)) & (time < DT.datetime(2020, 12, 1)))
subset1 = files[recent[0][0]:].copy()


colormap = plt.cm.gist_ncar

ax[1].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subset1)))))
labels = []
for i in range(len(subset1)):

    data = getBathy(os.path.join(geomorphdir, subset1[i]), lower=720, upper=750)
    temp = subset1[i].split('_')
    ax[1].plot(data['x'], data['z'], label=temp[1])

ax[1].legend(loc='upper right')
ax[1].set_xlim([50, 800])
ax[1].set_ylim([-8, 4])
ax[1].set_title('Profile variability in 2018')





# plt.rcParams.update({
#     "lines.color": "white",
#     "patch.edgecolor": "white",
#     "text.color": "white",
#     "axes.facecolor": "black",
#     "axes.edgecolor": "lightgray",
#     "axes.labelcolor": "white",
#     "xtick.color": "white",
#     "ytick.color": "white",
#     "grid.color": "lightgray",
#     "figure.facecolor": "black",
#     "figure.edgecolor": "black",
#     "savefig.facecolor": "black",
#     "savefig.edgecolor": "black"})


plt.style.use('dark_background')


fig = plt.figure()
ax1 = plt.subplot2grid((1,1),(0,0),colspan=1,rowspan=1)
i = -1
data = getBathy(os.path.join(geomorphdir, subset1[i]), lower=720, upper=750)
temp = subset1[i].split('_')
legendLabel = temp[1][0:4] + '-' + temp[1][4:6] + '-' + temp[1][6:9]
ax1.plot(data['x'], data['z'], label=legendLabel,color='orange')
ax1.legend(loc='upper right')
ax1.set_xlim([50, 600])
ax1.set_ylim([-6, 4])
ax1.set_xlabel('cross-shore (m)')
ax1.set_ylabel('elevation (m)')
#ax1.set_title('Profile variability in 2018')




#recent2 = np.nonzero((time > DT.datetime(2016, 6, 1)) & (time < DT.datetime(2018, 1, 1)))
recent2 = np.nonzero((time > DT.datetime(2020, 5, 1)) & (time < DT.datetime(2020, 12, 1)))

subset2 = files[recent2[0][0]:].copy()
del subset2[37:]
fig2 = plt.figure()
ax2 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
ax2.set_prop_cycle(plt.cycler('color', plt.cm.rainbow(np.linspace(0, 1, len(subset2)))))

labels = []
for i in range(len(subset2)):

    data = getBathy(os.path.join(geomorphdir, subset2[i]), lower=720, upper=750)
    temp = subset2[i].split('_')
    ax2.plot(data['x'], data['z'], label=temp[1])

ax2.legend(loc='upper right')
ax2.set_xlim([50, 600])
ax2.set_ylim([-6, 4])
ax2.set_xlabel('cross-shore (m)')
ax2.set_ylabel('elevation (m)')
ax2.set_title('Profile Variability since May 2020')
