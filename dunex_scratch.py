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


subset = files[41:49].copy()
#from palettable.colorbrewer.sequential import Blues_8
#ax.imshow(data, cmap=Blues_8.mpl_colormap

colormap = plt.cm.gist_ncar

fig, ax = plt.subplots(2,1)

ax[0].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subset)))))
labels = []
for i in range(len(subset)):

    data = getBathy(os.path.join(geomorphdir, subset[i]), lower=270, upper=280)
    temp = subset[i].split('_')
    ax[0].plot(data['x'], data['z'], label=temp[1])

ax[0].legend(loc='upper right')
ax[0].set_xlim([50, 800])
ax[0].set_ylim([-8, 4])
ax[0].set_title('South of the pier: yFRF=275')

subset1 = files[41:49].copy()

colormap = plt.cm.gist_ncar

ax[1].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subset1)))))
labels = []
for i in range(len(subset1)):

    #data = getBathy(os.path.join(geomorphdir, subset1[i]), lower=1070, upper=1100)
    data = getBathy(os.path.join(geomorphdir, subset1[i]), lower=760, upper=780)

    temp = subset1[i].split('_')
    ax[1].plot(data['x'], data['z'], label=temp[1])
    #labels.append(temp[1])

ax[1].legend(loc='upper right')
ax[1].set_xlim([50, 800])
ax[1].set_ylim([-8, 4])
ax[1].set_title('North of the pier: yFRF=775')


#
# subset1 = files[15:47].copy()
#
# colormap = plt.cm.gist_ncar
#
# ax[2].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subset1)))))
# labels = []
# for i in range(len(subset1)):
#
#     data = getBathy(os.path.join(geomorphdir, subset1[i]), lower=1070, upper=1100)
#     temp = subset1[i].split('_')
#     ax[2].plot(data['x'], data['z'], label=temp[1])
#
# ax[2].legend(loc='upper right')
# ax[2].set_xlim([50, 800])
# ax[2].set_ylim([-8, 4])

#loc='upper center',
#           bbox_to_anchor=[0.5, 1.1],
#           columnspacing=1.0, labelspacing=0.0,
#           handletextpad=0.0, handlelength=1.5,
#           fancybox=True, shadow=True)


# go  = getDataFRF.getObs(DT.datetime(1900,1,1), DT.datetime(2019,12,1))
# survey = go.getBathyTransectFromNC(forceReturnAll=True)
# print(survey.keys())
# # what are my north and south lines of interest
# NorthIdx = survey['profileNumber'] == 1097
# southIdx = survey['profileNumber'] == 1
# crossShoreMax = 1500   # how far do we want to look in cross-shore

# startTime = '2019-09-15T00:00:00Z'           # project start time
# endTime = '2019-11-08T00:00:00Z'             # project End time
#
# d1 = DT.datetime.strptime(startTime, '%Y-%m-%dT%H:%M:%SZ')
# d2 = DT.datetime.strptime(endTime, '%Y-%m-%dT%H:%M:%SZ')
# server = 'CHL'
#
# gdTB = getDataFRF.getDataTestBed(d1, d2, THREDDS=server)
# gObs = getDataFRF.getObs(d1, d2, THREDDS=server)
#
# bathy = gObs.getBathyTransectFromNC(forceReturnAll=True)  # , ForcedSurveyDate=ForcedSurveyDate)

# x = bathy['xFRF']
# y = bathy['yFRF']
# z = bathy['profileNumber']
# label = 'survey'
# title = 'survey'
# fig3a, ax3a = plt.subplots(1, 1, figsize=(5, 5))
# sc3a = ax3a.scatter(x, y, c=z, vmin=0, vmax = 2000)
# cbar = plt.colorbar(sc3a, ax=ax3a)
# cbar.set_label(label)
# ax3a.set_title(title)












