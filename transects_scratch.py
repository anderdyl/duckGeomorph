#   Lets identify a transect to the north to check on how much nourishment evolution was caught
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


def interpBathy(xSub, zSub, x):

    f = interp1d(xSub, zSub, kind='linear', bounds_error=False)
    newz = f(x)
    newBathy = dict()
    newBathy['x'] = x
    newBathy['z'] = newz
    return newBathy


subset = files[0:969].copy()
#from palettable.colorbrewer.sequential import Blues_8
#ax.imshow(data, cmap=Blues_8.mpl_colormap

colormap = plt.cm.gist_ncar


labels = []
xinterp = np.arange(108, 608, 2.5)

bathy = dict()
#alllines = np.empty((len(xinterp),))
count = 0
count2 = 0
worstcount = 0
#fig = plt.figure(figsize=(10,10))
for i in range(len(subset)):
    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=1080, upper=1100)
    data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-2, upper=10)

    temp = subset[i].split('_')

    surveydate = DT.datetime.strptime(temp[1], '%Y%m%d')
    elevs = data['z']
    cross = data['x']
    xSub = np.ma.MaskedArray.filled(cross, np.nan)
    zSub = np.ma.MaskedArray.filled(elevs, np.nan)
    realValues = ~np.isnan(xSub)

    if any(realValues):

        newdata = interpBathy(xSub, zSub, xinterp)
        nanValues = np.isnan(newdata['z'])

        if any(nanValues):
            #print('Found a transect with nans {}'.format(i))
            print('Trying to extend the lower half of the profile with a linear fit: {}'.format(surveydate))
            if count2 == 0:
                badlines = newdata['z']
                count2 = count2+1
                badtime = surveydate
            else:
                badlines = np.vstack((badlines, newdata['z']))
                count2 =count2+1
                badtime = np.append(badtime,surveydate)

            #print('Threw out a survey with data due to nans: {}'.format(surveydate))
            #plt.plot(xSub,zSub)


            xindex = np.where((xSub>=650))

            if len(xindex[0]) > 3:

                xbot = xSub[xindex]
                zbot = zSub[xindex]
                mb = np.polyfit(xbot, zbot, 1)
                f = np.poly1d(mb)
                maxXind = np.where((xinterp > np.nanmax(xSub)))
                newX = xinterp[maxXind]
                newZ = f(newX)
                xSubNew = np.append(xSub, newX)
                zSubNew = np.append(zSub, newZ)
                moredata = interpBathy(xSubNew, zSubNew, xinterp)

                del xSubNew, zSubNew,newZ,newX,f,mb,zbot,xbot

                nanValues2 = np.isnan(moredata['z'])

                if any(nanValues2):
                    print('WHAT HAPPENED AT TRANSECT {}'.format(i))
                    #plt.plot(moredata['x'], moredata['z'], 'r-')
                else:
                    if count == 0:
                        alllines = moredata['z']
                        time = surveydate
                        count = count + 1
                    else:
                        alllines = np.vstack((alllines, moredata['z']))
                        time = np.append(time, surveydate)
                        count = count + 1
                    #plt.plot(moredata['x'], moredata['z'], 'k-')
                    del moredata



            elif len(xindex[0]) == 3:
                print('Could not do that: number of points being fit: {}'.format(len(xindex[0])))
                worstcount = worstcount+1
            else:
                print('Could not do that: maximum X value in transect is: {}'.format(np.nanmax(xSub)))
                worstcount = worstcount+1
                #plt.plot(xSub, zSub, 'r-')


        else:
            if count == 0:
                alllines = newdata['z']
                time = surveydate
                count = count+1
            else:
                alllines = np.vstack((alllines, newdata['z']))
                time = np.append(time,surveydate)
                count = count+1
    else:
        print('Survey with no data at this line {}'.format(i))



fig, ax = plt.subplots(2,1)
ax[0].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(alllines)))))
for i in range(len(alllines)):
    ax[0].plot(xinterp, alllines[i,:], label=time[i])

#ax[0].legend(loc='upper right')
ax[0].set_xlim([50, 850])
ax[0].set_ylim([-8, 4])
ax[0].set_title('Cross-shore profile variability')

tg, xg = np.meshgrid(time, xinterp)
ax[1].pcolor(xg,tg,alllines.T)


fig20 = plt.figure(figsize=(10,5))
plt.plot(time,alllines[:,-1])
plt.plot(time,alllines[:,-10])
plt.plot(time,alllines[:,-20])
plt.plot(time,alllines[:,-30])




from sklearn import preprocessing
scalerMorph = preprocessing.StandardScaler().fit(alllines)
X = preprocessing.scale(alllines)
#
# #from sklearn import joblib
# #jolib.dum

#X = scalerMorph.transform(alllines)


from sklearn.decomposition import PCA

skpca = PCA()
skpca.fit(X)

f, ax = plt.subplots(figsize=(5,5))
ax.plot(skpca.explained_variance_ratio_[0:10]*100)
ax.plot(skpca.explained_variance_ratio_[0:10]*100,'ro')
ax.set_title("% of variance explained", fontsize=14)
ax.grid()

PCs = skpca.transform(X)
EOFs = skpca.components_
variance = skpca.explained_variance_
n_components = np.shape(PCs)[0]
n_features = np.shape(EOFs)[1]
pred_mean = alllines.mean(axis=0)
pred_std = alllines.std(axis=0)
scalerPCs = PCs/np.sqrt(variance[0])

import xarray as xr

morphPCA = xr.Dataset(
        {
            'PCs': (('time', 'n_components'), PCs),
            'EOFs': (('n_components','n_features'), EOFs),
            'variance': (('n_components',), variance),

            'pred_mean': (('n_features',), pred_mean),
            'pred_std': (('n_features',), pred_std),
            'pred_time': (('time',), time),

        })


from morph_utils import KMA_simple
numClusters = 20

KMA = KMA_simple(morphPCA, numClusters, repres=1)

from morph_utils import sort_cluster_gen_corr_end as sortClusters

sorted = sortClusters(KMA.centroids, numClusters)

from matplotlib import gridspec
fig2 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(4, 5, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
for i in range(numClusters):
    #getind = sorted[i]
    profile = KMA.centroids[i,:] * pred_std + pred_mean

    ax = plt.subplot(gs[gr, gc])

    ax.plot(xinterp,profile)
    ax.set_xlim([80, 820])
    ax.set_ylim([-8, 2.25])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    ax.text(400,0, KMA.group_size.values[i], fontweight='bold')

    if gc > 0:
        ax.set_yticks([])

    if gr < (5-1):
        ax.set_xticks([])
    #  counter
    gc += 1
    if gc >= 5:
        gc = 0
        gr += 1

import scipy.signal as sig
profiles = np.zeros((np.shape(KMA.centroids)))
deviation = np.zeros((np.shape(KMA.centroids)))

for i in range(numClusters):
    profiles[i,:] = KMA.centroids[i,:] * pred_std + pred_mean
    deviation[i,:] = profiles[i,:] - pred_mean


offshorePeaks = np.zeros((numClusters,))
inshorePeaks = np.zeros((numClusters,))
inshoreDepth = np.zeros((numClusters,))
offshoreDepth = np.zeros((numClusters,))


fig3 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(4, 5, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
for i in range(numClusters):
    #getind = sorted[i]
    dev = KMA.centroids[i,:] * pred_std
    true = KMA.centroids[i,:] * pred_std + pred_mean

    peaks = sig.find_peaks(x=(true), prominence=0.05)

    if len(peaks[0]) > 0:
        offshorePeaks[i] = np.max(peaks[0])
        offshoreDepth[i] = true[np.max(peaks[0])]

    if len(peaks[0]) > 1:
        inshorePeaks[i] = np.min(peaks[0])
        inshoreDepth[i] = true[np.min(peaks[0])]

    ax = plt.subplot(gs[gr, gc])

    ax.plot(xinterp,true)
    ax.plot(xinterp[peaks[0]],true[peaks[0]],'ro')
    ax.set_xlim([80, 820])
    ax.set_ylim([-8, 2.25])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    ax.text(400,0, KMA.group_size.values[i], fontweight='bold')

    if gc > 0:
        ax.set_yticks([])

    if gr < (5-1):
        ax.set_xticks([])
    #  counter
    gc += 1
    if gc >= 5:
        gc = 0
        gr += 1



# totalDepths = offshoreDepth+inshoreDepth
sortedPeaks = np.sort(offshorePeaks)
sortedPeakInd = np.argsort(offshorePeaks)

# doubleIndex = np.where((inshorePeaks > 0))
# singleIndex = np.where((inshorePeaks == 0))
# sortedInPeaks = np.argsort(offshorePeaks[doubleIndex[0]])
# sortedOffPeaks = np.argsort(offshorePeaks[singleIndex[0]])
# sortedPeakInd = np.append(singleIndex[0][sortedOffPeaks],doubleIndex[0][sortedInPeaks])


fig3 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(10, 2, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, numClusters))

for i in range(numClusters):
    #getind = sorted[i]
    dev = KMA.centroids[sortedPeakInd[i],:] * pred_std
    true = KMA.centroids[sortedPeakInd[i],:] * pred_std + pred_mean

    peaks = sig.find_peaks(x=(true), prominence=0.05)

    #if len(peaks[0]) > 0:
    #    offshorePeaks[i] = np.max(peaks[0])

    ax = plt.subplot(gs[gr, gc])

    ax.plot(xinterp,true,color=colors[i])
    #ax.plot(xinterp[peaks[0]],true[peaks[0]],'ro')
    ax.set_xlim([80, 720])
    ax.set_ylim([-8, 2.25])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    ax.text(400,0, KMA.group_size.values[i], fontweight='bold')

    if gc > 0:
        ax.set_yticks([])

    if gr < (10-1):
        ax.set_xticks([])
    #  counter
    gr += 1
    if gr >= 10:
        gc += 1
        gr = 0


bmu = np.zeros((len(time),))

for i in range(numClusters):
    bmuIndex = np.where((KMA.bmus == sortedPeakInd[i]))
    bmu[bmuIndex] = i

fig5 = plt.figure(figsize=(10,5))
plt.plot(bmu,time,'o')




fig6 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(4, 5, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
for i in range(numClusters):
    #getind = sorted[i]
    bmuIndex = np.where((KMA.bmus == sortedPeakInd[i]))

    subset = alllines[bmuIndex,:]

    ax = plt.subplot(gs[gr, gc])
    for x in range(len(bmuIndex[0])):
        ax.plot(xinterp,subset[0,x,:])

    ax.set_xlim([80, 820])
    ax.set_ylim([-8, 2.25])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    ax.text(400,0, KMA.group_size.values[sortedPeakInd[i]], fontweight='bold')

    if gc > 0:
        ax.set_yticks([])

    if gr < (5-1):
        ax.set_xticks([])
    #  counter
    gc += 1
    if gc >= 5:
        gc = 0
        gr += 1




innerBar = np.nan * np.ones((alllines.shape[0],))
outerBar = np.nan * np.ones((alllines.shape[0],))

zOuterBar = np.nan * np.ones((alllines.shape[0],))
zInnerBar = np.nan * np.ones((alllines.shape[0],))

innerTrough = np.nan * np.ones((alllines.shape[0],))
outerTrough = np.nan * np.ones((alllines.shape[0],))
zOuterTrough = np.nan * np.ones((alllines.shape[0],))
zInnerTrough = np.nan * np.ones((alllines.shape[0],))


for tt in range(alllines.shape[0]):
    bathy = alllines[tt,:]
    bathyX = xinterp
    if bathy[~np.isnan(bathy)].any():
        fname = "/home/dylananderson/projects/duckGeomorph/sandBarTool/BarId_{}.png".format(time[tt].strftime("%Y%m%d"))
        #xFRFbar, xFRFtrough = mL.findSandBarAndTrough1D(bathyX, bathy, plotFname=fname, smoothLengthScale=10, profileTrend=np.mean(alllines,axis=0))
        xFRFbar, xFRFtrough = mL.findSandBarAndTrough1D(bathyX, bathy, plotFname=None, smoothLengthScale=20, profileTrend=np.mean(alllines,axis=0))

        if xFRFbar is not None:
            #for sandbarX in xFRFbar:
            #    ax[1].plot(sandbarX, time[tt], 'ro', label='bar')

            if len(xFRFbar) > 1:
                outerBar[tt] = np.max(xFRFbar)
                zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == bathyX))
                zOuterBar[tt] = bathy[zOuterInd[0]]

                innerBar[tt] = np.min(xFRFbar)
                zInnerInd = np.where((np.round(2*np.min(xFRFbar))/2 == bathyX))
                zInnerBar[tt] = bathy[zInnerInd[0]]

            else:
                outerBar[tt] = xFRFbar
                zOuterInd = np.where((np.round(2*xFRFbar)/2 == bathyX))
                zOuterBar[tt] = bathy[zOuterInd[0]]

        if xFRFtrough is not None:
            #for troughX in xFRFtrough:
            #    ax[1].plot(troughX, time[tt], 'bd', label='trough')
            if len(xFRFtrough) > 1:
                outerTrough[tt] = np.max(xFRFtrough)
                zOuterInd = np.where((np.round(2*np.max(xFRFtrough))/2 == bathyX))
                zOuterTrough[tt] = bathy[zOuterInd[0]]

                innerTrough[tt] = np.min(xFRFtrough)
                zInnerInd = np.where((np.round(2*np.min(xFRFtrough))/2 == bathyX))
                zInnerTrough[tt] = bathy[zInnerInd[0]]
            else:
                outerTrough[tt] = xFRFtrough
                zOuterInd = np.where((np.round(2*xFRFtrough)/2 == bathyX))
                zOuterTrough[tt] = bathy[zOuterInd[0]]
        #barPatch = mpatches.Patch(color='red', label='sandbar')
        #troughPatch = mpatches.Patch(color='blue', label='trough')


fig2 = plt.figure(figsize=(8,8))
plt.plot(time,innerBar,'bo')
plt.plot(time,outerBar,'ro')
plt.show()


fig11, ax11 = plt.subplots(2,1)
#ax10 = fig.add_subplot(111)
#ax10.set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, numClusters))))
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, numClusters))
colors2 = cm.gray(np.linspace(0, 1, numClusters))
for i in range(numClusters):
    bmuIndex = np.where((KMA.bmus == sortedPeakInd[i]))
    ax11[0].scatter(time[bmuIndex], outerBar[bmuIndex], label=time[i], color = colors[i])
    ax11[1].scatter(time[bmuIndex], innerBar[bmuIndex], color = colors[i])

ax11[0].set_ylabel('xFRF')
ax11[1].set_ylabel('xFRF')

plt.show()



fig3 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(10, 10, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, numClusters))

for i in range(numClusters):
    #getind = sorted[i]
    dev = KMA.centroids[sortedPeakInd[i],:] * pred_std
    true = KMA.centroids[sortedPeakInd[i],:] * pred_std + pred_mean

    peaks = sig.find_peaks(x=(true), prominence=0.05)

    #if len(peaks[0]) > 0:
    #    offshorePeaks[i] = np.max(peaks[0])

    ax = plt.subplot(gs[gr, gc])

    ax.plot(xinterp,true,color=colors[i])
    #ax.plot(xinterp[peaks[0]],true[peaks[0]],'ro')
    ax.set_xlim([80, 720])
    ax.set_ylim([-8, 2.25])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    ax.text(400,0, KMA.group_size.values[i], fontweight='bold')

    if gc > 0:
        ax.set_yticks([])

    if gr < (10-1):
        ax.set_xticks([])
    #  counter
    gr += 1
    if gr >= 10:
        gc += 1
        gr = 0


ax15 = plt.subplot(gs[0:5, 3:])
ax16 = plt.subplot(gs[6:,3:])
for i in range(numClusters):
    bmuIndex = np.where((KMA.bmus == sortedPeakInd[i]))
    ax15.scatter(time[bmuIndex], outerBar[bmuIndex], label=time[i], color = colors[i])
    ax16.scatter(time[bmuIndex], innerBar[bmuIndex], color = colors[i])

ax15.set_ylabel('xFRF')
ax16.set_ylabel('xFRF')
ax15.set_title('Offshore Bar')
ax16.set_title('Onshore Bar')

import xarray as xr

encoding = {'alllines':{'dtype':'float32','_FillValue':-9999}}

da_argus6 = xr.Dataset({'alllines': (['t','x'], alllines)},
                           coords={'t': time,
                                   'x': xinterp})

da_argus6.to_netcdf('alllinesSouth.nc',encoding=encoding)




#PCs_std = scaler_PCs.transform(PCs)

#scaler_alllines = preprocessing.scale(alllines)



# from scipy.signal import detrend
#
# f, ax = plt.subplots(figsize=(10,4))
# PCdf.ix[:,0].plot(ax=ax, color='k', label='PC1')
# ax.axhline(0, c='0.8')
# #ax.set_xlabel('period', fontsize=18)
# ax.plot(PCdf.index, detrend(PCdf.ix[:,0].values), 'r',  label='PC1 (trend removed)')
# ax.grid('off')
# ax.legend(loc=1);







#
#
# subset1 = files[0:-1].copy()
#
# colormap = plt.cm.gist_ncar
#
# ax[1].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subset1)))))
# labels = []
# for i in range(len(subset1)):
#
#     data = getBathy(os.path.join(geomorphdir, subset1[i]), lower=-5, upper=10)
#     temp = subset1[i].split('_')
#     ax[1].plot(data['x'], data['z'], label=temp[1])
#     #labels.append(temp[1])
#
# ax[1].legend(loc='upper right')
# ax[1].set_xlim([50, 800])
# ax[1].set_ylim([-8, 4])
# ax[1].set_title('Profile variability post-nourishment')


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












