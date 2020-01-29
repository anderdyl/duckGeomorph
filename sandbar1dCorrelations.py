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
#from downloads import pyeemd



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


def interpBathy(xSub, zSub, x):

    f = interp1d(xSub, zSub, kind='linear', bounds_error=False)
    newz = f(x)
    newBathy = dict()
    newBathy['x'] = x
    newBathy['z'] = newz
    return newBathy


subset = files[0:960].copy()
#from palettable.colorbrewer.sequential import Blues_8
#ax.imshow(data, cmap=Blues_8.mpl_colormap

colormap = plt.cm.gist_ncar


labels = []
xinterp = np.arange(100, 650, 2.5)

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


            xindex = np.where((xSub>=500))

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



fig, ax = plt.subplots(1,2)
ax[0].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(alllines)))))
for i in range(len(alllines)):
    ax[0].plot(xinterp, alllines[i,:], label=time[i])

#ax[0].legend(loc='upper right')
ax[0].set_xlim([50, 750])
ax[0].set_ylim([-8, 4])
ax[0].set_title('Cross-shore profile variability')

tg, xg = np.meshgrid(time, xinterp)
ax[1].pcolor(xg,tg,alllines.T)

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
        xFRFbar, xFRFtrough = mL.findSandBarAndTrough1D(bathyX, bathy, plotFname=fname, smoothLengthScale=15, profileTrend=np.mean(alllines,axis=0))
        #xFRFbar, xFRFtrough = mL.findSandBarAndTrough1D(bathyX, bathy, plotFname=fname, smoothLengthScale=10)

        if xFRFbar is not None:
            for sandbarX in xFRFbar:
                ax[1].plot(sandbarX, time[tt], 'ro', label='bar')

            if len(xFRFbar) == 3:
                outerBar[tt] = np.max(xFRFbar)
                zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == bathyX))
                zOuterBar[tt] = bathy[zOuterInd[0]]

                innerBar[tt] = np.median(xFRFbar)
                zInnerInd = np.where((np.round(2*np.median(xFRFbar))/2 == bathyX))
                zInnerBar[tt] = bathy[zInnerInd[0]]

            elif len(xFRFbar) == 2:
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
            for troughX in xFRFtrough:
                ax[1].plot(troughX, time[tt], 'bd', label='trough')
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



fig2,ax2 = plt.subplots(3,1)#(figsize=(8,8))
ax2[0].plot(time,innerBar,'go')
ax2[1].plot(time,outerBar,'mo')
ax2[2].plot(time[1:],np.diff(outerBar))

Ab = zOuterBar-zOuterTrough
Ab2 = zInnerBar-zInnerTrough

fig3 = plt.figure(figsize=(8,8))
plt.plot(time,Ab,'-')
