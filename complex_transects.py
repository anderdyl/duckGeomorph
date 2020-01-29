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

def interpBathy2(xSub,zSub,x):
    from Loess import Loess

    loess = Loess(xSub, zSub)
    zloess = np.nan * np.ones(len(x),)
    count = 0
    for xxx in x:
        zloess[count] = loess.estimate(xxx, window=10, use_matrix=False, degree=1)
        count = count+1
    newBathy = dict()
    newBathy['x'] = x
    newBathy['z'] = zloess
    return newBathy


subset = files[0:970].copy()
#from palettable.colorbrewer.sequential import Blues_8
#ax.imshow(data, cmap=Blues_8.mpl_colormap

colormap = plt.cm.gist_ncar

import scipy.stats as spstats

labels = []
xinterp = np.arange(110, 560, 3)

bathy = dict()
#alllines = np.empty((len(xinterp),))
count = 0
count2 = 0
count3 = 0
count4 = 0
count5 = 0
worstcount = 0
#fig = plt.figure(figsize=(10,10))
for i in range(len(subset)):
    data = getBathy(os.path.join(geomorphdir, subset[i]), lower=1070, upper=1100)
    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=750, upper=950)
    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-10, upper=20)

    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-2, upper=10)

    temp = subset[i].split('_')

    surveydate = DT.datetime.strptime(temp[1], '%Y%m%d')
    elevs = data['z']
    cross = data['x']
    crossind = np.argsort(data['x'])
    crossS = cross[crossind]
    elevsS = elevs[crossind]

    xSub = np.ma.MaskedArray.filled(crossS, np.nan)
    zSub = np.ma.MaskedArray.filled(elevsS, np.nan)

    #realValues = ~np.isnan(xSub)

    #binnedz, bin_edges, binnumber = spstats.binned_statistic(xSub,zSub, statistic='mean',bins=np.arange(106,566,4))
    #bin_width = (bin_edges[1] - bin_edges[0])
    #bin_centers = bin_edges[1:] - bin_width / 2
    #xMSub = bin_centers
    #zMSub = binnedz
    #xSubNew = xMSub[~np.isnan(zMSub)]
    #zSubNew = zMSub[~np.isnan(zMSub)]

    realValues = ~np.isnan(xSub)
    xSubNew = xSub[~np.isnan(zSub)]
    zSubNew = zSub[~np.isnan(zSub)]



    if any(realValues):

        if np.nanmax(xSubNew) > 500:

            if np.nanmin(xSubNew) < 125:

                if len(xSubNew) > 10:
                    #newdata = interpBathy(xSubNew, zSubNew, xinterp)
                    #newdata = interpBathy2(newdata['x'], newdata['z'], xinterp)
                    newdata = interpBathy2(xSubNew, zSubNew, xinterp)
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
                            count2 = count2+1
                            badtime = np.append(badtime, surveydate)

                        #print('Threw out a survey with data due to nans: {}'.format(surveydate))
                        #plt.plot(xSub,zSub)


                        xindex = np.where((xSubNew>=700))

                        if len(xindex[0]) > 3:

                            xbot = xSubNew[xindex]
                            zbot = zSubNew[xindex]
                            mb = np.polyfit(xbot, zbot, 1)
                            f = np.poly1d(mb)
                            maxXind = np.where((xinterp > np.nanmax(xSubNew)))
                            newX = xinterp[maxXind]
                            newZ = f(newX)
                            xSubNew2 = np.append(xSubNew, newX)
                            zSubNew2 = np.append(zSubNew, newZ)
                            moredata = interpBathy(xSubNew2, zSubNew2, xinterp)
                            moredata = interpBathy(moredata['x'], moredata['z'], xinterp)
                            #moredata = interpBathy(xSubNew2, moredata['z'], xinterp)

                            del xSubNew2, zSubNew2,newZ,newX,f,mb,zbot,xbot

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


                                if count3 == 0:
                                    extrapolatedlines = newdata['z']
                                    count3 = count3 + 1
                                    extrapolatedtime = surveydate
                                else:
                                    extrapolatedlines = np.vstack((extrapolatedlines, newdata['z']))
                                    count3 = count3 + 1
                                    extrapolatedtime = np.append(extrapolatedtime, surveydate)
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
                    print('Data is less than 10 points long at line {}'.format(i))
            else:
                print('No data onshore of 125 meters at line {}'.format(i))
                print('Most onshore point at {}'.format(np.nanmin(xSubNew)))
                newdata = interpBathy2(xSubNew, zSubNew, xinterp)

                # if count4 == 0:
                #     noneOnshore = newdata['z']
                #     noOnshoreTime = surveydate
                #     count4 = count4 + 1
                # else:
                #     noneOnshore = np.vstack((noneOnshore, newdata['z']))
                #     noOnshoreTime = np.append(noOnshoreTime, surveydate)
                #     count4 = count4 + 1
        else:
            print('No data deeper than 500 meters at line {}'.format(i))
            print('Most offshore point at {}'.format(np.nanmax(xSubNew)))
            newdata = interpBathy2(xSubNew, zSubNew, xinterp)

            # if count5 == 0:
            #     noneOffshore = newdata['z']
            #     noOffshoreTime = surveydate
            #     count5 = count5 + 1
            # else:
            #     noneOffshore = np.vstack((noneOffshore, newdata['z']))
            #     noOffshoreTime = np.append(noOffshoreTime, surveydate)
            #     count4 = count5 + 1
    else:
        print('Survey with no data at this line {}'.format(i))

        #plt.plot(xSubNew,zSubNew)


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
#
#
# #fig20 = plt.figure(figsize=(10,5))
# #plt.plot(time,alllines[:,-1])
# #plt.plot(time,alllines[:,-10])
# #plt.plot(time,alllines[:,-20])
# #plt.plot(time,alllines[:,-30])
# #plt.plot(time,alllines[:,-40])

#file = '/home/dylananderson/projects/alllinesSouth.nc'
#bathy = Dataset(file)
#alllines = bathy.variables['alllines'][:]
#alllines = np.ma.MaskedArray.filled(alllines)

#time = bathy.variables['t'][:]
#     zs_bathy = bathy.variables['elevation'][:]

demean = alllines - np.mean(alllines,axis=0)

from scipy.signal import hilbert
#from scipy.fftpack import hilbert
data = (hilbert(demean.T))

#datareal = np.imag(hilbert(alllines))
#dataimag = np.real(hilbert(alllines))
data = data.T
c = np.matmul(np.conj(data).T,data)/np.shape(data)[0]


import scipy.linalg as la
import numpy.linalg as npla

lamda, loadings = la.eigh(c)

lamda2, loadings2 = npla.eig(c)

ind = np.argsort(lamda[::-1])

lamda[::-1].sort()

loadings = loadings[:,ind]

pcs = np.dot(data, loadings)# / np.sqrt(lamda)
loadings = loadings# * np.sqrt(lamda)
pcsreal = np.real(pcs[:,0:200])
pcsimag = np.imag(pcs[:,0:200])
eofreal = np.real(loadings[:,0:200])
eofimag = np.imag(loadings[:,0:200])
S = np.power(loadings*np.conj(loadings),0.5) * np.sqrt(lamda)

theta = np.arctan2(eofimag,eofreal)
theta2 = theta*180/np.pi

Rt = np.power(pcs*np.conj(pcs),0.5) / np.sqrt(lamda)

phit = np.arctan2(pcsimag,pcsreal)
phit2 = phit*180/np.pi

mode = 0

fig, ax = plt.subplots(2,2)

ax[0,0].plot(xinterp, S[:,mode],'o')
ax[0,1].plot(xinterp, theta2[:,mode],'o')
ax[1,0].plot(time,Rt[:,mode],'o')
ax[1,1].plot(time,phit2[:,mode],'o')

PC1 = Rt[:, mode]*np.sin(phit[:, mode]) + Rt[:, mode]*np.cos(phit[:, mode])
PC1a = Rt[:, mode]*np.sin(phit[:, mode])
PC1b = Rt[:, mode]*np.cos(phit[:, mode])
#
# plt.figure()
# plt.plot(time,PC1,'.')
# plt.plot(time,PC1a,'.')
# plt.plot(time,PC1b,'.')


totalV = np.sum(lamda)
percentV = lamda / totalV

ztemp = 0*np.ones(len(xinterp),)
timestep = 200
for mode in range(5):
    ztemp = ztemp + Rt[timestep,mode]*np.sin(phit[timestep,mode]) * S[:,mode]*np.sin(theta[:,mode]) + Rt[timestep,mode]*np.cos(phit[timestep,mode]) * S[:,mode]*np.cos(theta[:,mode])
fig2 = plt.figure()
#plt.plot(xinterp,np.mean(alllines,axis=0)+ztemp+ztemp1+ztemp2+ztemp3+ztemp4+ztemp5+ztemp6+ztemp7+ztemp8)
plt.plot(xinterp,np.mean(alllines,axis=0)+ztemp)
#plt.plot(xinterp,np.mean(alllines,axis=0)+ztemp2)
plt.plot(xinterp,np.mean(alllines,axis=0))
plt.plot(xinterp,alllines[timestep,:])

#plt.plot(xinterp,ztemp)
#plt.plot(xinterp,ztemp1)
#plt.plot(xinterp,ztemp2)
#plt.plot(xinterp,ztemp+ztemp1+ztemp2)


def P2R(radii, angles):
    return radii * np.exp(1j*angles)

def R2P(x):
    return np.abs(x), np.angle(x)



timeind = np.arange(0,len(time))
#timeind = np.arange(3, 50)
RtSubset = Rt[timeind, :]
phitSubset = phit[timeind, :]
phit2Subset = phit2[timeind, :]
timeSubset = time[timeind]
alllinesSubset = alllines[timeind, :]

eofPred = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred2 = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred3 = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred4 = np.nan * np.ones((np.shape(alllinesSubset)))

for timestep in range(len(timeind)):
    mode = 0
    eofPred[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
            timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    mode = 1
    eofPred2[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
            timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    mode = 2
    eofPred3[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
            timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    mode = 3
    eofPred4[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
            timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    # ztemp = 0 * np.ones(len(xinterp), )
    # for mode in range(8):
    #     ztemp = ztemp + RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
    #         timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])


t1 = 500
t2 = -1

fig, ax = plt.subplots(1,5)
#plt.set_cmap('RdBu')#bwr')
plt.set_cmap('bwr')

tg, xg = np.meshgrid(time, xinterp)
plt0 = ax[0].pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
fig.colorbar(plt0, ax=ax[0], orientation='horizontal')
ax[0].set_ylim([time[t1], time[t2]])
ax[0].set_title('Surveys (dev.)')

plt1 = ax[1].pcolor(xg,tg,eofPred.T, vmin=-1.05, vmax=1.05)
ax[1].set_ylim([time[t1], time[t2]])
fig.colorbar(plt1, ax=ax[1], orientation='horizontal')
ax[1].set_title('CEOF1 {:.2f}'.format(percentV[0]))
ax[1].get_yaxis().set_ticks([])

plt2 = ax[2].pcolor(xg,tg,eofPred2.T, vmin=-.65, vmax=.65)
ax[2].set_ylim([time[t1], time[t2]])
fig.colorbar(plt2, ax=ax[2], orientation='horizontal')
ax[2].set_title('CEOF2 {:.2f}'.format(percentV[1]))
ax[2].get_yaxis().set_ticks([])

plt3 = ax[3].pcolor(xg,tg,eofPred3.T, vmin=-.65, vmax=.65)
ax[3].set_ylim([time[t1], time[t2]])
ax[3].set_title('CEOF3 {:.2f}'.format(percentV[2]))
ax[3].get_yaxis().set_ticks([])

fig.colorbar(plt3, ax=ax[3], orientation='horizontal')

plt4 = ax[4].pcolor(xg,tg,eofPred4.T, vmin=-.65, vmax=.65)
ax[4].set_ylim([time[t1], time[t2]])
ax[4].set_title('CEOF4 {:.2f}'.format(percentV[3]))
ax[4].get_yaxis().set_ticks([])
fig.colorbar(plt4, ax=ax[4], orientation='horizontal')

#plt.tight_layout(pad=0.5)



wavedir = '/media/dylananderson/Elements/WIS_ST63218/'

# Need to sort the files to ensure correct temporal order...
files = os.listdir(wavedir)
files.sort()
files_path = [os.path.join(os.path.abspath(wavedir), x) for x in files]

wis = Dataset(files_path[0])

def getWIS(file):
    waves = Dataset(file)

    waveHs = waves.variables['waveHs'][:]
    waveTp = waves.variables['waveTp'][:]
    waveTm = waves.variables['waveTm'][:]
    waveTm1 = waves.variables['waveTm1'][:]
    waveTm2 = waves.variables['waveTm2'][:]
    waveMeanDirection = waves.variables['waveMeanDirection'][:]
    waveMeanDirectionSwell = waves.variables['waveMeanDirectionSwell'][:]
    timeW = waves.variables['time'][:]


    output = dict()
    output['waveHs'] = waveHs
    output['waveTp'] = waveTp
    output['waveTm'] = waveTm
    output['waveTm1'] = waveTm1
    output['waveTm2'] = waveTm2
    output['waveMeanDirection'] = waveMeanDirection
    output['waveMeanDirectionSwell'] = waveMeanDirectionSwell
    output['t'] = timeW

    return output


Hs = []
Tm = []
Dm = []
timeWave = []
for i in files_path:
    waves = getWIS(i)
    Hs = np.append(Hs,waves['waveHs'])
    Tm = np.append(Tm,waves['waveTm'])
    Dm = np.append(Dm,waves['waveMeanDirection'])
    timeWave = np.append(timeWave,waves['t'])





# waveNorm = Dm - 72
# neg = np.where((waveNorm > 180))
# waveNorm[neg[0]] = waveNorm[neg[0]]-360
# offpos = np.where((waveNorm>90))
# offneg = np.where((waveNorm<-90))
# waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
# waveNorm[offneg[0]] = waveNorm[offneg[0]]*0
#
# LWP_means = 1025*np.square(Hs)*Tm*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
tWave = [DT.datetime.fromtimestamp(x) for x in timeWave]

#fig = plt.figure(figsize=(10,10))
#plt.plot(tWave,Hs)
#plt.plot(time,np.ones(len(time),),'.')

HsSurvey = np.nan * np.ones(len(time))
HsMaxSurvey = np.nan * np.ones(len(time))
TmSurvey = np.nan * np.ones(len(time))
DmSurvey = np.nan * np.ones(len(time))
WPSurvey = np.nan * np.ones(len(time))
EFSurvey = np.nan * np.ones(len(time))
EFLSurvey = np.nan * np.ones(len(time))

timeDiff = np.nan * np.ones(len(time))
counter = 0
for tt in range(len(time)-1):
    #within = [date for date in tWave if time[tt] < date < time[tt+1]]
    withinIndex = np.where((tWave > np.datetime64(time[tt])) & (tWave < np.datetime64(time[tt+1])))
    if len(withinIndex[0]) < 2:
        print('index is not working after time {}'.format(time[tt]))
        print('which is index = {}'.format(tt))
    else:

        diff = time[tt+1]-time[tt]
        timeDiff[counter] = diff.days
        HsSurvey[counter] = np.nanmean(Hs[withinIndex[0]])
        HsMaxSurvey[counter] = np.nanmax(Hs[withinIndex[0]])
        TmSurvey[counter] = np.nanmean(Tm[withinIndex[0]])
        dirs = Dm[withinIndex[0]]-72
        neg = np.where((dirs > 180))
        dirs[neg[0]] = dirs[neg[0]]-360
        waveNorm = dirs
        offpos = np.where((waveNorm>90))
        offneg = np.where((waveNorm<-90))
        waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
        waveNorm[offneg[0]] = waveNorm[offneg[0]]*0
        DmSurvey[counter] = np.nanmean(dirs)

        c0 = 9.81*Tm[withinIndex[0]]/(2 * np.pi)
        HsSquare = np.square(Hs[withinIndex[0]])
        EF = 1025*9.81*HsSquare*c0/16
        EFlongshore = (1025*9.81*HsSquare*c0/16)*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
        EFSurvey[counter] = np.sum(EF)
        EFLSurvey[counter] = np.nanmean(np.abs(EFlongshore))
        WPSurvey[counter] = np.nanmean(HsSquare*Tm[withinIndex[0]])


        # fig10 = plt.figure(figsize=(12, 12))
        # ax100 = plt.subplot2grid((4, 3), (0, 0), colspan=3, rowspan=1)
        # ax100.plot(tWave[withinIndex[0][0]:withinIndex[0][-1]], Hs[withinIndex[0][0:-1]])
        # ax100.set_ylabel('Hs (m)')
        # ax100.set_title('Duration = {} days'.format(diff.days))
        # ax100.text(tWave[withinIndex[0][0]], HsSurvey[counter], 'Average Hs = {:.2f}'.format(HsSurvey[counter]))
        # ax100.plot([tWave[withinIndex[0][0]], tWave[withinIndex[0][-1]]],[HsSurvey[counter],HsSurvey[counter]])
        # #     ax7.text(320, 0, '{}'.format(timeSubset[timestep]))
        # ax101 = plt.subplot2grid((4, 3), (1, 0), colspan=3, rowspan=1)
        # ax101.plot(tWave[withinIndex[0][0]:withinIndex[0][-1]], Tm[withinIndex[0][0:-1]])
        # ax101.text(tWave[withinIndex[0][0]], TmSurvey[counter], 'Average Tm = {:.2f}'.format(TmSurvey[counter]))
        # ax101.plot([tWave[withinIndex[0][0]], tWave[withinIndex[0][-1]]],[TmSurvey[counter],TmSurvey[counter]])
        # ax101.set_ylabel('Tm (s)')
        # ax103 = plt.subplot2grid((4, 3), (2, 0), colspan=3, rowspan=1)
        # ax103.plot(tWave[withinIndex[0][0]:withinIndex[0][-1]], dirs[0:-1])
        # ax103.text(tWave[withinIndex[0][0]], DmSurvey[counter], 'Average Dm = {:.2f}'.format(DmSurvey[counter]))
        # ax103.plot([tWave[withinIndex[0][0]], tWave[withinIndex[0][-1]]],[DmSurvey[counter],DmSurvey[counter]])
        # ax103.set_ylabel('Dm (deg N)')
        # ax102 = plt.subplot2grid((4, 3), (3, 0), colspan=3, rowspan=1)
        # ax102.plot(xinterp,alllines[tt])
        # ax102.plot(xinterp,alllines[tt+1])
        #
        #
        #
        #
        # if tt < 10:
        #     plt.savefig('Conditions00'+str(tt)+'.png')
        #
        # elif tt < 100:
        #     plt.savefig('Conditions0'+str(tt)+'.png')
        #
        # else:
        #     plt.savefig('Conditions'+str(tt)+'.png')
        # plt.close()
        counter = counter+1



#HsReg = HsSurvey.copy()
HsReg = HsSurvey[~np.isnan(HsSurvey)]
HsRegMax = HsMaxSurvey[~np.isnan(HsSurvey)]
TmReg = TmSurvey[~np.isnan(HsSurvey)]
DmReg = DmSurvey[~np.isnan(HsSurvey)]
WPReg =WPSurvey[~np.isnan(HsSurvey)]
EFReg =EFSurvey[~np.isnan(HsSurvey)]
EFLReg =EFLSurvey[~np.isnan(HsSurvey)]

timeDiff = timeDiff[~np.isnan(HsSurvey)]

wavetime = time[1:len(HsReg)+1]






cpc1 = P2R(Rt[1:len(HsReg)+1, 0], phitSubset[1:len(HsReg)+1, 0])
cpc2 = P2R(Rt[1:len(HsReg)+1, 1], phitSubset[1:len(HsReg)+1, 1])
cpc3 = P2R(Rt[1:len(HsReg)+1, 2], phitSubset[1:len(HsReg)+1, 2])

precpc1 = P2R(Rt[0:len(HsReg), 0], phitSubset[0:len(HsReg), 0])
precpc2 = P2R(Rt[0:len(HsReg), 1], phitSubset[0:len(HsReg), 1])
precpc3 = P2R(Rt[0:len(HsReg), 2], phitSubset[0:len(HsReg), 2])

dcpc1 = P2R(Rt[0:len(HsReg)+1, 0], phitSubset[0:len(HsReg)+1, 0])
dcpc2 = P2R(Rt[0:len(HsReg)+1, 1], phitSubset[0:len(HsReg)+1, 1])
dcpc3 = P2R(Rt[0:len(HsReg)+1, 2], phitSubset[0:len(HsReg)+1, 2])

realD1Dt = np.real(dcpc1[1:])-np.real(dcpc1[0:-1])
imagD1Dt = np.imag(dcpc1[1:])-np.imag(dcpc1[0:-1])
realD2Dt = np.real(dcpc2[1:])-np.real(dcpc2[0:-1])
imagD2Dt = np.imag(dcpc2[1:])-np.imag(dcpc2[0:-1])
realD3Dt = np.real(dcpc3[1:])-np.real(dcpc3[0:-1])
imagD3Dt = np.imag(dcpc3[1:])-np.imag(dcpc3[0:-1])
# fig4, ax4 = plt.subplots(2,1)
# ax4[0].plot(time[1:-1],HsSurvey[0:-2])
# ax4[1].plot(time, np.real(cpc1),'-')
# ax4[1].plot(time, np.imag(cpc1),'-')
#
# ax4[0].set_xlim([time[t1], time[t2]])
# ax4[1].set_xlim([time[t1], time[t2]])

PC1 = Rt[0:len(HsReg)+1, 0]*np.sin(phit[0:len(HsReg)+1, 0]) + Rt[0:len(HsReg)+1, 0]*np.cos(phit[0:len(HsReg)+1, 0])
PC2 = Rt[0:len(HsReg)+1, 1]*np.sin(phit[0:len(HsReg)+1, 1]) + Rt[0:len(HsReg)+1, 1]*np.cos(phit[0:len(HsReg)+1, 1])
PC3 = Rt[0:len(HsReg)+1, 2]*np.sin(phit[0:len(HsReg)+1, 2]) + Rt[0:len(HsReg)+1, 2]*np.cos(phit[0:len(HsReg)+1, 2])

dPC1 = PC1[1:] - PC1[0:-1]
dPC2 = PC2[1:] - PC2[0:-1]
dPC3 = PC3[1:] - PC3[0:-1]

prePC1 = Rt[0:len(HsReg), 0]*np.sin(phit[0:len(HsReg), 0]) + Rt[0:len(HsReg), 0]*np.cos(phit[0:len(HsReg), 0])
prePC2 = Rt[0:len(HsReg), 1]*np.sin(phit[0:len(HsReg), 1]) + Rt[0:len(HsReg), 1]*np.cos(phit[0:len(HsReg), 1])
prePC3 = Rt[0:len(HsReg), 2]*np.sin(phit[0:len(HsReg), 2]) + Rt[0:len(HsReg), 2]*np.cos(phit[0:len(HsReg), 2])

t1 = 300
t2 = -1

fig2, ax2 = plt.subplots(1,4)
ax2[0].plot(EFReg/timeDiff,wavetime)
ax2[1].pcolor(xg,tg,eofPred2.T, vmin=-1.25, vmax=1.25)
ax2[2].plot(Rt[1:len(HsReg)+1,1],wavetime)
ax2[3].plot(phit[1:len(HsReg)+1,1],wavetime)

ax2[0].set_ylim([wavetime[t1], wavetime[t2]])
ax2[1].set_ylim([wavetime[t1], wavetime[t2]])
ax2[2].set_ylim([wavetime[t1], wavetime[t2]])
ax2[3].set_ylim([wavetime[t1], wavetime[t2]])



stackRI = np.vstack((realD1Dt, imagD1Dt, np.real(cpc1), np.imag(cpc1))).T

stackRI2 = np.vstack((realD2Dt, imagD2Dt)).T

from sklearn.linear_model import LinearRegression

model1 = LinearRegression()

#model1.fit(HsReg.reshape((-1, 1)),stackRI)
model1.fit(stackRI, WPReg.reshape((-1, 1)))


r_sq = model1.score(stackRI, WPReg.reshape((-1, 1)))

X = np.vstack((realD1Dt, imagD1Dt, np.real(cpc1), np.imag(cpc1))).T
Y = WPReg

#X2 = np.vstack((np.real(precpc2), np.imag(precpc2), HsReg, TmReg, WPReg)).T
X2 = np.vstack((EFReg/timeDiff, EFLReg/timeDiff, np.real(prePC3))).T
Y2 = np.abs(np.real(dPC3))

from statsmodels import api

X2 = api.add_constant(X2)

model2 = api.OLS(Y2, X2).fit() ## sm.OLS(output, input)
predictions = model2.predict(X2)
model2.summary()

#
# # For making a bunch of figures with 6 of the CEOFs in magnitude/space dimensions
# timeind = np.arange(0,len(time))
# #timeind = np.arange(3, 50)
# RtSubset = Rt[timeind, :]
# phitSubset = phit[timeind, :]
# phit2Subset = phit2[timeind, :]
# timeSubset = time[timeind]
# alllinesSubset = alllines[timeind, :]
# for timestep in range(len(timeind)):
#
#     fig10 = plt.figure(figsize=(12, 12))
#
#     ax1 = plt.subplot2grid((8, 8), (0, 0), colspan=4, rowspan=2)
#     ax1.plot(time, Rt[:, 0], '.')
#     ax1.plot(timeSubset[timestep], RtSubset[timestep, 0], 'ro')
#     ax1.set_ylabel('Amplitude')
#
#     ax2 = plt.subplot2grid((8, 8), (2, 0), colspan=4, rowspan=2)
#     ax2.plot(time, phit2[:, 0], '.')
#     ax2.plot(timeSubset[timestep], phit2Subset[timestep, 0], 'ro')
#     ax2.set_ylabel('Phase')
#
#     ax3 = plt.subplot2grid((8, 8), (0, 4), colspan=4, rowspan=2)
#     eof1 = RtSubset[timestep, 0] * np.sin(phitSubset[timestep, 0]) * S[:, mode] * np.sin(theta[:, 0]) + RtSubset[
#         timestep, 0] * np.cos(phitSubset[timestep, 0]) * S[:, 0] * np.cos(theta[:, 0])
#     ax3.plot(xinterp,eof1)
#     ax3.set_ylim([-1, 1])
#     ax3.set_xlim([100, 675])
#     ax3.set_title('CEOF1')
#
#     ax4 = plt.subplot2grid((8, 8), (2, 4), colspan=4, rowspan=2)
#     cpc1 = P2R(RtSubset[0:timestep+1, 0], phitSubset[0:timestep+1, 0])
#     ax4.plot(np.real(cpc1),np.imag(cpc1),'-')
#     ax4.plot(np.real(cpc1[-1]),np.imag(cpc1[-1]),'o')
#     ax4.set_xlim([-2, 2])
#     ax4.set_ylim([-2, 2])
#
#     ztemp = 0 * np.ones(len(xinterp), )
#     for mode in range(10):
#         ztemp = ztemp + RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
#             timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
#
#     ax7 = plt.subplot2grid((8, 8), (4, 0), colspan=8, rowspan=4)
#     ax7.plot(xinterp, np.mean(alllines, axis=0) + ztemp, label='prediction')
#     ax7.plot(xinterp, alllinesSubset[timestep, :], label='surveyed')
#     ax7.legend()
#     ax7.set_ylim([-7, 1])
#     ax7.set_xlim([100, 675])
#     ax7.text(320, 0.5, 'Summation of CEOFs 1-10')
#     ax7.text(320, 0, '{}'.format(timeSubset[timestep]))
#
#     plt.tight_layout(pad=2.0)
#
#     if timestep < 10:
#         plt.savefig('frame00'+str(timestep)+'.png')
#
#     elif timestep < 100:
#         plt.savefig('frame0'+str(timestep)+'.png')
#
#     else:
#         plt.savefig('frame'+str(timestep)+'.png')
#     plt.close()



#
# # For making a bunch of figures with 6 of the CEOFs in magnitude/space dimensions
# for timestep in range(1):
#
#     fig10 = plt.figure(figsize=(12, 12))
#     ax1 = plt.subplot2grid((4, 3), (0, 0), colspan=1, rowspan=1)
#     eof1 = Rt[timestep, 0] * np.sin(phit[timestep, 0]) * S[:, mode] * np.sin(theta[:, 0]) + Rt[
#         timestep, 0] * np.cos(phit[timestep, 0]) * S[:, 0] * np.cos(theta[:, 0])
#     cpc1 = P2R(Rt[0:timestep,0],phit[0:timestep,0])
#     #ax1.plot(xinterp,eof1)
#     ax1.plot(np.real(cpc1),np.imag(cpc1),'-')
#     ax1.plot(np.real(cpc1[-1]),np.image(cpc1[-1]),'o')
#     ax1.set_ylim([-1, 1])
#     ax1.set_xlim([100, 675])
#     ax1.set_title('CEOF1')
#
#     ax2 = plt.subplot2grid((4, 3), (0, 1), colspan=1, rowspan=1)
#     eof2 = Rt[timestep, 1] * np.sin(phit[timestep, 1]) * S[:, mode] * np.sin(theta[:, 1]) + Rt[
#         timestep, 1] * np.cos(phit[timestep, 1]) * S[:, 1] * np.cos(theta[:, 1])
#     ax2.plot(xinterp,eof2)
#     ax2.set_ylim([-1, 1])
#     ax2.set_xlim([100, 675])
#     ax2.set_title('CEOF2')
#
#     ax3 = plt.subplot2grid((4, 3), (0, 2), colspan=1, rowspan=1)
#     eof3 = Rt[timestep, 2] * np.sin(phit[timestep, 2]) * S[:, mode] * np.sin(theta[:, 2]) + Rt[
#         timestep, 2] * np.cos(phit[timestep, 2]) * S[:, 2] * np.cos(theta[:, 2])
#     ax3.plot(xinterp,eof3)
#     ax3.set_ylim([-1, 1])
#     ax3.set_xlim([100, 675])
#     ax3.set_title('CEOF3')
#
#     ax4 = plt.subplot2grid((4, 3), (1, 0), colspan=1, rowspan=1)
#     eof4 = Rt[timestep, 3] * np.sin(phit[timestep, 3]) * S[:, mode] * np.sin(theta[:, 3]) + Rt[
#         timestep, 3] * np.cos(phit[timestep, 3]) * S[:, 3] * np.cos(theta[:, 3])
#     ax4.plot(xinterp,eof4)
#     ax4.set_ylim([-1, 1])
#     ax4.set_xlim([100, 675])
#     ax4.set_title('CEOF4')
#
#     ax5 = plt.subplot2grid((4, 3), (1, 1), colspan=1, rowspan=1)
#     eof5 = Rt[timestep, 4] * np.sin(phit[timestep, 4]) * S[:, mode] * np.sin(theta[:, 4]) + Rt[
#         timestep, 4] * np.cos(phit[timestep, 4]) * S[:, 4] * np.cos(theta[:, 4])
#     ax5.plot(xinterp,eof5)
#     ax5.set_ylim([-1, 1])
#     ax5.set_xlim([100, 675])
#     ax5.set_title('CEOF5')
#
#     ax6 = plt.subplot2grid((4, 3), (1, 2), colspan=1, rowspan=1)
#     eof6 = Rt[timestep, 5] * np.sin(phit[timestep, 5]) * S[:, mode] * np.sin(theta[:, 5]) + Rt[
#         timestep, 5] * np.cos(phit[timestep, 5]) * S[:, 5] * np.cos(theta[:, 5])
#     ax6.plot(xinterp,eof6)
#     ax6.set_ylim([-1, 1])
#     ax6.set_xlim([100, 675])
#     ax6.set_title('CEOF6')
#
#
#
#     ztemp = 0 * np.ones(len(xinterp), )
#     for mode in range(10):
#         ztemp = ztemp + Rt[timestep, mode] * np.sin(phit[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + Rt[
#             timestep, mode] * np.cos(phit[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
#
#     ax7 = plt.subplot2grid((4, 3), (2, 0), colspan=3, rowspan=2)
#     ax7.plot(xinterp, np.mean(alllines, axis=0) + ztemp)
#     ax7.set_ylim([-7, 1])
#     ax7.set_xlim([100, 675])
#     ax7.text(320, 0.5, 'Summation of CEOFs 1-10')
#     ax7.text(320, 0, '{}'.format(time[timestep]))
#     if timestep < 10:
#         plt.savefig('frame00'+str(timestep)+'.png')
#
#     elif timestep < 100:
#         plt.savefig('frame0'+str(timestep)+'.png')
#
#     else:
#         plt.savefig('frame'+str(timestep)+'.png')
#     plt.close()

#
# # For making a bunch of figures with 6 of the CEOFs in magnitude/space dimensions
# for timestep in range(1):
#
#     fig10 = plt.figure(figsize=(12, 12))
#     ax1 = plt.subplot2grid((4, 3), (0, 0), colspan=1, rowspan=1)
#     eof1 = Rt[timestep, 0] * np.sin(phit[timestep, 0]) * S[:, mode] * np.sin(theta[:, 0]) + Rt[
#         timestep, 0] * np.cos(phit[timestep, 0]) * S[:, 0] * np.cos(theta[:, 0])
#     ax1.plot(xinterp,eof1)
#     ax1.set_ylim([-1, 1])
#     ax1.set_xlim([100, 675])
#     ax1.set_title('CEOF1')
#
#     ax2 = plt.subplot2grid((4, 3), (0, 1), colspan=1, rowspan=1)
#     eof2 = Rt[timestep, 1] * np.sin(phit[timestep, 1]) * S[:, mode] * np.sin(theta[:, 1]) + Rt[
#         timestep, 1] * np.cos(phit[timestep, 1]) * S[:, 1] * np.cos(theta[:, 1])
#     ax2.plot(xinterp,eof2)
#     ax2.set_ylim([-1, 1])
#     ax2.set_xlim([100, 675])
#     ax2.set_title('CEOF2')
#
#     ax3 = plt.subplot2grid((4, 3), (0, 2), colspan=1, rowspan=1)
#     eof3 = Rt[timestep, 2] * np.sin(phit[timestep, 2]) * S[:, mode] * np.sin(theta[:, 2]) + Rt[
#         timestep, 2] * np.cos(phit[timestep, 2]) * S[:, 2] * np.cos(theta[:, 2])
#     ax3.plot(xinterp,eof3)
#     ax3.set_ylim([-1, 1])
#     ax3.set_xlim([100, 675])
#     ax3.set_title('CEOF3')
#
#     ax4 = plt.subplot2grid((4, 3), (1, 0), colspan=1, rowspan=1)
#     eof4 = Rt[timestep, 3] * np.sin(phit[timestep, 3]) * S[:, mode] * np.sin(theta[:, 3]) + Rt[
#         timestep, 3] * np.cos(phit[timestep, 3]) * S[:, 3] * np.cos(theta[:, 3])
#     ax4.plot(xinterp,eof4)
#     ax4.set_ylim([-1, 1])
#     ax4.set_xlim([100, 675])
#     ax4.set_title('CEOF4')
#
#     ax5 = plt.subplot2grid((4, 3), (1, 1), colspan=1, rowspan=1)
#     eof5 = Rt[timestep, 4] * np.sin(phit[timestep, 4]) * S[:, mode] * np.sin(theta[:, 4]) + Rt[
#         timestep, 4] * np.cos(phit[timestep, 4]) * S[:, 4] * np.cos(theta[:, 4])
#     ax5.plot(xinterp,eof5)
#     ax5.set_ylim([-1, 1])
#     ax5.set_xlim([100, 675])
#     ax5.set_title('CEOF5')
#
#     ax6 = plt.subplot2grid((4, 3), (1, 2), colspan=1, rowspan=1)
#     eof6 = Rt[timestep, 5] * np.sin(phit[timestep, 5]) * S[:, mode] * np.sin(theta[:, 5]) + Rt[
#         timestep, 5] * np.cos(phit[timestep, 5]) * S[:, 5] * np.cos(theta[:, 5])
#     ax6.plot(xinterp,eof6)
#     ax6.set_ylim([-1, 1])
#     ax6.set_xlim([100, 675])
#     ax6.set_title('CEOF6')
#
#
#
#     ztemp = 0 * np.ones(len(xinterp), )
#     for mode in range(10):
#         ztemp = ztemp + Rt[timestep, mode] * np.sin(phit[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + Rt[
#             timestep, mode] * np.cos(phit[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
#
#     ax7 = plt.subplot2grid((4, 3), (2, 0), colspan=3, rowspan=2)
#     ax7.plot(xinterp, np.mean(alllines, axis=0) + ztemp)
#     ax7.set_ylim([-7, 1])
#     ax7.set_xlim([100, 675])
#     ax7.text(320, 0.5, 'Summation of CEOFs 1-10')
#     ax7.text(320, 0, '{}'.format(time[timestep]))
#     if timestep < 10:
#         plt.savefig('frame00'+str(timestep)+'.png')
#
#     elif timestep < 100:
#         plt.savefig('frame0'+str(timestep)+'.png')
#
#     else:
#         plt.savefig('frame'+str(timestep)+'.png')
#     plt.close()






# # For making a bunch of figures with only 3 of the ceofs
# for timestep in range(500):
#
#     fig10 = plt.figure(figsize=(12, 12))
#     ax1 = plt.subplot2grid((6, 6), (0, 0), colspan=2, rowspan=2)
#     eof1 = Rt[timestep, 0] * np.sin(phit[timestep, 0]) * S[:, mode] * np.sin(theta[:, 0]) + Rt[
#         timestep, 0] * np.cos(phit[timestep, 0]) * S[:, 0] * np.cos(theta[:, 0])
#     ax1.plot(xinterp,eof1)
#     ax1.set_ylim([-1, 1])
#     ax1.set_xlim([100, 675])
#     ax1.set_title('CEOF1')
#
#     ax2 = plt.subplot2grid((6, 6), (0, 2), colspan=2, rowspan=2)
#     eof2 = Rt[timestep, 1] * np.sin(phit[timestep, 1]) * S[:, mode] * np.sin(theta[:, 1]) + Rt[
#         timestep, 1] * np.cos(phit[timestep, 1]) * S[:, 1] * np.cos(theta[:, 1])
#     ax2.plot(xinterp,eof2)
#     ax2.set_ylim([-1, 1])
#     ax2.set_xlim([100, 675])
#     ax2.set_title('CEOF2')
#
#     ax3 = plt.subplot2grid((6, 6), (0, 4), colspan=2, rowspan=2)
#     eof3 = Rt[timestep, 2] * np.sin(phit[timestep, 2]) * S[:, mode] * np.sin(theta[:, 2]) + Rt[
#         timestep, 2] * np.cos(phit[timestep, 2]) * S[:, 2] * np.cos(theta[:, 2])
#     ax3.plot(xinterp,eof3)
#     ax3.set_ylim([-1, 1])
#     ax3.set_xlim([100, 675])
#     ax3.set_title('CEOF3')
#     ztemp = 0 * np.ones(len(xinterp), )
#     for mode in range(10):
#         ztemp = ztemp + Rt[timestep, mode] * np.sin(phit[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + Rt[
#             timestep, mode] * np.cos(phit[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
#
#     ax4 = plt.subplot2grid((6, 6), (2, 0), colspan=6, rowspan=4)
#     ax4.plot(xinterp, np.mean(alllines, axis=0) + ztemp)
#     ax4.set_ylim([-7, 1])
#     ax4.set_xlim([100, 675])
#     ax4.text(320, 0.5, 'Summation of CEOFs 1-10')
#     ax4.text(320, 0, '{}'.format(time[timestep]))
#     if timestep < 10:
#         plt.savefig('frame00'+str(timestep)+'.png')
#
#     elif timestep < 100:
#         plt.savefig('frame0'+str(timestep)+'.png')
#
#     else:
#         plt.savefig('frame'+str(timestep)+'.png')
#     plt.close()








# from sklearn.cluster import KMeans
# num_clusters = 20
# PCsub = np.hstack((np.real(Rt[:,0:10]),phit[:,0:10]))
#
# #  KMEANS
# kma = KMeans(n_clusters=num_clusters, n_init=2000).fit(PCsub)
#
# # groupsize
# _, group_size = np.unique(kma.labels_, return_counts=True)
#
# # groups
# d_groups = {}
# for k in range(num_clusters):
#     d_groups['{0}'.format(k)] = np.where(kma.labels_ == k)

#
# # centroids
# centroids = np.dot(kma.cluster_centers_, EOFsub)
#
# # km, x and var_centers
# km = np.multiply(
#     centroids,
#     np.tile(var_anom_std, (num_clusters, 1))
# ) + np.tile(var_anom_mean, (num_clusters, 1))
#
# # sort kmeans
# kma_order = np.argsort(np.mean(-km, axis=1))
#
# # reorder clusters: bmus, km, cenEOFs, centroids, group_size
# sorted_bmus = np.zeros((len(kma.labels_),), ) * np.nan
# for i in range(num_clusters):
#     posc = np.where(kma.labels_ == kma_order[i])
#     sorted_bmus[posc] = i
# sorted_km = km[kma_order]
# sorted_cenEOFs = kma.cluster_centers_[kma_order]
# sorted_centroids = centroids[kma_order]
# sorted_group_size = group_size[kma_order]




