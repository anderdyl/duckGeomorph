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


subset = files[0:972].copy()
#from palettable.colorbrewer.sequential import Blues_8
#ax.imshow(data, cmap=Blues_8.mpl_colormap

colormap = plt.cm.gist_ncar

import scipy.stats as spstats

labels = []
# xinterp = np.arange(108, 608, 2.5)
#xinterp = np.arange(108, 618, 2.5)
xinterp = np.arange(103, 608, 4)
deeperThan = 593
# bathy = dict()
# #alllines = np.empty((len(xinterp),))
# count = 0
# count2 = 0
# count3 = 0
# count4 = 0
# count5 = 0
# worstcount = 0
# #fig = plt.figure(figsize=(10,10))
# for i in range(len(subset)):
#     #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=1070, upper=1100)
#     #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=750, upper=950)
#     #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-10, upper=20)
#     data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-10, upper=20)
#     #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=990, upper=1100)
#
#     #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-2, upper=10)
#
#     temp = subset[i].split('_')
#
#     surveydate = DT.datetime.strptime(temp[1], '%Y%m%d')
#     elevs = data['z']
#     cross = data['x']
#     crossind = np.argsort(data['x'])
#     crossS = cross[crossind]
#     elevsS = elevs[crossind]
#
#     xSub = np.ma.MaskedArray.filled(crossS, np.nan)
#     zSub = np.ma.MaskedArray.filled(elevsS, np.nan)
#
#     #realValues = ~np.isnan(xSub)
#
#     binnedz, bin_edges, binnumber = spstats.binned_statistic(xSub,zSub, statistic='mean',bins=np.arange(98,616,4))
#     bin_width = (bin_edges[1] - bin_edges[0])
#     bin_centers = bin_edges[1:] - bin_width / 2
#     xMSub = bin_centers
#     zMSub = binnedz
#     xSubNewB = xMSub[~np.isnan(zMSub)]
#     zSubNewB = zMSub[~np.isnan(zMSub)]
#     realValues = ~np.isnan(xSubNewB)
#     zSubNew = zSubNewB[~np.isnan(xSubNewB)]
#     xSubNew = xSubNewB[~np.isnan(xSubNewB)]
#
#     # realValues = ~np.isnan(xSub)
#     # xSubNew = xSub[~np.isnan(xSub)]
#     # zSubNew = zSub[~np.isnan(xSub)]
#
#
#
#     if any(realValues):
#
#         if np.nanmax(xSubNew) > deeperThan:
#
#             if np.nanmin(xSubNew) < 125:
#
#                 if len(xSubNew) > 10:
#                     #newdata = interpBathy(xSubNew, zSubNew, xinterp)
#                     #newdata = interpBathy2(newdata['x'], newdata['z'], xinterp)
#                     newdata = interpBathy2(xSubNew, zSubNew, xinterp)
#                     nanValues = np.isnan(newdata['z'])
#
#                     if any(nanValues):
#                         #print('Found a transect with nans {}'.format(i))
#                         print('Trying to extend the lower half of the profile with a linear fit: {}'.format(surveydate))
#                         if count2 == 0:
#                             badlines = newdata['z']
#                             count2 = count2+1
#                             badtime = surveydate
#                         else:
#                             badlines = np.vstack((badlines, newdata['z']))
#                             count2 = count2+1
#                             badtime = np.append(badtime, surveydate)
#
#                         #print('Threw out a survey with data due to nans: {}'.format(surveydate))
#                         #plt.plot(xSub,zSub)
#
#
#                         xindex = np.where((xSubNew>=deeperThan))
#
#                         if len(xindex[0]) > 3:
#
#                             xbot = xSubNew[xindex]
#                             zbot = zSubNew[xindex]
#                             mb = np.polyfit(xbot, zbot, 1)
#                             f = np.poly1d(mb)
#                             maxXind = np.where((xinterp > np.nanmax(xSubNew)))
#                             newX = xinterp[maxXind]
#                             newZ = f(newX)
#                             xSubNew2 = np.append(xSubNew, newX)
#                             zSubNew2 = np.append(zSubNew, newZ)
#                             moredata = interpBathy(xSubNew2, zSubNew2, xinterp)
#                             #moredata = interpBathy(moredata['x'], moredata['z'], xinterp)
#                             #moredata = interpBathy(xSubNew2, moredata['z'], xinterp)
#
#                             del xSubNew2, zSubNew2,newZ,newX,f,mb,zbot,xbot
#
#                             nanValues2 = np.isnan(moredata['z'])
#
#                             if any(nanValues2):
#                                 print('WHAT HAPPENED AT TRANSECT {}'.format(i))
#                                 #plt.plot(moredata['x'], moredata['z'], 'r-')
#                             else:
#                                 if count == 0:
#                                     alllines = moredata['z']
#                                     time = surveydate
#                                     count = count + 1
#                                 else:
#                                     alllines = np.vstack((alllines, moredata['z']))
#                                     time = np.append(time, surveydate)
#                                     count = count + 1
#
#
#                                 if count3 == 0:
#                                     extrapolatedlines = newdata['z']
#                                     count3 = count3 + 1
#                                     extrapolatedtime = surveydate
#                                 else:
#                                     extrapolatedlines = np.vstack((extrapolatedlines, newdata['z']))
#                                     count3 = count3 + 1
#                                     extrapolatedtime = np.append(extrapolatedtime, surveydate)
#                                 #plt.plot(moredata['x'], moredata['z'], 'k-')
#                                 del moredata
#
#
#
#                         elif len(xindex[0]) == 3:
#                             print('Could not do that: number of points being fit: {}'.format(len(xindex[0])))
#                             worstcount = worstcount+1
#                         else:
#                             print('Could not do that: maximum X value in transect is: {}'.format(np.nanmax(xSub)))
#                             worstcount = worstcount+1
#                             #plt.plot(xSub, zSub, 'r-')
#
#
#                     else:
#                         if count == 0:
#                             alllines = newdata['z']
#                             time = surveydate
#                             count = count+1
#                         else:
#                             alllines = np.vstack((alllines, newdata['z']))
#                             time = np.append(time,surveydate)
#                             count = count+1
#
#                         # fig = plt.figure(figsize=(10,10))
#                         # plt.plot(cross,elevs,'.')
#                         # plt.plot(xinterp,newdata['z'],'.-')
#                         # plt.xlim([0, 600])
#                         # plt.ylim([-7,2])
#                         # plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
#                         # plt.savefig('/home/dylananderson/projects/duckGeomorph/interpolatedProfiles/Survey_{}'.format(i))
#                         # plt.close()
#                 else:
#                     print('Data is less than 10 points long at line {}'.format(i))
#             else:
#                 print('No data onshore of 125 meters at line {}'.format(i))
#                 print('Most onshore point at {}'.format(np.nanmin(xSubNew)))
#                 newdata = interpBathy2(xSubNew, zSubNew, xinterp)
#
#                 # if count4 == 0:
#                 #     noneOnshore = newdata['z']
#                 #     noOnshoreTime = surveydate
#                 #     count4 = count4 + 1
#                 # else:
#                 #     noneOnshore = np.vstack((noneOnshore, newdata['z']))
#                 #     noOnshoreTime = np.append(noOnshoreTime, surveydate)
#                 #     count4 = count4 + 1
#         else:
#             print('No data deeper than 500 meters at line {}'.format(i))
#             print('Most offshore point at {}'.format(np.nanmax(xSubNew)))
#             #newdata = interpBathy2(xSubNew, zSubNew, xinterp)
#
#             # if count5 == 0:
#             #     noneOffshore = newdata['z']
#             #     noOffshoreTime = surveydate
#             #     count5 = count5 + 1
#             # else:
#             #     noneOffshore = np.vstack((noneOffshore, newdata['z']))
#             #     noOffshoreTime = np.append(noOffshoreTime, surveydate)
#             #     count4 = count5 + 1
#     else:
#         print('Survey with no data at this line {}'.format(i))

        #plt.plot(xSubNew,zSubNew)
plt.style.use('dark_background')


bathy = dict()
#alllines = np.empty((len(xinterp),))
count = 0
count2 = 0
worstcount = 0
#fig = plt.figure(figsize=(10,10))
for i in range(len(subset)):
    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=1080, upper=1100)
    data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-10, upper=20)

    temp = subset[i].split('_')

    surveydate = DT.datetime.strptime(temp[1], '%Y%m%d')

    # fig = plt.figure(figsize=(10,10))
    # plt.plot(data['x'],data['y'],'o')
    # plt.plot([500,500],[-100,1200],'w--')
    # plt.xlim([0, 800])
    # plt.ylim([-100,1200])
    # plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
    # plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyExtents/Survey_{}'.format(i))
    # plt.close()

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

# Plotting all of the lines
# fig, ax = plt.subplots(2,1)
# ax[0].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(alllines)))))
# for i in range(len(alllines)):
#     ax[0].plot(xinterp, alllines[i,:], label=time[i])
#
# #ax[0].legend(loc='upper right')
# ax[0].set_xlim([50, 850])
# ax[0].set_ylim([-8, 4])
# ax[0].set_title('Cross-shore profile variability')

tg, xg = np.meshgrid(time, xinterp)
# ax[1].pcolor(xg,tg,alllines.T)
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
data = (hilbert(demean.T))

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



# plt.style.use('default')
# mode = 0
# fig = plt.figure(figsize=(11,6))
# ax1 = plt.subplot2grid((2,2),(0,0),rowspan=1,colspan=1)
# ax1.plot(xinterp, S[:,mode],'.')
# ax1.set_ylabel('amplitude')
# ax1.set_xlabel('cross-shore (m)')
# text = ax1.text(-0.15,1.05, "A)", transform=ax1.transAxes)
#
# ax2 = plt.subplot2grid((2,2),(1,0),rowspan=1,colspan=1)
# ax2.plot(xinterp, theta2[:,mode],'.')
# ax2.set_ylabel('phase')
# ax2.set_xlabel('cross-shore (m)')
# text2 = ax2.text(-0.15,1.05, "B)", transform=ax2.transAxes)
#
# ax3 = plt.subplot2grid((2,2),(0,1),rowspan=1,colspan=1)
# ax3.plot(time,Rt[:,mode],'.')
# ax3.set_ylabel('magnitude')
# ax3.set_xlabel('time')
# text3 = ax3.text(-0.15,1.05, "C)", transform=ax3.transAxes)
#
# ax4 = plt.subplot2grid((2,2),(1,1),rowspan=1,colspan=1)
# ax4.plot(time,phit2[:,mode],'.')
# ax4.set_ylabel('phase')
# ax4.set_xlabel('time')
# text4 = ax4.text(-0.15,1.05, "D)", transform=ax4.transAxes)
#

#
# ax1.set_xlim([120,600])
# ax1.set_title('Cross-shore Surveys')
# ax1.set_ylabel('Years')
# ax1.set_xlabel('cross-shore distance (m)')
# fig.subplots_adjust(right=0.84)
# cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
# cbar = fig.colorbar(plt0, cax=cbar_ax)
# # cbar.set_label('WE (Hs$^2$Tp)')
# cbar.set_label('Deviation from mean profile (m)')




mode = 1

# fig, ax = plt.subplots(2,2)
#
# ax[0,0].plot(xinterp, S[:,mode],'o')
# ax[0,1].plot(xinterp, theta2[:,mode],'o')
# ax[1,0].plot(time,Rt[:,mode],'o')
# ax[1,1].plot(time,phit2[:,mode],'o')



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
for mode in range(8):
    ztemp = ztemp + Rt[timestep,mode]*np.sin(phit[timestep,mode]) * S[:,mode]*np.sin(theta[:,mode]) + Rt[timestep,mode]*np.cos(phit[timestep,mode]) * S[:,mode]*np.cos(theta[:,mode])

# Plotting a comparison
# fig2 = plt.figure()
# #plt.plot(xinterp,np.mean(alllines,axis=0)+ztemp+ztemp1+ztemp2+ztemp3+ztemp4+ztemp5+ztemp6+ztemp7+ztemp8)
# plt.plot(xinterp,np.mean(alllines,axis=0)+ztemp)
# #plt.plot(xinterp,np.mean(alllines,axis=0)+ztemp2)
# plt.plot(xinterp,np.mean(alllines,axis=0))
# plt.plot(xinterp,alllines[timestep,:])

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



t1 = 0
t2 = 300
#
# plt.style.use('default')
#
# fig = plt.figure(figsize=(7,11))
# ax1 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
# plt.set_cmap('RdBu_r')
#
# tg, xg = np.meshgrid(time, xinterp)
# plt0 = ax1.pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1, vmax=1)
# # cb = fig.colorbar(plt0, ax=ax1, fraction=0.045, pad=0.1)
# ax1.set_ylim([time[t1], time[t2]])
# ax1.set_xlim([120,600])
# ax1.set_title('Cross-shore Surveys')
# ax1.set_ylabel('Years')
# ax1.set_xlabel('cross-shore distance (m)')
# fig.subplots_adjust(right=0.84)
# cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
# cbar = fig.colorbar(plt0, cax=cbar_ax)
# # cbar.set_label('WE (Hs$^2$Tp)')
# cbar.set_label('Deviation from mean profile (m)')



t1 =495
t2 = -24
#
# plt.style.use('default')
#
# fig = plt.figure()
# ax1 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
# #plt.set_cmap('RdBu')#bwr')
# plt.set_cmap('bwr')
#
# tg, xg = np.meshgrid(time, xinterp)
# plt0 = ax1.pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.6, vmax=1.6)
# fig.colorbar(plt0, ax=ax1, orientation='horizontal')
# ax1.set_ylim([time[t1], time[t2]])
# ax1.set_title('Cross-shore Surveys (deviation from mean profile)')
#
#
# fig, ax = plt.subplots(1,5)
# #plt.set_cmap('RdBu')#bwr')
# plt.set_cmap('RdBu_r')
#
# tg, xg = np.meshgrid(time, xinterp)
# plt0 = ax[0].pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
# fig.colorbar(plt0, ax=ax[0], orientation='horizontal')
# ax[0].set_ylim([time[t1], time[t2]])
# ax[0].set_title('Surveys (dev.)')
#
# plt1 = ax[1].pcolor(xg,tg,eofPred.T, vmin=-1.05, vmax=1.05)
# ax[1].set_ylim([time[t1], time[t2]])
# fig.colorbar(plt1, ax=ax[1], orientation='horizontal')
# ax[1].set_title('CEOF1 {:.2f}'.format(percentV[0]))
# ax[1].get_yaxis().set_ticks([])
#
# plt2 = ax[2].pcolor(xg,tg,eofPred2.T, vmin=-.85, vmax=.85)
# ax[2].set_ylim([time[t1], time[t2]])
# fig.colorbar(plt2, ax=ax[2], orientation='horizontal')
# ax[2].set_title('CEOF2 {:.2f}'.format(percentV[1]))
# ax[2].get_yaxis().set_ticks([])
#
# plt3 = ax[3].pcolor(xg,tg,eofPred3.T, vmin=-.65, vmax=.65)
# ax[3].set_ylim([time[t1], time[t2]])
# ax[3].set_title('CEOF3 {:.2f}'.format(percentV[2]))
# ax[3].get_yaxis().set_ticks([])
#
# fig.colorbar(plt3, ax=ax[3], orientation='horizontal')
#
# plt4 = ax[4].pcolor(xg,tg,eofPred4.T, vmin=-.65, vmax=.65)
# ax[4].set_ylim([time[t1], time[t2]])
# ax[4].set_title('CEOF4 {:.2f}'.format(percentV[3]))
# ax[4].get_yaxis().set_ticks([])
# fig.colorbar(plt4, ax=ax[4], orientation='horizontal')
#

#
# t1 = 200
# t2 = 400
# fig, ax = plt.subplots(1,5)
# #plt.set_cmap('RdBu')#bwr')
# plt.set_cmap('RdBu_r')
#
# tg, xg = np.meshgrid(time, xinterp)
# plt0 = ax[0].pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
# fig.colorbar(plt0, ax=ax[0], orientation='horizontal')
# ax[0].set_ylim([time[t1], time[t2]])
# ax[0].set_title('Surveys (dev.)')
#
# plt1 = ax[1].pcolor(xg,tg,eofPred.T, vmin=-1.05, vmax=1.05)
# ax[1].set_ylim([time[t1], time[t2]])
# fig.colorbar(plt1, ax=ax[1], orientation='horizontal')
# ax[1].set_title('CEOF1 {:.2f}'.format(percentV[0]))
# ax[1].get_yaxis().set_ticks([])
#
# plt2 = ax[2].pcolor(xg,tg,eofPred2.T, vmin=-.85, vmax=.85)
# ax[2].set_ylim([time[t1], time[t2]])
# fig.colorbar(plt2, ax=ax[2], orientation='horizontal')
# ax[2].set_title('CEOF2 {:.2f}'.format(percentV[1]))
# ax[2].get_yaxis().set_ticks([])
#
# plt3 = ax[3].pcolor(xg,tg,eofPred3.T, vmin=-.65, vmax=.65)
# ax[3].set_ylim([time[t1], time[t2]])
# ax[3].set_title('CEOF3 {:.2f}'.format(percentV[2]))
# ax[3].get_yaxis().set_ticks([])
#
# fig.colorbar(plt3, ax=ax[3], orientation='horizontal')
#
# plt4 = ax[4].pcolor(xg,tg,eofPred4.T, vmin=-.65, vmax=.65)
# ax[4].set_ylim([time[t1], time[t2]])
# ax[4].set_title('CEOF4 {:.2f}'.format(percentV[3]))
# ax[4].get_yaxis().set_ticks([])
# fig.colorbar(plt4, ax=ax[4], orientation='horizontal')

#plt.tight_layout(pad=0.5)

plt.show()


cpc1 = P2R(Rt[:, 0], phitSubset[:, 0])
cpc2 = P2R(Rt[:, 1], phitSubset[:, 1])
cpc3 = P2R(Rt[:, 2], phitSubset[:, 2])
cpc4 = P2R(Rt[:, 3], phitSubset[:, 3])
cpc5 = P2R(Rt[:, 4], phitSubset[:, 4])
cpc6 = P2R(Rt[:, 5], phitSubset[:, 5])
cpc7 = P2R(Rt[:, 6], phitSubset[:, 6])
cpc8 = P2R(Rt[:, 7], phitSubset[:, 7])

modes = 8

CPCs = np.vstack((np.real(cpc1), np.imag(cpc1), np.real(cpc2), np.imag(cpc2), np.real(cpc3), np.imag(cpc3),
                  np.real(cpc4), np.imag(cpc4), np.real(cpc5), np.imag(cpc5), np.real(cpc6), np.imag(cpc6),
                  np.real(cpc7), np.imag(cpc7), np.real(cpc8), np.imag(cpc8),))
CPCs = CPCs.T

var_explained = np.array((percentV[0], percentV[0], percentV[1], percentV[1], percentV[2], percentV[2],
                          percentV[3], percentV[3], percentV[4], percentV[4], percentV[5], percentV[5],
                          percentV[6], percentV[6], percentV[7], percentV[7]))

#
# fig = plt.figure()
# plt.plot(np.real(cpc1),np.imag(cpc1))
# plt.plot(np.real(cpc2),np.imag(cpc2))
# plt.plot(np.real(cpc3),np.imag(cpc3))
# plt.plot(np.real(cpc4),np.imag(cpc4))
# plt.plot(np.real(cpc5),np.imag(cpc5))
# plt.plot(np.real(cpc6),np.imag(cpc6))
# plt.plot(np.real(cpc7),np.imag(cpc7))
#
# fig = plt.figure()
# plt.plot(np.real(cpc1)*var_explained[0],np.imag(cpc1)*var_explained[0])
# plt.plot(np.real(cpc2)*var_explained[2],np.imag(cpc2)*var_explained[2])
# plt.plot(np.real(cpc3)*var_explained[4],np.imag(cpc3)*var_explained[4])
# plt.plot(np.real(cpc4)*var_explained[6],np.imag(cpc4)*var_explained[6])
# plt.plot(np.real(cpc5)*var_explained[8],np.imag(cpc5)*var_explained[8])

temp = CPCs*var_explained
# fig = plt.figure()
# plt.plot(temp[:,0], temp[:,1])
# plt.plot(temp[:,2], temp[:,3])
# plt.plot(temp[:,4], temp[:,5])
# plt.plot(temp[:,6], temp[:,7])
# plt.plot(temp[:,8], temp[:,9])


from sklearn.cluster import KMeans
numClusters = 12

PCsub = CPCs*var_explained

#  KMEANS
kma = KMeans(n_clusters=numClusters, n_init=2000).fit(PCsub)

# groupsize
_, group_size = np.unique(kma.labels_, return_counts=True)

# groups
d_groups = {}
for k in range(numClusters):
    d_groups['{0}'.format(k)] = np.where(kma.labels_ == k)

# # centroids
# centroids = np.dot(kma.cluster_centers_, EOFsub)
# centroids = kma.cluster_centers_
centroids = np.zeros((numClusters, PCsub.shape[1]))
for k in range(numClusters):
    print(PCsub[d_groups['{0}'.format(k)],0:1])
    centroids[k,:] = np.mean(PCsub[d_groups['{0}'.format(k)],:], axis=1)/var_explained

bmus = kma.labels_

    # # reorder clusters: bmus, km, cenEOFs, centroids, group_size
    # sorted_bmus = np.zeros((len(kma.labels_),),)*np.nan
    # for i in range(numClusters):
    #     posc = np.where(kma.labels_ == kma_order[i])
    #     sorted_bmus[posc] = i


from matplotlib import gridspec
# fig2 = plt.figure(figsize=(10,10))
# gs = gridspec.GridSpec(4, 5, wspace=0.0, hspace=0.0)
# gr, gc = 0, 0
profiles = np.zeros((numClusters,len(xinterp)))
magPhaseCentroid = np.zeros((np.shape(centroids)))

for i in range(numClusters):
    #getind = sorted[i]
    # cpc1RT, cpc1Phi = R2P(complex(centroids[i, 0], centroids[i, 1]))
    # cpc2RT, cpc2Phi = R2P(complex(centroids[i, 2], centroids[i, 3]))
    # cpc3RT, cpc3Phi = R2P(complex(centroids[i, 4], centroids[i, 5]))
    # cpc4RT, cpc4Phi = R2P(complex(centroids[i, 6], centroids[i, 7]))
    # cpc5RT, cpc5Phi = R2P(complex(centroids[i, 8], centroids[i, 9]))
    # cpc6RT, cpc6Phi = R2P(complex(centroids[i, 10], centroids[i, 11]))

    cpcRt = np.zeros((modes,))
    cpcPhi = np.zeros((modes,))

    cpcRt[0], cpcPhi[0] = R2P(complex(centroids[i, 0], centroids[i, 1]))
    cpcRt[1], cpcPhi[1] = R2P(complex(centroids[i, 2], centroids[i, 3]))
    cpcRt[2], cpcPhi[2] = R2P(complex(centroids[i, 4], centroids[i, 5]))
    cpcRt[3], cpcPhi[3] = R2P(complex(centroids[i, 6], centroids[i, 7]))
    cpcRt[4], cpcPhi[4] = R2P(complex(centroids[i, 8], centroids[i, 9]))
    cpcRt[5], cpcPhi[5] = R2P(complex(centroids[i, 10], centroids[i, 11]))

    cpcRt[6], cpcPhi[6] = R2P(complex(centroids[i, 8], centroids[i, 9]))
    cpcRt[7], cpcPhi[7] = R2P(complex(centroids[i, 10], centroids[i, 11]))

    magPhaseCentroid[i,0:16] = [cpcRt[0], cpcPhi[0], cpcRt[1], cpcPhi[1], cpcRt[2], cpcPhi[2],
                                cpcRt[3], cpcPhi[3], cpcRt[4], cpcPhi[4], cpcRt[5], cpcPhi[5],
                                cpcRt[6], cpcPhi[6], cpcRt[7], cpcPhi[7]]


    profile = 0 * np.ones(len(xinterp), )
    for mode in range(6):
        profile = profile + cpcRt[mode] * np.sin(cpcPhi[mode]) * S[:, mode] * np.sin(theta[:, mode]) + cpcRt[mode]* np.cos(cpcPhi[mode]) * S[:, mode] * np.cos(theta[:, mode])

    profiles[i,:] = profile+np.mean(alllines,axis=0)
    #profile = kma.centroids[i,:] * pred_std + pred_mean

    # ax = plt.subplot(gs[gr, gc])
    #
    # ax.plot(xinterp,profile+np.mean(alllines,axis=0))
    # ax.set_xlim([80, 820])
    # ax.set_ylim([-8, 2.25])
    # #ax.set_title('{}'.format(KMA.group_size.values[i]))
    # ax.text(400,0, group_size[i], fontweight='bold')
    #
    # if gc > 0:
    #     ax.set_yticks([])
    #
    # if gr < (5-1):
    #     ax.set_xticks([])
    # #  counter
    # gc += 1
    # if gc >= 5:
    #     gc = 0
    #     gr += 1


import scipy.signal as sig
deviation = np.zeros((np.shape(profiles)))

for i in range(numClusters):
    deviation[i,:] = profiles[i,:] - np.mean(alllines,axis=0)


magEOF1cen = magPhaseCentroid[:,0]
phaseEOF1cen = magPhaseCentroid[:,1]*180/np.pi


zeroShoreline = np.zeros((numClusters,))
offshorePeaks = np.zeros((numClusters,))
inshorePeaks = np.zeros((numClusters,))
inshoreDepth = np.zeros((numClusters,))
offshoreDepth = np.zeros((numClusters,))


# fig3 = plt.figure(figsize=(10,10))
# gs = gridspec.GridSpec(4, 5, wspace=0.0, hspace=0.0)
# gr, gc = 0, 0
for i in range(numClusters):
    #getind = sorted[i]
    dev = deviation[i,:]
    true = profiles[i,:]

    indexTesting = np.where((np.abs(true) == np.min(np.abs(true))))
    zeroShoreline[i] = indexTesting[0]
    peaks = sig.find_peaks(x=(true), prominence=0.01)

    if len(peaks[0]) > 0:
        offshorePeaks[i] = np.max(peaks[0])
        offshoreDepth[i] = true[np.max(peaks[0])]

    if len(peaks[0]) > 1:
        inshorePeaks[i] = np.min(peaks[0])
        inshoreDepth[i] = true[np.min(peaks[0])]

    # ax = plt.subplot(gs[gr, gc])
    #
    # ax.plot(xinterp,true)
    # ax.plot(xinterp[peaks[0]],true[peaks[0]],'ro')
    # ax.set_xlim([80, 820])
    # ax.set_ylim([-8, 2.25])
    # #ax.set_title('{}'.format(KMA.group_size.values[i]))
    # ax.text(400,0, group_size[i], fontweight='bold')
    # #ax.text(400,0, phaseEOF1cen[i], fontweight='bold')
    #
    # if gc > 0:
    #     ax.set_yticks([])
    #
    # if gr < (5-1):
    #     ax.set_xticks([])
    # #  counter
    # gc += 1
    # if gc >= 5:
    #     gc = 0
    #     gr += 1






# # totalDepths = offshoreDepth+inshoreDepth
# shorelinePeakDiff = offshorePeaks-zeroShoreline
# sortedPeaks = np.sort(shorelinePeakDiff)
# sortedPeakInd = np.argsort(shorelinePeakDiff)

sortedPeaks = np.sort(phaseEOF1cen)
sortedPeakInd = np.argsort(phaseEOF1cen)

sortedPhases = phaseEOF1cen[sortedPeakInd]

# doubleIndex = np.where((inshorePeaks > 0))
# singleIndex = np.where((inshorePeaks == 0))
# sortedInPeaks = np.argsort(offshorePeaks[doubleIndex[0]])
# sortedOffPeaks = np.argsort(offshorePeaks[singleIndex[0]])
# sortedPeakInd = np.append(singleIndex[0][sortedOffPeaks],doubleIndex[0][sortedInPeaks])

#
# fig3 = plt.figure(figsize=(10,10))
# gs = gridspec.GridSpec(numClusters, 1, wspace=0.0, hspace=0.0)
# gr, gc = 0, 0
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, numClusters))
#
# for i in range(numClusters):
#     #getind = sorted[i]
#     dev = deviation[sortedPeakInd[i],:]
#     true = profiles[sortedPeakInd[i],:]
#
#     peaks = sig.find_peaks(x=(true), prominence=0.05)
#
#     #if len(peaks[0]) > 0:
#     #    offshorePeaks[i] = np.max(peaks[0])
#
#     ax = plt.subplot(gs[gr, gc])
#
#     ax.plot(xinterp,true,color=colors[i])
#     #ax.plot(xinterp[peaks[0]],true[peaks[0]],'ro')
#     ax.set_xlim([80, 720])
#     ax.set_ylim([-8, 2.25])
#     #ax.set_title('{}'.format(KMA.group_size.values[i]))
#     #ax.text(400,0, group_size[sortedPeakInd[i]], fontweight='bold')
#     ax.text(400,0, sortedPhases[i], fontweight='bold')
#
#     if gc > 0:
#         ax.set_yticks([])
#
#     if gr < (numClusters-1):
#         ax.set_xticks([])
#     #  counter
#     gr += 1
#     if gr >= numClusters:
#         gc += 1
#         gr = 0

#
# fig3 = plt.figure(figsize=(10,10))
# gs = gridspec.GridSpec(1, numClusters, wspace=0.0, hspace=0.0)
# gr, gc = 0, 0
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, numClusters))
#
# for i in range(numClusters):
#     #getind = sorted[i]
#     dev = deviation[sortedPeakInd[i],:]
#     true = profiles[sortedPeakInd[i],:]
#
#     peaks = sig.find_peaks(x=(true), prominence=0.05)
#
#     #if len(peaks[0]) > 0:
#     #    offshorePeaks[i] = np.max(peaks[0])
#
#     ax = plt.subplot(gs[gr, gc])
#     ax.plot(xinterp,np.mean(alllines,axis=0),'w-')
#
#     ax.plot(xinterp,true,color=colors[i])
#     #ax.plot(xinterp[peaks[0]],true[peaks[0]],'ro')
#     ax.set_xlim([80, 650])
#     ax.set_ylim([-7, 1])
#     #ax.set_title('{}'.format(KMA.group_size.values[i]))
#     #ax.text(400,0, group_size[sortedPeakInd[i]], fontweight='bold')
#     ax.text(400,0, sortedPhases[i], fontweight='bold')
#
#     if gr > 0:
#         ax.set_yticks([])
#
#     # if gc < numClusters:
#     #     ax.set_xticks([])
#     #  counter
#     gc += 1
#     if gc >= numClusters:
#         gr += 1
#         gc = 0


bmu = np.zeros((len(time),))

for i in range(numClusters):
    bmuIndex = np.where((bmus == sortedPeakInd[i]))
    bmu[bmuIndex] = i
#
# fig5 = plt.figure(figsize=(10,5))
# plt.plot(bmu,time,'o')
# plt.show()
#
# fig15 = plt.figure()
# plt.plot(centroids[:,0],centroids[:,1],'wo')
# plt.scatter(CPCs[:,0],CPCs[:,1],5,bmu,cmap='gist_rainbow')
# plt.show()
#
# fig16 = plt.figure()
# plt.plot(centroids[:,2],centroids[:,3],'wo')
# plt.scatter(CPCs[:,2],CPCs[:,3],5,bmu,cmap='gist_rainbow')
# plt.show()
#
# fig6 = plt.figure(figsize=(10,10))
# gs = gridspec.GridSpec(4, 5, wspace=0.0, hspace=0.0)
# gr, gc = 0, 0
# for i in range(numClusters):
#     #getind = sorted[i]
#     bmuIndex = np.where((bmus == sortedPeakInd[i]))
#
#     subset = alllines[bmuIndex,:]
#
#     ax = plt.subplot(gs[gr, gc])
#     for x in range(len(bmuIndex[0])):
#         ax.plot(xinterp,subset[0,x,:])
#
#     ax.set_xlim([80, 820])
#     ax.set_ylim([-8, 2.25])
#     #ax.set_title('{}'.format(KMA.group_size.values[i]))
#     ax.text(400,0, group_size[sortedPeakInd[i]], fontweight='bold')
#
#     if gc > 0:
#         ax.set_yticks([])
#
#     if gr < (5-1):
#         ax.set_xticks([])
#     #  counter
#     gc += 1
#     if gc >= 5:
#         gc = 0
#         gr += 1
#


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


# fig2 = plt.figure(figsize=(8,8))
# plt.plot(time,innerBar,'bo')
# plt.plot(time,outerBar,'ro')


#
# fig11, ax11 = plt.subplots(2,1)
# #ax10 = fig.add_subplot(111)
# #ax10.set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, numClusters))))
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, numClusters))
colors2 = cm.gray(np.linspace(0, 1, numClusters))
# for i in range(numClusters):
#     bmuIndex = np.where((bmus == sortedPeakInd[i]))
#     ax11[0].scatter(time[bmuIndex], outerBar[bmuIndex], label=time[i], color = colors[i])
#     ax11[1].scatter(time[bmuIndex], innerBar[bmuIndex], color = colors[i])
#
# ax11[0].set_ylabel('xFRF')
# ax11[1].set_ylabel('xFRF')



#
#
# fig3 = plt.figure(figsize=(10,10))
# gs = gridspec.GridSpec(10, 10, wspace=0.0, hspace=0.0)
# gr, gc = 0, 0
# import matplotlib.cm as cm
# colors = cm.rainbow(np.linspace(0, 1, numClusters))
#
# for i in range(numClusters):
#     #getind = sorted[i]
#     dev = deviation[sortedPeakInd[i],:]
#     true = profiles[sortedPeakInd[i],:]
#
#     peaks = sig.find_peaks(x=(true), prominence=0.05)
#
#     #if len(peaks[0]) > 0:
#     #    offshorePeaks[i] = np.max(peaks[0])
#
#     ax = plt.subplot(gs[gr, gc])
#
#     ax.plot(xinterp,true,color=colors[i])
#     #ax.plot(xinterp[peaks[0]],true[peaks[0]],'ro')
#     ax.set_xlim([80, 720])
#     ax.set_ylim([-8, 2.25])
#     #ax.set_title('{}'.format(KMA.group_size.values[i]))
#     ax.text(400,0, group_size[i], fontweight='bold')
#
#     if gc > 0:
#         ax.set_yticks([])
#
#     if gr < (10-1):
#         ax.set_xticks([])
#     #  counter
#     gr += 1
#     if gr >= 10:
#         gc += 1
#         gr = 0
#
#
# ax15 = plt.subplot(gs[0:5, 3:])
# ax16 = plt.subplot(gs[6:,3:])
# for i in range(numClusters):
#     bmuIndex = np.where((bmus == sortedPeakInd[i]))
#     ax15.scatter(time[bmuIndex], outerBar[bmuIndex], label=time[i], color = colors[i])
#     ax16.scatter(time[bmuIndex], innerBar[bmuIndex], color = colors[i])
#
# ax15.set_ylabel('xFRF')
# ax16.set_ylabel('xFRF')
# ax15.set_title('Offshore Bar')
# ax16.set_title('Onshore Bar')

sorted_bmus = np.tile(0,(len(kma.labels_),), )
sorted_time = np.tile(0,(len(kma.labels_),), )

for i in range(numClusters):
    posc = np.where(kma.labels_ == sortedPeakInd[i])
    sorted_bmus[posc] = int(i)
    #sorted_time[posc] = time[posc]

# fig5 = plt.figure(figsize=(10,5))
# plt.plot(sorted_bmus,time,'o')
# plt.show()

def transition_matrix(transitions):
    n = 1+ max(transitions) #number of states
    M = [[0]*n for _ in range(n)]
    M2 = [[0]*n for _ in range(n)]
    for (i,j) in zip(transitions,transitions[1:]):
        M[i][j] += 1
        M2[i][j] += 1
    #now convert to probabilities:
    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    return M, M2

# def transitionDictionary(bmus,dates):

bins = dict()
date = dict()
nextbin = dict()
nextdate = dict()
prevbin = dict()
prevdate = dict()
for xx in range(numClusters):
    binind = np.where((sorted_bmus==xx))
    bins[xx] = binind
    date[xx] = time[binind]
    nextbinind = np.where((sorted_bmus[0:-1]==xx))
    K = 1
    res = [x + K for x in nextbinind]
    nextbin[xx] = sorted_bmus[res]
    nextdate[xx] = time[res]
    prevbinind = np.where((sorted_bmus==xx))
    res2 = [x - K for x in prevbinind]
    prevbin[xx] = sorted_bmus[res2]
    prevdate[xx] = time[res2]



m,m2 = transition_matrix(sorted_bmus)
for row in m: print(' '.join('{0:.2f}'.format(x) for x in row))

flat_list = [item for sublist in m for item in sublist]
flatarray = np.asarray(flat_list)
flatarray.resize(numClusters, numClusters)

flat_list2 = [item for sublist in m2 for item in sublist]
flatarray2 = np.asarray(flat_list2)
flatarray2.resize(numClusters, numClusters)

# fig = plt.figure(figsize=(10,10))
# ax = plt.subplot2grid((2,2),(0,0),colspan=2,rowspan=2)
# ax.imshow(flatarray,vmin=0,vmax=0.3,cmap='Oranges')
# plt.show()
#
# import h5py
#
# waveFile = '/home/dylananderson/projects/duckGeomorph/8m_edat.mat'
# # waveFile = '/home/dylananderson/projects/duckGeomorph/17m_2D_edat.mat'
#
# waves = dict()
# with h5py.File(waveFile, 'r') as f:
#     for k in f.keys():
#         print(k)
#     waves['hs'] = f['edat/hs'][:]
#     waves['tp'] = f['edat/tp'][:]
#     waves['dm'] = f['edat/dm'][:]
#     waves['time'] = f['edat/time'][:]
#
# from datetime import datetime
# from datetime import timedelta
# def datenum_to_datetime(datenum):
#     """
#     Convert Matlab datenum into Python datetime.
#     :param datenum: Date in datenum format
#     :return:        Datetime object corresponding to datenum.
#     """
#     days = datenum % 1
#     hours = days % 1 * 24
#     minutes = hours % 1 * 60
#     seconds = minutes % 1 * 60
#     return datetime.fromordinal(int(datenum)) \
#            + timedelta(days=int(days)) \
#            + timedelta(hours=int(hours)) \
#            + timedelta(minutes=int(minutes)) \
#            + timedelta(seconds=round(seconds)) \
#            - timedelta(days=366)
#
# timeWave = [datenum_to_datetime(x) for x in waves['time'].flatten()]
#
# # python_datetime = datetime.fromordinal(int(matlab_datenum)) + timedelta(days=matlab_datenum%1) - timedelta(days = 366)
#
#
#
# def getArray(file):
#     waves = Dataset(file)
#
#     waveHs = waves.variables['waveHs'][:]
#     waveTp = waves.variables['waveTp'][:]
#     #waveTm = waves.variables['waveTm'][:]
#     #waveTm1 = waves.variables['waveTm1'][:]
#     #waveTm2 = waves.variables['waveTm2'][:]
#     waveMeanDirection = waves.variables['waveMeanDirection'][:]
#     #waveMeanDirectionSwell = waves.variables['waveMeanDirectionSwell'][:]
#     timeW = waves.variables['time'][:]
#
#
#     output = dict()
#     output['waveHs'] = waveHs
#     output['waveTp'] = waveTp
#     #output['waveTm'] = waveTm
#     #output['waveTm1'] = waveTm1
#     #output['waveTm2'] = waveTm2
#     output['waveMeanDirection'] = waveMeanDirection
#     #output['waveMeanDirectionSwell'] = waveMeanDirectionSwell
#     output['t'] = timeW
#
#     return output
#
# wavedir = '/media/dylananderson/Elements/8mArray/'
# # Need to sort the files to ensure correct temporal order...
# files = os.listdir(wavedir)
# files.sort()
# files_path = [os.path.join(os.path.abspath(wavedir), x) for x in files]
# array8m = Dataset(files_path[0])
# Hs8m = []
# Tp8m = []
# Dm8m = []
# timeWave8m = []
# for i in files_path:
#     waves8m = getArray(i)
#     Hs8m = np.append(Hs8m,waves8m['waveHs'])
#     Tp8m = np.append(Tp8m,waves8m['waveTp'])
#     Dm8m = np.append(Dm8m,waves8m['waveMeanDirection'])
#     timeWave8m = np.append(timeWave8m,waves8m['t'])
#
# ind = np.where((timeWave8m > 1000000))
# hs8m = Hs8m[ind]
# tp8m = Tp8m[ind]
# dm8m = Dm8m[ind]
# t8m = timeWave8m[ind]
# tWave8m = [DT.datetime.fromtimestamp(x) for x in t8m]
#
#
#
# wavedir = '/media/dylananderson/Elements/17mArray/'
# # Need to sort the files to ensure correct temporal order...
# files = os.listdir(wavedir)
# files.sort()
# files_path = [os.path.join(os.path.abspath(wavedir), x) for x in files]
# array17m = Dataset(files_path[0])
# Hs17m = []
# Tp17m = []
# Dm17m = []
# timeWave17m = []
# for i in files_path:
#     waves17m = getArray(i)
#     Hs17m = np.append(Hs17m,waves17m['waveHs'])
#     Tp17m = np.append(Tp17m,waves17m['waveTp'])
#     Dm17m = np.append(Dm17m,waves17m['waveMeanDirection'])
#     timeWave17m = np.append(timeWave17m,waves17m['t'])
#
# ind = np.where((Hs17m > 0))
# hs17m = Hs17m[ind]
# tp17m = Tp17m[ind]
# dm17m = Dm17m[ind]+1000
# t17m = timeWave17m[ind]
# tWave17m = [DT.datetime.fromtimestamp(x) for x in t17m]
#
#
# # plt.figure(figsize=(10,10))
# # ax = plt.subplot2grid((3,3,),(0,0),rowspan=1,colspan=3)
# # ax.plot(timeWave,waves['hs'],label='reconstruction')
# # ax.plot(tWave17m,hs17m,label='17mArray')
# # ax.plot(tWave8m,hs8m,label='8mArray')
# # plt.legend()
# # ax2 = plt.subplot2grid((3,3,),(1,0),rowspan=1,colspan=3)
# # ax2.plot(timeWave,waves['tp'])
# # ax2.plot(tWave17m,tp17m)
# # ax2.plot(tWave8m,tp8m)
# # ax3 = plt.subplot2grid((3,3,),(2,0),rowspan=1,colspan=3)
# # ax3.plot(timeWave,waves['dm'])
# # ax3.plot(tWave17m,dm17m)
# # ax3.plot(tWave8m,dm8m)
# #
# # plt.show()
#
# HsArrays = np.append(hs17m,hs8m)
# HsCombined = np.append(waves['hs'], HsArrays)
# TpArrays = np.append(tp17m,tp8m)
# TpCombined = np.append(waves['tp'], TpArrays)
# DmArrays = np.append(dm17m,dm8m)
# DmCombined = np.append(waves['dm'], DmArrays)
# TimeArrays = np.append(tWave17m,tWave8m)
# TimeCombined = np.append(timeWave, TimeArrays)
#
# tC = TimeCombined[np.argsort(TimeCombined)]
# hsC = HsCombined[np.argsort(TimeCombined)]
# tpC = TpCombined[np.argsort(TimeCombined)]
# dmC = DmCombined[np.argsort(TimeCombined)]
#
# badDirs = np.where((dmC > 360))
# dmC[badDirs] = dmC[badDirs]*np.nan
#
# waveNorm = dmC - 72
# neg = np.where((waveNorm > 180))
# waveNorm[neg[0]] = waveNorm[neg[0]]-360
# offpos = np.where((waveNorm>90))
# offneg = np.where((waveNorm<-90))
# waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
# waveNorm[offneg[0]] = waveNorm[offneg[0]]*0
#
# lwpC = 1025*np.square(hsC)*tpC*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
# weC = np.square(hsC)*tpC
# # tWave = [DT.datetime.fromtimestamp(x) for x in timeWave]
#
# # plt.figure(figsize=(10,10))
# # ax = plt.subplot2grid((3,3,),(0,0),rowspan=1,colspan=3)
# # ax.plot(tC,hsC,'-',label='Combined')
# # plt.legend()
# # ax2 = plt.subplot2grid((3,3,),(1,0),rowspan=1,colspan=3)
# # ax2.plot(tC,tpC,'-')
# # ax3 = plt.subplot2grid((3,3,),(2,0),rowspan=1,colspan=3)
# # ax3.plot(tC,dmC,'.')
# # ax3.set_ylim([0, 180])
# # plt.show()
#
# # overlap = np.in1d(tWave17m,tWave8m) # length of datetime1, gives True for those elements which exist in datetime2
# # ntWave17m = np.delete(tWave17m,overlap)


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

from datetime import datetime
from datetime import timedelta
def datenum_to_datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    days = datenum % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60
    return datetime.fromordinal(int(datenum)) \
           + timedelta(days=int(days)) \
           + timedelta(hours=int(hours)) \
           + timedelta(minutes=int(minutes)) \
           + timedelta(seconds=round(seconds)) \
           - timedelta(days=366)



Hs = []
Tp = []
Dm = []
timeWave = []
for i in files_path:
    waves = getWIS(i)
    Hs = np.append(Hs,waves['waveHs'])
    Tp = np.append(Tp,waves['waveTp'])
    Dm = np.append(Dm,waves['waveMeanDirection'])
    #timeTemp = [datenum_to_datetime(x) for x in waves['t'].flatten()]
    timeWave = np.append(timeWave,waves['t'].flatten())


hsC = Hs
tpC = Tp
dmC = Dm
badDirs = np.where((dmC > 360))
dmC[badDirs] = dmC[badDirs]*np.nan

waveNorm = dmC - 72
neg = np.where((waveNorm > 180))
waveNorm[neg[0]] = waveNorm[neg[0]]-360
offpos = np.where((waveNorm>90))
offneg = np.where((waveNorm<-90))
waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
waveNorm[offneg[0]] = waveNorm[offneg[0]]*0

lwpC = 1025*np.square(hsC)*tpC*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
weC = np.square(hsC)*tpC

tWave = [DT.datetime.fromtimestamp(x) for x in timeWave]

tC = np.array(tWave)








# bins = dict()
# date = dict()
# nextbin = dict()
# nextdate = dict()
# prevbin = dict()
# prevdate = dict()
wHs = []
wTp = []
wDm = []
wLWP = []
wWE = []
wT = []
for xx in range(numClusters):
    innerListHs = []
    innerListTp = []
    innerListDm = []
    innerListLWP = []
    innerListWE = []
    innerListT = []
    for yy in range(numClusters):
        # if both are equal then the wave conditions didn't cause a transition and
        # everything between the last two dates should be added to a distribution

        # if yy is different...
        wInd = np.where(prevbin[xx]==yy)
        if len(wInd[0])>0:
            tempHs = []
            tempTp = []
            tempDm = []
            tempLWP = []
            tempWE = []
            tempTi = []
            for tt in range(len(wInd[0])):
                tempT = np.where((tC < date[xx][wInd[0][tt]]) & (tC > prevdate[xx][wInd[0][tt]]))
                tempHs = np.append(tempHs,hsC[tempT])
                tempTp = np.append(tempTp,tpC[tempT])
                tempDm = np.append(tempDm,waveNorm[tempT])
                tempLWP = np.append(tempLWP,lwpC[tempT])
                tempWE = np.append(tempWE,weC[tempT])
                tempTi = np.append(tempTi,tC[tempT])

        else:
            tempHs = []
            tempTp = []
            tempDm = []
            tempLWP =[]
            tempWE= []
            tempTi = []

        innerListHs.append(tempHs)
        innerListTp.append(tempTp)
        innerListDm.append(tempDm)
        innerListLWP.append(tempLWP)
        innerListWE.append(tempWE)
        innerListT.append(tempTi)

    wHs.append(innerListHs)
    wTp.append(innerListTp)
    wDm.append(innerListDm)
    wLWP.append(innerListLWP)
    wWE.append(innerListWE)
    wT.append(innerListT)



# Ok, first lets load a DWT struct

import scipy.io as sio
from scipy.io.matlab.mio5_params import mat_struct
import numpy as np

def ReadMatfile(p_mfile):
    'Parse .mat file to nested python dictionaries'

    def RecursiveMatExplorer(mstruct_data):
        # Recursive function to extrat mat_struct nested contents

        if isinstance(mstruct_data, mat_struct):
            # mstruct_data is a matlab structure object, go deeper
            d_rc = {}
            for fn in mstruct_data._fieldnames:
                d_rc[fn] = RecursiveMatExplorer(getattr(mstruct_data, fn))
            return d_rc

        else:
            # mstruct_data is a numpy.ndarray, return value
            return mstruct_data

    # base matlab data will be in a dict
    mdata = sio.loadmat(p_mfile, squeeze_me=True, struct_as_record=False)
    mdata_keys = [x for x in mdata.keys() if x not in
                  ['__header__','__version__','__globals__']]

    # use recursive function
    dout = {}
    for k in mdata_keys:
        dout[k] = RecursiveMatExplorer(mdata[k])
    return dout


DWT = ReadMatfile('/media/dylananderson/Elements1/NC_climate/Nags_Head_DWTS_25_w50minDates_plus5TCs.mat')

PCA = ReadMatfile('/media/dylananderson/Elements1/NC_climate/Nags_Head_SLPS_1degree_memory.mat')

mycols = ReadMatfile('/media/dylananderson/Elements1/shusin6_contents/codes/mycmap_col.mat')
mycmap = mycols['mycmap_col']
# Need to get the dates for the bmus into the correct format (datetime)
def datevec2datetime(d_vec):
    '''
    Returns datetime list from a datevec matrix
    d_vec = [[y1 m1 d1 H1 M1],[y2 ,2 d2 H2 M2],..]
    '''
    return [datetime(d[0], d[1], d[2], d[3], d[4]) for d in d_vec]

dwtDates = np.array(datevec2datetime(DWT['DWT']['dates']))
dwtBmus = DWT['DWT']['bmus']
removeTooSoon = np.where(dwtDates < tC[-45])

dwtTimes = dwtDates[removeTooSoon]
dwtBMUS = dwtBmus[removeTooSoon]
# Next, we need to split the waves up into each DWT
numDWTs = 30
dwtHs = []
dwtMaxHs = []
dwtTp = []
dwtDm = []
dwtLWP = []
dwtWE = []
dwtT = []
for xx in range(numDWTs):
    tempHs = []
    tempMaxHs = []
    tempTp = []
    tempDm = []
    tempLWP = []
    tempWE = []
    tempTi = []

    wInd = np.where((dwtBMUS[0:-1] == (xx+1)))
    for tt in range(len(wInd[0])):
        tempT = np.where((tC < dwtTimes[wInd[0][tt]+1]) & (tC > dwtTimes[wInd[0][tt]]))
        if len(tempT[0]) > 0:
            tempMaxHs = np.append(tempMaxHs, np.nanmax(hsC[tempT]))

        tempHs = np.append(tempHs, hsC[tempT])
        tempTp = np.append(tempTp, tpC[tempT])
        tempDm = np.append(tempDm, waveNorm[tempT])
        tempLWP = np.append(tempLWP, lwpC[tempT])
        tempWE = np.append(tempWE, weC[tempT])
        tempTi = np.append(tempTi, tC[tempT])

    dwtHs.append(tempHs)
    dwtMaxHs.append(tempMaxHs)
    dwtTp.append(tempTp)
    dwtDm.append(tempDm)
    dwtLWP.append(tempLWP)
    dwtWE.append(tempWE)
    dwtT.append(tempTi)




order = DWT['DWT']['order']
meanDWTHs = np.zeros((np.shape(order)))
for xx in range(numDWTs):
    data = dwtMaxHs[xx]
    meanDWTHs[xx] = np.nanmean(data)

newOrder = np.argsort(meanDWTHs)

orderedDWTs = np.zeros((np.shape(dwtBMUS)))
for xx in range(numDWTs):
    dateInd = np.where((dwtBMUS == xx))
    orderedDWTs[dateInd] = newOrder[xx]

etcolors = cm.rainbow(np.linspace(0, 1, numDWTs-5))
tccolors = np.flipud(cm.gray(np.linspace(0,1,6)))

dwtcolors = np.vstack((etcolors,tccolors[1:,:]))


import matplotlib.cm as cm
import matplotlib.colors as mcolors
plt.style.use('dark_background')

from scipy.stats.kde import gaussian_kde
dist_space = np.linspace(0, 4, 80)
fig = plt.figure(figsize=(10,10))
colorparam = np.zeros((numDWTs,))
counter = 0
plotIndx = 0
plotIndy = 0
for xx in range(numDWTs):

    dwtInd = order[xx]-1
    #dwtInd = newOrder[xx]

    ax = plt.subplot2grid((6, 5), (plotIndx, plotIndy), rowspan=1, colspan=1)
    # normalize = mcolors.Normalize(vmin=np.min(colorparams), vmax=np.max(colorparams))
    normalize = mcolors.Normalize(vmin=.5, vmax=1.8)

    ax.set_xlim([0, 3])
    ax.set_ylim([0, 2])
    data = dwtHs[dwtInd]
    if len(data) > 0:
        kde = gaussian_kde(data)
        colorparam[counter] = np.nanmean(data)
        colormap = cm.Reds
        color = colormap(normalize(colorparam[counter]))
        ax.plot(dist_space, kde(dist_space), linewidth=1, color=color)
        ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
        ax.spines['top'].set_color([0.5, 0.5, 0.5])
        ax.spines['right'].set_color([0.5, 0.5, 0.5])
        ax.spines['left'].set_color([0.5, 0.5, 0.5])
        # ax.text(1.8, 1, np.round(colorparam*100)/100, fontweight='bold')

    else:
        ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
        ax.spines['top'].set_color([0.3, 0.3, 0.3])
        ax.spines['right'].set_color([0.3, 0.3, 0.3])
        ax.spines['left'].set_color([0.3, 0.3, 0.3])

    if plotIndx < 5:
        ax.xaxis.set_ticks([])
        ax.xaxis.set_ticklabels([])
    if plotIndy > 0:
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
    counter = counter + 1
    if plotIndy < 4:
        plotIndy = plotIndy + 1
    else:
        plotIndy = 0
        plotIndx = plotIndx + 1
    print(plotIndy, plotIndx)

plt.show()
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map.set_array(colorparam)
fig.subplots_adjust(right=0.86)
cbar_ax = fig.add_axes([0.89, 0.15, 0.02, 0.7])
cbar = fig.colorbar(s_map, cax=cbar_ax)
cbar.set_label('Mean Hs (m)')



repmatDesviacion = np.tile(DWT['PCA']['Desviacion'], (25,1))
repmatMedia = np.tile(DWT['PCA']['Media'], (25,1))
Km_ = np.multiply(DWT['KMA']['centroids'],repmatDesviacion) + repmatMedia

[mK, nK] = np.shape(Km_)

Km_slp = Km_[:,0:int(nK/2)]
Km_grd = Km_[:,int(nK/2):]
X_B = DWT['X_B']
Y_B = DWT['Y_B']
SLP_C = DWT['SLP_C']
Km_slp = Km_slp[:,0:len(X_B)]
Km_grd = Km_grd[:,1:len(X_B)]

Xs = np.arange(np.min(X_B),np.max(X_B))
#Xs = np.arange(-(360-np.min(X_B)),(np.max(X_B)-360))
Ys = np.arange(np.min(Y_B),np.max(Y_B))
lenXB = len(X_B)
[XR,YR] = np.meshgrid(Xs,Ys)
#sea_nodes = np.arange(0,len(X_B))
#sea_nodes = np.zeros((np.shape(X_B)))
sea_nodes = []
for qq in range(lenXB-1):
    sea_nodes.append(np.where((XR == X_B[qq]) & (YR == Y_B[qq])))
    #print(indexTest[0])

flat_list = [item for sublist in sea_nodes for item in sublist]


plt.style.use('default')
from mpl_toolkits.basemap import Basemap
fig2 = plt.figure(figsize=(8,10))
gs1 = gridspec.GridSpec(6, 5)
gs1.update(wspace=0.00, hspace=0.00) # set the spacing between axes.
counter = 0
plotIndx = 0
plotIndy = 0
for qq in range(numDWTs):
#qq = 0

    #ax = plt.subplot2grid((7, 7), (plotIndx, plotIndy), rowspan=1, colspan=1)
    ax = plt.subplot(gs1[qq])
    #num = order[qq]
    num = newOrder[qq]
    wt = np.ones((np.shape(XR))) * np.nan

    if qq < 25:
        num_index = np.where((dwtBMUS == num))
        temp = Km_slp[(num-1),:]/100 - np.nanmean(SLP_C,axis=0)/100
    else:
        num_index = np.where((dwtBmus == num))
        temp = np.nanmean(SLP_C[num_index,:],axis=1)/100 - np.nanmean(SLP_C,axis=0)/100
        temp= temp.flatten()
    for tt in range(len(sea_nodes)):
        wt[sea_nodes[tt]] = temp[tt]
#reshaped = np.tile(wt,(np.shape(XR)))


#m = Basemap(projection='cyl',lon_0=320,lat_0=0)
    m = Basemap(projection='merc',llcrnrlat=-40,urcrnrlat=50,llcrnrlon=275,urcrnrlon=370,lat_ts=10,resolution='c')
#m = Basemap(projection='stere',lon_0=-120,lat_0=20,llcrnrlat=-10,llcrnrlon=50,urcrnrlat=50, urcrnrlon=190)
#m = Basemap(width=12000000,height=8000000,resolution='l',projection='stere',lat_ts=30,lat_0=20,lon_0=-70.)
    m.fillcontinents(color=dwtcolors[qq])
    m.drawcoastlines()   #畫海岸線
#parallels = np.arange(-90.,90,30.)
#meridians = np.arange(0.,360.,20.)
#m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)  # 緯度度線、在左右兩邊加緯度標籤
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    clevels = np.arange(-17,17,1)
    cx,cy =m(XR,YR)  # convert to map projection coordinate
#CS = m.pcolormesh(cx,cy,wt,clevels,cmap=cm.jet,shading='gouraud')
    CS = m.contourf(cx,cy,wt,clevels,vmin=-7.5,vmax=7.5,cmap=cm.RdBu_r,shading='gouraud')

    #plt.colorbar(CS,orientation='horizontal')
    if plotIndx < 6:
        ax.xaxis.set_ticks([])
        ax.xaxis.set_ticklabels([])
    if plotIndy > 0:
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
    counter = counter + 1
    if plotIndy < 6:
        plotIndy = plotIndy + 1
    else:
        plotIndy = 0
        plotIndx = plotIndx + 1
#date=time_d[0].strftime('%Y/%m/%d')
#plt.title('Cylindrical, T at 1000 hPa, GFS'+date)
#plt.tight_layout()
colormap = cm.RdBu_r
normalize = mcolors.Normalize(vmin=-7.5, vmax=7.5)

s_map2 = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map2.set_array(colorparam)
fig.subplots_adjust(right=0.85)
cbar_ax2 = fig2.add_axes([0.91, 0.15, 0.02, 0.7])
cbar2 = fig2.colorbar(s_map2, cax=cbar_ax2)
cbar2.set_label('slp anom (mbar)')

plt.show()




### separate the DWTs into each transitions

transitionDWTs = []

for xx in range(numClusters):
    innerListDWT = []

    for yy in range(numClusters):
        # if both are equal then the wave conditions didn't cause a transition and
        # everything between the last two dates should be added to a distribution

        # if yy is different...
        wInd = np.where(prevbin[xx]==yy)
        if len(wInd[0])>0:
            tempDWT = []

            for tt in range(len(wInd[0])):
                tempT = np.where((dwtTimes < date[xx][wInd[0][tt]]) & (dwtTimes > prevdate[xx][wInd[0][tt]]))

                #tempDWT = np.append(tempDWT,dwtBMUS[tempT])
                tempDWT = np.append(tempDWT,orderedDWTs[tempT])

        else:
            tempDWT = []

        innerListDWT.append(tempDWT)

    transitionDWTs.append(innerListDWT)




def survey(results, category_names,category_colors,ax):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)


    #fig, ax = plt.subplots(figsize=(9.2, 5))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5,
                label=colname, color=color)
        xcenters = starts + widths / 2

        r, g, b, _ = color
        text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        # for y, (x, c) in enumerate(zip(xcenters, widths)):
        #     ax.text(x, y, str(int(c)), ha='center', va='center',
        #             color=text_color)
    # ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
    #           loc='lower left', fontsize='small')

    return ax
#
# import cartopy.crs as ccrs
#
# fig = plt.figure(figsize=(10, 5))
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
#
# lons = XR
# lats = YR
# data = wt
#
# ax.contourf(lons, lats, data,
#                 transform=ccrs.PlateCarree(),
#                 cmap='nipy_spectral')
# ax.coastlines()
# ax.set_global()
# plt.show()

# t1 = transitionDWTs[0][0]
# t1.sort()
#
# uni = np.unique(t1)
# uniInt = [int(num) for num in uni]
# results = dict()
# categories = list()
# results[str(uni[0])] = list()
# for qq in range(len(uni)):
#     categories.append(str(uni[qq]))
#     results[str(uni[0])].append(len((np.where((t1 == uni[qq])))[0]))
#
# category_colors = dwtcolors[uniInt,:]
#
# survey(results, categories,category_colors)


import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.stats.kde import gaussian_kde
dist_space = np.linspace(0, 5, 50)
fig10 = plt.figure(figsize=(10,10))
colorparam = np.zeros((numClusters*numClusters,))
counter = 0
for xx in range(numClusters):
    for yy in range(numClusters):
        ax = plt.subplot2grid((numClusters,numClusters), (yy,xx), rowspan=1, colspan=1)
        #normalize = mcolors.Normalize(vmin=np.min(colorparams), vmax=np.max(colorparams))
        #normalize = mcolors.Normalize(vmin=.5, vmax=1.5)

        #ax.set_xlim([0,3])
        #ax.set_ylim([0,2])
        data = transitionDWTs[xx][yy]
        data.sort()
        if len(data)>0:
            uni = np.unique(data)
            uniInt = [int(num) for num in uni]
            results = dict()
            categories = list()
            #results[str(uni[0])] = list()
            results['0'] = list()

            for qq in range(len(uni)):
                categories.append(str(uni[qq]))
                # results[str(uni[0])].append(len((np.where((t1 == uni[qq])))[0]))
                results['0'].append(len((np.where((data == uni[qq])))[0]))

            category_colors = dwtcolors[uniInt, :]

            survey(results, categories, category_colors,ax)

        # if len(data)>0:
        #     kde = gaussian_kde(data)
        #     colorparam[counter] = np.nanmean(data)
        #     colormap = cm.Reds
        #     color = colormap(normalize(colorparam[counter]))
        #     ax.plot(dist_space, kde(dist_space),linewidth=1,color=color)
        #     ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
        #     ax.spines['top'].set_color([0.5, 0.5, 0.5])
        #     ax.spines['right'].set_color([0.5, 0.5, 0.5])
        #     ax.spines['left'].set_color([0.5, 0.5, 0.5])
        #     #ax.text(1.8, 1, np.round(colorparam*100)/100, fontweight='bold')
        #
        # else:
        #     ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
        #     ax.spines['top'].set_color([0.3, 0.3, 0.3])
        #     ax.spines['right'].set_color([0.3, 0.3, 0.3])
        #     ax.spines['left'].set_color([0.3, 0.3, 0.3])
        #     if yy < 14:
        #         ax.xaxis.set_ticks([])
        # if yy < 14:
        #     ax.xaxis.set_ticklabels([])
        # ax.yaxis.set_ticklabels([])
        # ax.yaxis.set_ticks([])
        counter = counter+1
plt.show()
# s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
# s_map.set_array(colorparam)
# fig.subplots_adjust(right=0.92)
# cbar_ax = fig.add_axes([0.94, 0.10, 0.02, 0.7])
# cbar = fig.colorbar(s_map, cax=cbar_ax)
# cbar.set_label('Mean Hs (m)')


asdf





import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.stats.kde import gaussian_kde
dist_space = np.linspace(0, 5, 50)
fig = plt.figure(figsize=(10,10))
colorparam = np.zeros((numClusters*numClusters,))
counter = 0
for xx in range(numClusters):
    for yy in range(numClusters):
        ax = plt.subplot2grid((numClusters,numClusters), (yy,xx), rowspan=1, colspan=1)
        #normalize = mcolors.Normalize(vmin=np.min(colorparams), vmax=np.max(colorparams))
        normalize = mcolors.Normalize(vmin=.5, vmax=1.5)

        ax.set_xlim([0,3])
        ax.set_ylim([0,2])
        data = wHs[xx][yy]
        if len(data)>0:
            kde = gaussian_kde(data)
            colorparam[counter] = np.nanmean(data)
            colormap = cm.Reds
            color = colormap(normalize(colorparam[counter]))
            ax.plot(dist_space, kde(dist_space),linewidth=1,color=color)
            ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
            ax.spines['top'].set_color([0.5, 0.5, 0.5])
            ax.spines['right'].set_color([0.5, 0.5, 0.5])
            ax.spines['left'].set_color([0.5, 0.5, 0.5])
            #ax.text(1.8, 1, np.round(colorparam*100)/100, fontweight='bold')

        else:
            ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
            ax.spines['top'].set_color([0.3, 0.3, 0.3])
            ax.spines['right'].set_color([0.3, 0.3, 0.3])
            ax.spines['left'].set_color([0.3, 0.3, 0.3])
            if yy < 14:
                ax.xaxis.set_ticks([])
        if yy < 14:
            ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
        counter = counter+1
plt.show()
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map.set_array(colorparam)
fig.subplots_adjust(right=0.92)
cbar_ax = fig.add_axes([0.94, 0.10, 0.02, 0.7])
cbar = fig.colorbar(s_map, cax=cbar_ax)
cbar.set_label('Mean Hs (m)')
#cbar = ax.cax.colorbar(s_map, ticks=colorparam, format='%2.2g') # format='%2i' for integer





from scipy.stats.kde import gaussian_kde
dist_space = np.linspace(0, 18, 20)
fig = plt.figure(figsize=(10,10))

for xx in range(numClusters):
    for yy in range(numClusters):
        ax = plt.subplot2grid((15,15), (yy,xx), rowspan=1, colspan=1)
        ax.set_xlim([0,18])
        ax.set_ylim([0,0.3])

        data = wTp[xx][yy]
        if len(data)>0:
            kde = gaussian_kde(data)
            ax.plot(dist_space, kde(dist_space),linewidth=1)
            ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
            ax.spines['top'].set_color([0.5, 0.5, 0.5])
            ax.spines['right'].set_color([0.5, 0.5, 0.5])
            ax.spines['left'].set_color([0.5, 0.5, 0.5])
        else:
            ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
            ax.spines['top'].set_color([0.3, 0.3, 0.3])
            ax.spines['right'].set_color([0.3, 0.3, 0.3])
            ax.spines['left'].set_color([0.3, 0.3, 0.3])
            if yy < 14:
                ax.xaxis.set_ticks([])
        if yy < 14:
            ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
plt.show()



dist_space = np.linspace(-180, 180, 30)
fig = plt.figure(figsize=(10,10))

for xx in range(numClusters):
    for yy in range(numClusters):
        ax = plt.subplot2grid((15,15), (yy,xx), rowspan=1, colspan=1)
        ax.set_xlim([-90,90])
        normalize = mcolors.Normalize(vmin=-15, vmax=15)

        data = wDm[xx][yy]
        if len(data)>0:
            baddata = np.isnan((data))
            data = data[~baddata]
            if len(data)>0:
                kde = gaussian_kde(data)
                colorparam = np.nanmean(data)
                colormap = cm.bwr
                color = colormap(normalize(colorparam))
                ax.plot(dist_space, kde(dist_space),linewidth=1,color=color)
                ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
                ax.spines['top'].set_color([0.5, 0.5, 0.5])
                ax.spines['right'].set_color([0.5, 0.5, 0.5])
                ax.spines['left'].set_color([0.5, 0.5, 0.5])
            else:
                ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
                ax.spines['top'].set_color([0.3, 0.3, 0.3])
                ax.spines['right'].set_color([0.3, 0.3, 0.3])
                ax.spines['left'].set_color([0.3, 0.3, 0.3])
        else:
            ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
            ax.spines['top'].set_color([0.3, 0.3, 0.3])
            ax.spines['right'].set_color([0.3, 0.3, 0.3])
            ax.spines['left'].set_color([0.3, 0.3, 0.3])
            if yy < 14:
                ax.xaxis.set_ticks([])
        if yy < 14:
            ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
plt.show()
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map.set_array(colorparam)
fig.subplots_adjust(right=0.92)
cbar_ax = fig.add_axes([0.94, 0.10, 0.02, 0.7])
cbar = fig.colorbar(s_map, cax=cbar_ax)
cbar.set_label('Dir (deg)')



dist_space = np.linspace(-400, 400, 40)
fig = plt.figure(figsize=(10,10))

for xx in range(numClusters):
    for yy in range(numClusters):
        ax = plt.subplot2grid((15,15), (yy,xx), rowspan=1, colspan=1)
        ax.set_xlim([-400,400])
        ax.set_ylim([0, 0.008])
        normalize = mcolors.Normalize(vmin=0, vmax=200)

        data = wLWP[xx][yy]
        if len(data)>0:
            baddata = np.isnan((data))
            data = data[~baddata]
            if len(data)>0:
                kde = gaussian_kde(data)
                colorparam = np.abs(np.nanmean(data))
                colormap = cm.Reds
                color = colormap(normalize(colorparam))
                ax.plot(dist_space, kde(dist_space),linewidth=1,color=color)
                ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
                ax.spines['top'].set_color([0.5, 0.5, 0.5])
                ax.spines['right'].set_color([0.5, 0.5, 0.5])
                ax.spines['left'].set_color([0.5, 0.5, 0.5])
            else:
                ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
                ax.spines['top'].set_color([0.3, 0.3, 0.3])
                ax.spines['right'].set_color([0.3, 0.3, 0.3])
                ax.spines['left'].set_color([0.3, 0.3, 0.3])
        else:
            ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
            ax.spines['top'].set_color([0.3, 0.3, 0.3])
            ax.spines['right'].set_color([0.3, 0.3, 0.3])
            ax.spines['left'].set_color([0.3, 0.3, 0.3])
            if yy < 14:
                ax.xaxis.set_ticks([])
        if yy < 14:
            ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
plt.show()
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map.set_array(colorparam)
fig.subplots_adjust(right=0.92)
cbar_ax = fig.add_axes([0.94, 0.10, 0.02, 0.7])
cbar = fig.colorbar(s_map, cax=cbar_ax)
cbar.set_label('LWP (W$^2$)')




dist_space = np.linspace(0, 100, 40)
fig = plt.figure(figsize=(10,10))

for xx in range(numClusters):
    for yy in range(numClusters):
        ax = plt.subplot2grid((15,15), (yy,xx), rowspan=1, colspan=1)
        ax.set_xlim([-10, 100])
        #ax.set_ylim([0, 0.008])
        normalize = mcolors.Normalize(vmin=0, vmax=50)

        data = wWE[xx][yy]
        if len(data)>0:
            baddata = np.isnan((data))
            data = data[~baddata]
            if len(data)>0:
                kde = gaussian_kde(data)
                colorparam = np.abs(np.nanmean(data))
                colormap = cm.Reds
                color = colormap(normalize(colorparam))
                ax.plot(dist_space, kde(dist_space),linewidth=1,color=color)
                ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
                ax.spines['top'].set_color([0.5, 0.5, 0.5])
                ax.spines['right'].set_color([0.5, 0.5, 0.5])
                ax.spines['left'].set_color([0.5, 0.5, 0.5])
            else:
                ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
                ax.spines['top'].set_color([0.3, 0.3, 0.3])
                ax.spines['right'].set_color([0.3, 0.3, 0.3])
                ax.spines['left'].set_color([0.3, 0.3, 0.3])
        else:
            ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
            ax.spines['top'].set_color([0.3, 0.3, 0.3])
            ax.spines['right'].set_color([0.3, 0.3, 0.3])
            ax.spines['left'].set_color([0.3, 0.3, 0.3])
            if yy < 14:
                ax.xaxis.set_ticks([])
        if yy < 14:
            ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
plt.show()
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map.set_array(colorparam)
fig.subplots_adjust(right=0.92)
cbar_ax = fig.add_axes([0.94, 0.10, 0.02, 0.7])
cbar = fig.colorbar(s_map, cax=cbar_ax)
cbar.set_label('WE (Hs$^2$Tp)')


import matplotlib.dates as mdates
# specify a date to use for the times
zero =time[0]
timeDiff = time[1:]-time[0:-1]
timeNew = [zero + t for t in timeDiff]
# convert datetimes to numbers
zero = mdates.date2num(zero)
timeT = [t-zero for t in mdates.date2num(timeNew)]

# f = plt.figure()
# ax = f.add_subplot(1,1,1)
# ax.plot(timeT)
# ax.plot([0, 600],[7,7],'--')
# ax.plot([0, 600],[14,14],'--')
# ax.plot([0, 600],[21,21],'--')

#ax.yaxis_date()
#ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# add 10% margin on top (since ax.margins seems to not work here)
# ylim = ax.get_ylim()
# ax.set_ylim(None, ylim[1]+0.1*np.diff(ylim))

#plt.show()

# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()
# plt.close()

###

# lets find all of the first bins...
initialStates = np.where((sorted_bmus==0))
tempState = initialStates
jumps = tempState[0][1:]-tempState[0][0:-1]
bigJumps = np.where((jumps > 1))
allJumps = np.append(0,bigJumps[0]+1)

firstSequence = sorted_bmus[initialStates[0][allJumps[0]]:initialStates[0][allJumps[1]]]

matrix = flatarray2

#
# import plotly.graph_objects as go
#
# fig = go.Figure(go.Sankey(
#     arrangement = "snap",
#     node = {
#         "label": ["A", "B", "C", "D", "E", "F"],
#         "x": [0.2, 0.1, 0.5, 0.7, 0.3, 0.5],
#         "y": [0.7, 0.5, 0.2, 0.4, 0.2, 0.3],
#         'pad':10},  # 10 Pixels
#     link = {
#         "source": [0, 0, 1, 2, 5, 4, 3, 5],
#         "target": [5, 3, 4, 3, 0, 2, 2, 3],
#         "value": [1, 5, 1, 1, 1, 1, 1, 2]}))
#
# fig.write_image("fig2.png")

# import plotly.graph_objects as go
#
# fig = go.Figure(go.Sankey(
#     arrangement = "snap",
#     node = {
#         "label": ["A", "B", "C", "D", "E", "F"],
#         "x": [0.2, 0.1, 0.5, 0.7, 0.3, 0.5],
#         "y": [0.7, 0.5, 0.2, 0.4, 0.2, 0.3],
#         'pad':10},  # 10 Pixels
#     link = {
#         "source": [0, 0, 1, 2, 5, 4, 3, 5],
#         "target": [5, 3, 4, 3, 0, 2, 2, 3],
#         "value": [1, 5, 1, 1, 1, 1, 1, 2]}))
#
# fig.write_image("fig2.png")

### Chord option #2

# import holoviews as hv
# from holoviews import opts, dim
# import pandas as pd
# from bokeh.plotting import show
#
# source_list = sorted_bmus[0:-1]
# target_list = sorted_bmus[1:]
#
#
# # create a df from the data
# df_links = pd.DataFrame([source_list, target_list], ["source", "target"]).T
#
# # calculate the number of interactions between entities using pandas groupby
# # for now let's assume 1 - 1 flows
# df_links = df_links.groupby(["source", "target"]).apply(len)
#
# # convert the groupby into a dataframe
# df_links = df_links.to_frame().reset_index()
#
# # rename the 0 column with value
# df_links.rename(columns = {0:"value"}, inplace = True)
#
# # this is our data
# df_links.head()
# hv.extension('boken')
# mychordplt = hv.Chord(df_links)
# show(hv.render(mychordplt))
#

# ## Chord option #1
#
# # chord diagram
# import matplotlib.pyplot as plt
# from matplotlib.path import Path
# import matplotlib.patches as patches
#
# import numpy as np
#
# LW = 0.3
#
# def polar2xy(r, theta):
#     return np.array([r*np.cos(theta), r*np.sin(theta)])
#
# def hex2rgb(c):
#     return tuple(int(c[i:i+2], 16)/256.0 for i in (1, 3 ,5))
#
# def IdeogramArc(start=0, end=60, radius=1.0, width=0.2, ax=None, color=(1,0,0)):
#     # start, end should be in [0, 360)
#     if start > end:
#         start, end = end, start
#     start *= np.pi/180.
#     end *= np.pi/180.
#     # optimal distance to the control points
#     # https://stackoverflow.com/questions/1734745/how-to-create-circle-with-b%C3%A9zier-curves
#     opt = 4./3. * np.tan((end-start)/ 4.) * radius
#     inner = radius*(1-width)
#     verts = [
#         polar2xy(radius, start),
#         polar2xy(radius, start) + polar2xy(opt, start+0.5*np.pi),
#         polar2xy(radius, end) + polar2xy(opt, end-0.5*np.pi),
#         polar2xy(radius, end),
#         polar2xy(inner, end),
#         polar2xy(inner, end) + polar2xy(opt*(1-width), end-0.5*np.pi),
#         polar2xy(inner, start) + polar2xy(opt*(1-width), start+0.5*np.pi),
#         polar2xy(inner, start),
#         polar2xy(radius, start),
#         ]
#
#     codes = [Path.MOVETO,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.LINETO,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CLOSEPOLY,
#              ]
#
#     if ax == None:
#         return verts, codes
#     else:
#         path = Path(verts, codes)
#         print(color)
#         patch = patches.PathPatch(path, facecolor=color, edgecolor=color, lw=LW, alpha=0.5)
#         #patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.4,), lw=LW)
#
#         ax.add_patch(patch)
#
#
# def ChordArc(start1=0, end1=60, start2=180, end2=240, radius=1.0, chordwidth=0.7, ax=None, color=(1,0,0)):
#     # start, end should be in [0, 360)
#     if start1 > end1:
#         start1, end1 = end1, start1
#     if start2 > end2:
#         start2, end2 = end2, start2
#     start1 *= np.pi/180.
#     end1 *= np.pi/180.
#     start2 *= np.pi/180.
#     end2 *= np.pi/180.
#     opt1 = 4./3. * np.tan((end1-start1)/ 4.) * radius
#     opt2 = 4./3. * np.tan((end2-start2)/ 4.) * radius
#     rchord = radius * (1-chordwidth)
#     verts = [
#         polar2xy(radius, start1),
#         polar2xy(radius, start1) + polar2xy(opt1, start1+0.5*np.pi),
#         polar2xy(radius, end1) + polar2xy(opt1, end1-0.5*np.pi),
#         polar2xy(radius, end1),
#         polar2xy(rchord, end1),
#         polar2xy(rchord, start2),
#         polar2xy(radius, start2),
#         polar2xy(radius, start2) + polar2xy(opt2, start2+0.5*np.pi),
#         polar2xy(radius, end2) + polar2xy(opt2, end2-0.5*np.pi),
#         polar2xy(radius, end2),
#         polar2xy(rchord, end2),
#         polar2xy(rchord, start1),
#         polar2xy(radius, start1),
#         ]
#
#     codes = [Path.MOVETO,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              ]
#
#     if ax == None:
#         return verts, codes
#     else:
#         path = Path(verts, codes)
#         #patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.4,), lw=LW)
#         patch = patches.PathPatch(path, facecolor=color, edgecolor=color, lw=LW, alpha=0.5)
#
#         ax.add_patch(patch)
#
# def selfChordArc(start=0, end=60, radius=1.0, chordwidth=0.7, ax=None, color=(1,0,0)):
#     # start, end should be in [0, 360)
#     if start > end:
#         start, end = end, start
#     start *= np.pi/180.
#     end *= np.pi/180.
#     opt = 4./3. * np.tan((end-start)/ 4.) * radius
#     rchord = radius * (1-chordwidth)
#     verts = [
#         polar2xy(radius, start),
#         polar2xy(radius, start) + polar2xy(opt, start+0.5*np.pi),
#         polar2xy(radius, end) + polar2xy(opt, end-0.5*np.pi),
#         polar2xy(radius, end),
#         polar2xy(rchord, end),
#         polar2xy(rchord, start),
#         polar2xy(radius, start),
#         ]
#
#     codes = [Path.MOVETO,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              Path.CURVE4,
#              ]
#
#     if ax == None:
#         return verts, codes
#     else:
#         path = Path(verts, codes)
#         #patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.4,), lw=LW)
#         patch = patches.PathPatch(path, facecolor=color, edgecolor=color, lw=LW, alpha=0.5)
#
#         ax.add_patch(patch)
#
# def chordDiagram(X, ax, colors=None, width=0.1, pad=2, chordwidth=0.7):
#     """Plot a chord diagram
#     Parameters
#     ----------
#     X :
#         flux data, X[i, j] is the flux from i to j
#     ax :
#         matplotlib `axes` to show the plot
#     colors : optional
#         user defined colors in rgb format. Use function hex2rgb() to convert hex color to rgb color. Default: d3.js category10
#     width : optional
#         width/thickness of the ideogram arc
#     pad : optional
#         gap pad between two neighboring ideogram arcs, unit: degree, default: 2 degree
#     chordwidth : optional
#         position of the control points for the chords, controlling the shape of the chords
#     """
#     # X[i, j]:  i -> j
#     x = X.sum(axis = 1) # sum over rows
#     ax.set_xlim(-1.1, 1.1)
#     ax.set_ylim(-1.1, 1.1)
#
#     if colors is None:
#     # use d3.js category10 https://github.com/d3/d3-3.x-api-reference/blob/master/Ordinal-Scales.md#category10
#         colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
#                   '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
#         if len(x) > 10:
#             print('x is too large! Use x smaller than 10')
#         colors = [hex2rgb(colors[i]) for i in range(len(x))]
#
#     # find position for each start and end
#     y = x/np.sum(x).astype(float) * (360 - pad*len(x))
#
#     pos = {}
#     arc = []
#     nodePos = []
#     start = 0
#     for i in range(len(x)):
#         end = start + y[i]
#         arc.append((start, end))
#         angle = 0.5*(start+end)
#         #print(start, end, angle)
#         if -30 <= angle <= 210:
#             angle -= 90
#         else:
#             angle -= 270
#         nodePos.append(tuple(polar2xy(1.1, 0.5*(start+end)*np.pi/180.)) + (angle,))
#         z = (X[i, :]/x[i].astype(float)) * (end - start)
#         ids = np.argsort(z)
#         z0 = start
#         for j in ids:
#             pos[(i, j)] = (z0, z0+z[j])
#             z0 += z[j]
#         start = end + pad
#
#     for i in range(len(x)):
#         start, end = arc[i]
#         print(colors[i])
#         IdeogramArc(start=start, end=end, radius=1.0, ax=ax, color=colors[i], width=width)
#         start, end = pos[(i,i)]
#         selfChordArc(start, end, radius=1.-width, color=colors[i], chordwidth=chordwidth*0.7, ax=ax)
#         for j in range(i):
#             color = colors[i]
#             if X[i, j] > X[j, i]:
#                 color = colors[j]
#             start1, end1 = pos[(i,j)]
#             start2, end2 = pos[(j,i)]
#             ChordArc(start1, end1, start2, end2,
#                      radius=1.-width, color=colors[i], chordwidth=chordwidth, ax=ax)
#
#     #print(nodePos)
#     return nodePos
#
# ##################################
# #if __name__ == "__main__":
# fig = plt.figure(figsize=(6,6))
#     # flux = np.array([[11975,  5871, 8916, 2868],
#     #   [ 1951, 10048, 2060, 6171],
#     #   [ 8010, 16145, 8090, 8045],
#     #   [ 1013,   990,  940, 6907]
#     # ])
# flux = matrix
# ax = plt.axes([0,0,1,1])
#
#     #nodePos = chordDiagram(flux, ax, colors=[hex2rgb(x) for x in ['#666666', '#66ff66', '#ff6666', '#6666ff']])
# #nodePos = chordDiagram(flux, ax)
# tempColorsInd = np.where((colors==1))
# tempColorsInd2 = np.where((colors==0))
#
# tempColors = colors
# tempColors[tempColorsInd] = 0.9999
# tempColors[tempColorsInd2] = 0.0001
# tempColors = tempColors
# nodePos = chordDiagram(flux,ax,colors=tempColors[:,0:3])
# ax.axis('off')
# prop = dict(fontsize=16*0.8, ha='center', va='center')
#
# plt.show()
#

# nodes = ['1', '2', '3', '4', '5', '6', '7', '8']
# for i in range(8):
#     ax.text(nodePos[i][0], nodePos[i][1], nodes[i], rotation=nodePos[i][2], **prop)

    # plt.savefig("example.png", dpi=600,
    #         transparent=True,
    #         bbox_inches='tight', pad_inches=0.02)



# ## Chord option #4 This one works
#
# import numpy as np
# import plotly.graph_objs as go
# import colorlover as cl
# import pandas as pd
#
# def get_spaced_colors(n, randomized=False):
#     if n > 0:
#         max_value = 255
#         interval = max_value / n
#         hues = np.arange(0, max_value, interval)
#         return cl.to_rgb(["hsl(%d,80%%,40%%)" % i for i in hues])
#     else:
#         return None
#
#
# PI = np.pi
#
#
# def check_square(M):
#     d, n = M.shape
#     if d != n:
#         raise ValueError("Data array must be square.")
#     return n
#
#
# def moduloAB(x, a, b):
#     if a >= b:
#         raise ValueError('Incorrect inverval ends')
#     y = (x - a) % (b - a)
#     return y + b if y < 0 else y + a
#
#
# def test_2PI(x):
#     return 0 <= x < 2 * PI
#
#
# def get_ideogram_ends(ideaogram_len, gap):
#     ideo_ends = []
#     left = 0
#     for k in range(len(ideaogram_len)):
#         right = left + ideaogram_len[k]
#         ideo_ends.append([left, right])
#         left = right + gap
#     return ideo_ends
#
#
# def make_ideogram_arc(R, phi, a=50):
#     # R is the circle radius
#     # Phi is a list of the ends angle coordinates of an arc
#     # a is a parameter that controls the number of points to be evaluated
#     if not test_2PI(phi[0]) or not test_2PI(phi[1]):
#         phi = [moduloAB(t, 0, 2*PI) for t in phi]
#     length = (phi[1] - phi[0]) % 2 * PI
#     nr = 5 if length <= PI/4 else int(a * length / PI)
#     if phi[0] < phi[1]:
#         theta = np.linspace(phi[0], phi[1], nr)
#     else:
#         phi = [moduloAB(t, -PI, PI) for t in phi]
#         theta = np.linspace(phi[0], phi[1], nr)
#     return R * np.exp(1j*theta)
#
#
# def map_data(data_matrix, row_value, ideogram_length):
#     n = data_matrix.shape[0]  # square, so same as 1
#     mapped = np.zeros([n, n])
#     for j in range(n):
#         mapped[:, j] = ideogram_length * data_matrix[:, j] / row_value
#     return mapped
#
#
# def make_ribbon_ends(mapped_data, ideo_ends, idx_sort):
#     n = mapped_data.shape[0]
#     ribbon_boundary = np.zeros((n, n+1))
#     for k in range(n):
#         start = ideo_ends[k][0]
#         ribbon_boundary[k][0] = start
#         for j in range(1, n+1):
#             J = idx_sort[k][j-1]
#             ribbon_boundary[k][j] = start + mapped_data[k][J]
#             start = ribbon_boundary[k][j]
#     return [[(ribbon_boundary[k][j], ribbon_boundary[k][j+1])
#              for j in range(n)] for k in range(n)]
#
#
# def control_pts(angle, radius):
#     if len(angle) != 3:
#         raise ValueError('Angle must have len = 3')
#     b_cplx = np.array([np.exp(1j*angle[k]) for k in range(3)])
#     b_cplx[1] = radius * b_cplx[1]
#     return list(zip(b_cplx.real, b_cplx.imag))
#
#
# def ctrl_rib_chords(l, r, radius):
#     if len(l) != 2 or len(r) != 2:
#         raise ValueError('The arc ends must be elements in a list of len 2')
#     return [control_pts([l[j], (l[j]+r[j])/2, r[j]], radius) for j in range(2)]
#
#
# def make_q_bezier(b):
#     if len(b) != 3:
#         raise ValueError('Contaol polygon must have 3 points')
#     A, B, C = b
#     return 'M ' + str(A[0]) + "," + str(A[1]) + " " + "Q " + \
#            str(B[0]) + ", " + str(B[1]) + " " + \
#            str(C[0]) + ", " + str(C[1])
#
#
# def make_ribbon_arc(theta0, theta1):
#     if test_2PI(theta0) and test_2PI(theta1):
#         if theta0 < theta1:
#             theta0 = moduloAB(theta0, -PI, PI)
#             theta1 = moduloAB(theta1, -PI, PI)
#             if theta0 * theta1 > 0:
#                 raise ValueError('Incorrect angle coordinates for ribbon')
#         nr = int(40 * (theta0 - theta1) / PI)
#         if nr <= 2:
#             nr = 3
#         theta = np.linspace(theta0, theta1, nr)
#         pts = np.exp(1j * theta)
#         string_arc = ''
#         for k in range(len(theta)):
#             string_arc += "L " + str(pts.real[k]) + ", " + str(pts.imag[k])+' '
#         return string_arc
#     else:
#         raise ValueError('The angle coords for arc ribbon must be [0, 2*PI]')
#
#
# def make_layout(title):
#     xaxis = dict(showline=False,
#                  zeroline=False,
#                  showgrid=False,
#                  showticklabels=False,
#                  title='')
#     yaxis = {**xaxis, 'scaleanchor': 'x'}
#     return dict(title=title,
#                 xaxis=xaxis,
#                 yaxis=yaxis,
#                 showlegend=False,
#                 margin=dict(t=25, b=25, l=25, r=25),
#                 hovermode='closest',
#                 shapes=[])
#
#
# def make_ideo_shape(path, line_color, fill_color):
#     return dict(
#         line=go.Line(color=line_color, width=0.45),
#         path=path,
#         type='path',
#         fillcolor=fill_color,
#         layer='below'
#     )
#
#
# def make_ribbon(l, r, line_color, fill_color, radius=0.2):
#     poligon = ctrl_rib_chords(l, r, radius)
#     b, c = poligon
#     return dict(line=go.Line(color=line_color, width=0.5),
#                 path=make_q_bezier(b) + make_ribbon_arc(r[0], r[1]) +
#                 make_q_bezier(c[::-1]) + make_ribbon_arc(l[1], l[0]),
#                 type='path',
#                 fillcolor=fill_color,
#                 layer='below')
#
#
# def make_self_rel(l, line_color, fill_color, radius):
#     b = control_pts([l[0], (l[0]+l[1])/2, l[1]], radius)
#     return dict(
#         line=dict(color=line_color, width=0.5),
#         path=make_q_bezier(b) + make_ribbon_arc(l[1], l[0]),
#         type='path',
#         fillcolor=fill_color,
#         layer='below'
#     )
#
#
# def invPerm(perm):
#     inv = [0] * len(perm)
#     for i, s in enumerate(perm):
#         inv[s] = i
#     return inv
#
#
# def make_filled_chord(M):
#     n = M.shape[0]
#     labels = M.columns
#     M = M.T
#     matrix = M.as_matrix()
#     row_sum = [np.sum(matrix[k, :]) for k in range(n)]
#     gap = 2 * PI * 10e-8
#     ideogram_length = 2*PI*np.asarray(row_sum)/sum(row_sum) - gap*np.ones(n)
#     ideo_colors = [x[:3] + "a" + x[3:-1] + ",.75" + x[-1] for x in
#                    get_spaced_colors(len(labels))]
#     mapped_data = map_data(M.as_matrix(), row_sum, ideogram_length)
#     idx_sort = np.argsort(mapped_data, axis=1)
#     ideo_ends = get_ideogram_ends(ideogram_length, gap)
#     ribbon_ends = make_ribbon_ends(mapped_data, ideo_ends, idx_sort)
#     ribbon_color = [n * [ideo_colors[k]] for k in range(n)]
#     layout = make_layout(' ')
#     ribbon_info = []
#     radii_sribb = [0.2] * n
#     for k in range(n):
#         #print('k is {}'.format(k))
#         sigma = idx_sort[k]
#         sigma_inv = invPerm(sigma)
#         for j in range(k, n):
#             if M.iloc[k, j] == 0 and M.iloc[j, k] == 0:
#                 continue
#             eta = idx_sort[j]
#             eta_inv = invPerm(eta)
#             l = ribbon_ends[k][sigma_inv[j]]
#             if j == k:
#                 layout['shapes'].append(
#                     make_self_rel(l,
#                                   'rgb(175,175,175)',
#                                   ideo_colors[k],
#                                   radius=radii_sribb[k]))
#                 z = 0.9 * np.exp(1j * (l[0] + l[1]) / 2)
#                 #print('k is {}'.format(k))
#
#                 text = "{:d}".format(labels[k]) + " co-occurs with " + "{:d}".format(M.iloc[k, k]) + " of its own appearences"
#                 ribbon_info.append(
#                     go.Scatter(x=[z.real],
#                                y=[z.imag],
#                                mode='markers',
#                                text=text,
#                                hoverinfo="text",
#                                marker=dict(size=0.5,
#                                            color=ideo_colors[k])))
#             else:
#                 r = ribbon_ends[j][eta_inv[k]]
#                 zi = 0.9 * np.exp(1j * (l[0] + l[1]) / 2)
#                 zf = 0.9 * np.exp(1j * (r[0] + r[1]) / 2)
#                 texti = "{:d}".format(labels[k]) + " co-occurs with " + \
#                     "{:f}".format(matrix[k][j]) + " of the " + \
#                     "{:d}".format(labels[j]) + " appearences"
#                 textf = "{:d}".format(labels[j]) + " co-occurs with " + \
#                     "{:f}".format(matrix[j][k]) + " of the " + \
#                     "{:d}".format(labels[k]) + " appearences"
#                 ribbon_info.append(
#                     go.Scatter(x=[zi.real],
#                                y=[zi.imag],
#                                mode='markers',
#                                text=texti,
#                                hoverinfo="text",
#                                marker=dict(size=0.5,
#                                            color=ribbon_color[k][j])))
#                 ribbon_info.append(
#                     go.Scatter(x=[zf.real],
#                                y=[zf.imag],
#                                mode='markers',
#                                text=textf,
#                                hoverinfo="text",
#                                marker=dict(size=0.5,
#                                            color=ribbon_color[j][k])))
#                 r = (r[1], r[0])
#                 if matrix[k][j] > matrix[j][k]:
#                     color_of_highest = ribbon_color[k][j]
#                 else:
#                     color_of_highest = ribbon_color[j][k]
#                 layout['shapes'].append(
#                     make_ribbon(l, r, 'rgb(175, 175, 175)',
#                                 color_of_highest))
#     ideograms = []
#     for k in range(len(ideo_ends)):
#         z = make_ideogram_arc(1.1, ideo_ends[k])
#         zi = make_ideogram_arc(1.0, ideo_ends[k])
#         m = len(z)
#         n = len(zi)
#         ideograms.append(
#             go.Scatter(x=z.real,
#                        y=z.imag,
#                        mode='lines',
#                        line=dict(color=ideo_colors[k],
#                                  shape='spline',
#                                  width=0.25),
#                        text="{:d}".format(labels[k])+'<br>'+'{:d}'.format(row_sum[k]),
#                        hoverinfo='text'))
#         path = 'M '
#         for s in range(m):
#             path += str(z.real[s]) + ', ' + str(z.imag[s]) + ' L '
#         Zi = np.array(zi.tolist()[::-1])
#         for s in range(m):
#             path += str(Zi.real[s]) + ", " + str(Zi.imag[s]) + ' L '
#         path += str(z.real[0]) + ' ,' + str(z.imag[0])
#         layout['shapes'].append(make_ideo_shape(path,
#                                                 'rgb(150,150,150)',
#                                                 ideo_colors[k]))
#     data = ideograms + ribbon_info
#     fig = {
#         "data": data,
#         "layout": layout
#     }
#     return fig
#
#
#
# mat2 = matrix
# mat2[0,0] = 0
# mat2[1,1] = 0
# mat2[2,2] = 0
# mat2[3,3] = 0
# mat2[4,4] = 0
# mat2[5,5] = 0
#
# Matrix = pd.DataFrame(mat2)
#
# figtemp = make_filled_chord(Matrix)
#
# figgo = go.Figure(figtemp)
#
# figgo.write_image("fig3.png")

#
# ## Chord option #3
#
# import plotly.graph_objs as go
#
#
#
# def check_data(data_matrix):
#     L, M=data_matrix.shape
#     if L!=M:
#         raise ValueError('Data array must have (n,n) shape')
#     return L
#
# L=check_data(matrix)
#
# PI=np.pi
#
# def moduloAB(x, a, b): #maps a real number onto the unit circle identified with
#                        #the interval [a,b), b-a=2*PI
#         if a>=b:
#             raise ValueError('Incorrect interval ends')
#         y=(x-a)%(b-a)
#         return y+b if y<0 else y+a
#
# def test_2PI(x):
#     return 0<= x <2*PI
#
# row_sum=[np.sum(matrix[k,:]) for k in range(L)]
#
# #set the gap between two consecutive ideograms
# gap=2*PI*0.005
# ideogram_length=2*PI*np.asarray(row_sum)/sum(row_sum)-gap*np.ones(L)
# def get_ideogram_ends(ideogram_len, gap):
#     ideo_ends=[]
#     left=0
#     for k in range(len(ideogram_len)):
#         right=left+ideogram_len[k]
#         ideo_ends.append([left, right])
#         left=right+gap
#     return ideo_ends
#
# ideo_ends=get_ideogram_ends(ideogram_length, gap)
# ideo_ends
#
# def make_ideogram_arc(R, phi, a=50):
#     # R is the circle radius
#     # phi is the list of ends angle coordinates of an arc
#     # a is a parameter that controls the number of points to be evaluated on an arc
#     if not test_2PI(phi[0]) or not test_2PI(phi[1]):
#         phi=[moduloAB(t, 0, 2*PI) for t in phi]
#     length=(phi[1]-phi[0])% 2*PI
#     nr=5 if length<=PI/4 else int(a*length/PI)
#
#     if phi[0] < phi[1]:
#         theta=np.linspace(phi[0], phi[1], nr)
#     else:
#         phi=[moduloAB(t, -PI, PI) for t in phi]
#         theta=np.linspace(phi[0], phi[1], nr)
#     return R*np.exp(1j*theta)
#
# z=make_ideogram_arc(1.3, [11*PI/6, PI/17])
#
# labels=['1', '2', '3', '4', '5']#, '6','7','8','9','10','11','12','13','14','15']
# ideo_colors=['rgba(244, 109, 67, 0.75)',
#              'rgba(253, 174, 97, 0.75)',
#              'rgba(254, 224, 139, 0.75)',
#              'rgba(217, 239, 139, 0.75)',
#              'rgba(166, 217, 106, 0.75)']#,
#              # 'rgba(253, 174, 97, 0.75)',
#              # 'rgba(254, 224, 139, 0.75)',
#              # 'rgba(217, 239, 139, 0.75)',
#              # 'rgba(166, 217, 106, 0.75)',
#              # 'rgba(217, 239, 139, 0.75)',
#              # 'rgba(166, 217, 106, 0.75)',
#              # 'rgba(253, 174, 97, 0.75)',
#              # 'rgba(254, 224, 139, 0.75)',
#              # 'rgba(217, 239, 139, 0.75)',
#              # 'rgba(166, 217, 106, 0.75)']#brewer colors with alpha set on 0.75
#
# def map_data(data_matrix, row_value, ideogram_length):
#     mapped=np.zeros(data_matrix.shape)
#     for j  in range(L):
#         mapped[:, j]=ideogram_length*data_matrix[:,j]/row_value
#     return mapped
#
# mapped_data=map_data(matrix, row_sum, ideogram_length)
# mapped_data
#
# idx_sort=np.argsort(mapped_data, axis=1)
# idx_sort
#
# def make_ribbon_ends(mapped_data, ideo_ends,  idx_sort):
#     L=mapped_data.shape[0]
#     ribbon_boundary=np.zeros((L,L+1))
#     for k in range(L):
#         start=ideo_ends[k][0]
#         ribbon_boundary[k][0]=start
#         for j in range(1,L+1):
#             J=idx_sort[k][j-1]
#             ribbon_boundary[k][j]=start+mapped_data[k][J]
#             start=ribbon_boundary[k][j]
#     return [[(ribbon_boundary[k][j],ribbon_boundary[k][j+1] ) for j in range(L)] for k in range(L)]
#
# ribbon_ends=make_ribbon_ends(mapped_data, ideo_ends,  idx_sort)
#
# def control_pts(angle, radius):
#     #angle is a  3-list containing angular coordinates of the control points b0, b1, b2
#     #radius is the distance from b1 to the  origin O(0,0)
#
#     # if len(angle)!=3:
#     #     raise InvalidInputError('angle must have len =3')
#     b_cplx=np.array([np.exp(1j*angle[k]) for k in range(3)])
#     b_cplx[1]=radius*b_cplx[1]
#     return zip(b_cplx.real, b_cplx.imag)
#
# def ctrl_rib_chords(l, r, radius):
#     # this function returns a 2-list containing control poligons of the two quadratic Bezier
#     #curves that are opposite sides in a ribbon
#     #l (r) the list of angular variables of the ribbon arc ends defining
#     #the ribbon starting (ending) arc
#     # radius is a common parameter for both control polygons
#     if len(l)!=2 or len(r)!=2:
#         raise ValueError('the arc ends must be elements in a list of len 2')
#     return [control_pts([l[j], (l[j]+r[j])/2, r[j]], radius) for j in range(2)]
#
# ribbon_color=[L*[ideo_colors[k]] for k in range(L)]
# ribbon_color[0][4]=ideo_colors[4]
# ribbon_color[1][2]=ideo_colors[2]
# ribbon_color[2][3]=ideo_colors[3]
# ribbon_color[2][4]=ideo_colors[4]
#
#
# def make_q_bezier(b):  # defines the Plotly SVG path for a quadratic Bezier curve defined by the list if its control points
#     # if len(b) != 3:
#     #      raise valueError('control poligon must have 3 points')
#     A, B, C = b
#     return 'M ' + str(A[0]) + ',' + str(A[1]) + ' ' + 'Q ' + \
#        str(B[0]) + ', ' + str(B[1]) + ' ' + \
#        str(C[0]) + ', ' + str(C[1])
#
# b = [(1, 4), (-0.5, 2.35), (3.745, 1.47)]
#
# make_q_bezier(b)
# #print 'ribbon ends starting from the ideogram[2]\n', ribbon_ends[2]
# def make_ribbon_arc(theta0, theta1):
#
#     if test_2PI(theta0) and test_2PI(theta1):
#         if theta0 < theta1:
#             theta0= moduloAB(theta0, -PI, PI)
#             theta1= moduloAB(theta1, -PI, PI)
#             if theta0*theta1>0:
#                 raise ValueError('incorrect angle coordinates for ribbon')
#
#         nr=int(40*(theta0-theta1)/PI)
#         if nr<=2: nr=3
#         theta=np.linspace(theta0, theta1, nr)
#         pts=np.exp(1j*theta)# points on arc in polar complex form
#
#         string_arc=''
#         for k in range(len(theta)):
#             string_arc+='L '+str(pts.real[k])+', '+str(pts.imag[k])+' '
#         return   string_arc
#     else:
#         raise ValueError('the angle coordinates for an arc side of a ribbon must be in [0, 2*pi]')
#
# make_ribbon_arc(np.pi/3, np.pi/6)
#
#
# #import plotly.plotly as py
# #import plotly.figure_factory as ff
#
# #import chart_studio.plotly as py
# #import plot, iplot
# #import chart_studio.grid_objs as go
# #Replace plotly.presentation_objs with chart_studio.presentation_objs
# #    Replace plotly.widgets with chart_studio.widgets
#
#
#
#
# def make_layout(title, plot_size):
#     axis=dict(showline=False, # hide axis line, grid, ticklabels and  title
#           zeroline=False,
#           showgrid=False,
#           showticklabels=False,
#           title=''
#           )
#
#     return go.Layout(title=title,
#                   xaxis=dict(axis),
#                   yaxis=dict(axis),
#                   showlegend=False,
#                   width=plot_size,
#                   height=plot_size,
#                   margin=dict(t=25, b=25, l=25, r=25),
#                   hovermode='closest',
#                   shapes= []# to this list one appends below the dicts defining the ribbon,
#                            #respectively the ideogram shapes
#                  )
#
# def make_ideo_shape(path, line_color, fill_color):
#     #line_color is the color of the shape boundary
#     #fill_collor is the color assigned to an ideogram
#     return  dict(
#                   line=dict(
#                   color=line_color,
#                   width=0.45
#                  ),
#
#             path=  path,
#             type='path',
#             fillcolor=fill_color,
#             layer='below'
#         )
#
# def make_ribbon(l, r, line_color, fill_color, radius=0.2):
#     #l=[l[0], l[1]], r=[r[0], r[1]]  represent the opposite arcs in the ribbon
#     #line_color is the color of the shape boundary
#     #fill_color is the fill color for the ribbon shape
#     poligon=ctrl_rib_chords(l,r, radius)
#     b,c =poligon
#
#     return  dict(
#                 line=dict(
#                 color=line_color, width=0.5
#             ),
#             path=  make_q_bezier(b)+make_ribbon_arc(r[0], r[1])+
#                    make_q_bezier(c[::-1])+make_ribbon_arc(l[1], l[0]),
#             type='path',
#             fillcolor=fill_color,
#             layer='below'
#         )
#
# def make_self_rel(l, line_color, fill_color, radius):
#     #radius is the radius of Bezier control point b_1
#     b=control_pts([l[0], (l[0]+l[1])/2, l[1]], radius)
#     return  dict(
#                 line=dict(
#                 color=line_color, width=0.5
#             ),
#             path=  make_q_bezier(b)+make_ribbon_arc(l[1], l[0]),
#             type='path',
#             fillcolor=fill_color,
#             layer='below'
#         )
#
# def invPerm(perm):
#     # function that returns the inverse of a permutation, perm
#     inv = [0] * len(perm)
#     for i, s in enumerate(perm):
#         inv[s] = i
#     return inv
#
# layout=make_layout('Chord diagram', 400)
#
# radii_sribb=[0.4, 0.30, 0.35, 0.39, 0.12]#,0.4, 0.30, 0.35, 0.39, 0.12,0.4, 0.30, 0.35, 0.39, 0.12]
#
# ribbon_info=[]
# for k in range(L):
#
#     sigma=idx_sort[k]
#     sigma_inv=invPerm(sigma)
#     for j in range(k, L):
#         if matrix[k][j]==0 and matrix[j][k]==0: continue
#         eta=idx_sort[j]
#         eta_inv=invPerm(eta)
#         l=ribbon_ends[k][sigma_inv[j]]
#
#         if j==k:
#             layout['shapes'].append(make_self_rel(l, 'rgb(175,175,175)',ideo_colors[k], radius=radii_sribb[k]))
#             z=0.9*np.exp(1j*(l[0]+l[1])/2)
#             #the text below will be displayed when hovering the mouse over the ribbon
#             text=labels[k]+' commented on '+ '{:d}'.format(matrix[k][k])+' of '+ 'herself Fb posts',
#             ribbon_info.append(go.Scatter(x=[z.real],
#                                        y=[z.imag],
#                                        mode='markers',
#                                        marker=dict(size=0.5, color=ideo_colors[k]),
#                                        text=text,
#                                        hoverinfo='text'
#                                        )
#                               )
#         else:
#             r=ribbon_ends[j][eta_inv[k]]
#             zi=0.9*np.exp(1j*(l[0]+l[1])/2)
#             zf=0.9*np.exp(1j*(r[0]+r[1])/2)
#             #texti and textf are the strings that will be displayed when hovering the mouse
#             #over the two ribbon ends
#             texti=labels[k]+' commented on '+ '{:d}'.format(matrix[k][j])+' of '+\
#                   labels[j]+ ' Fb posts',
#
#             textf=labels[j]+' commented on '+ '{:d}'.format(matrix[j][k])+' of '+\
#             labels[k]+ ' Fb posts',
#             ribbon_info.append(go.Scatter(x=[zi.real],
#                                        y=[zi.imag],
#                                        mode='markers',
#                                        marker=dict(size=0.5, color=ribbon_color[k][j]),
#                                        text=texti,
#                                        hoverinfo='text'
#                                        )
#                               ),
#             ribbon_info.append(go.Scatter(x=[zf.real],
#                                        y=[zf.imag],
#                                        mode='markers',
#                                        marker=dict(size=0.5, color=ribbon_color[k][j]),
#                                        text=textf,
#                                        hoverinfo='text'
#                                        )
#                               )
#             r=(r[1], r[0])#IMPORTANT!!!  Reverse these arc ends because otherwise you get
#                           # a twisted ribbon
#             #append the ribbon shape
#             layout['shapes'].append(make_ribbon(l, r, 'rgb(175,175,175)' , ribbon_color[k][j]))
#
#
# ideograms=[]
# for k in range(len(ideo_ends)):
#     z= make_ideogram_arc(1.1, ideo_ends[k])
#     zi=make_ideogram_arc(1.0, ideo_ends[k])
#     m=len(z)
#     n=len(zi)
#     ideograms.append(go.Scatter(x=z.real,
#                              y=z.imag,
#                              mode='lines',
#                              line=dict(color=ideo_colors[k], shape='spline', width=0.25),
#                              text=labels[k]+'<br>'+'{:d}'.format(row_sum[k]),
#                              hoverinfo='text'
#                              )
#                      )
#
#
#     path='M '
#     for s in range(m):
#         path+=str(z.real[s])+', '+str(z.imag[s])+' L '
#
#     Zi=np.array(zi.tolist()[::-1])
#
#     for s in range(m):
#         path+=str(Zi.real[s])+', '+str(Zi.imag[s])+' L '
#     path+=str(z.real[0])+' ,'+str(z.imag[0])
#
#     layout['shapes'].append(make_ideo_shape(path,'rgb(150,150,150)' , ideo_colors[k]))
#
# data = go.Data(ideograms+ribbon_info)
# fig = go.Figure(data=data, layout=layout)
#
# import plotly.offline as off
# off.init_notebook_mode()
#
# off.iplot(fig, filename='chord-diagram-Fb')

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



# wavedir = '/media/dylananderson/Elements/WIS_ST63218/'
#
# # Need to sort the files to ensure correct temporal order...
# files = os.listdir(wavedir)
# files.sort()
# files_path = [os.path.join(os.path.abspath(wavedir), x) for x in files]
#
# wis = Dataset(files_path[0])
#
# def getWIS(file):
#     waves = Dataset(file)
#
#     waveHs = waves.variables['waveHs'][:]
#     waveTp = waves.variables['waveTp'][:]
#     waveTm = waves.variables['waveTm'][:]
#     waveTm1 = waves.variables['waveTm1'][:]
#     waveTm2 = waves.variables['waveTm2'][:]
#     waveMeanDirection = waves.variables['waveMeanDirection'][:]
#     waveMeanDirectionSwell = waves.variables['waveMeanDirectionSwell'][:]
#     timeW = waves.variables['time'][:]
#
#
#     output = dict()
#     output['waveHs'] = waveHs
#     output['waveTp'] = waveTp
#     output['waveTm'] = waveTm
#     output['waveTm1'] = waveTm1
#     output['waveTm2'] = waveTm2
#     output['waveMeanDirection'] = waveMeanDirection
#     output['waveMeanDirectionSwell'] = waveMeanDirectionSwell
#     output['t'] = timeW
#
#     return output
#
#
# Hs = []
# Tm = []
# Dm = []
# timeWave = []
# for i in files_path:
#     waves = getWIS(i)
#     Hs = np.append(Hs,waves['waveHs'])
#     Tm = np.append(Tm,waves['waveTm'])
#     Dm = np.append(Dm,waves['waveMeanDirection'])
#     timeWave = np.append(timeWave,waves['t'])
#
#
#
#
#
# # waveNorm = Dm - 72
# # neg = np.where((waveNorm > 180))
# # waveNorm[neg[0]] = waveNorm[neg[0]]-360
# # offpos = np.where((waveNorm>90))
# # offneg = np.where((waveNorm<-90))
# # waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
# # waveNorm[offneg[0]] = waveNorm[offneg[0]]*0
# #
# # LWP_means = 1025*np.square(Hs)*Tm*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
# tWave = [DT.datetime.fromtimestamp(x) for x in timeWave]
#
# #fig = plt.figure(figsize=(10,10))
# #plt.plot(tWave,Hs)
# #plt.plot(time,np.ones(len(time),),'.')
#
# HsSurvey = np.nan * np.ones(len(time))
# HsMaxSurvey = np.nan * np.ones(len(time))
# TmSurvey = np.nan * np.ones(len(time))
# DmSurvey = np.nan * np.ones(len(time))
# WPSurvey = np.nan * np.ones(len(time))
# EFSurvey = np.nan * np.ones(len(time))
# EFLSurvey = np.nan * np.ones(len(time))
#
# timeDiff = np.nan * np.ones(len(time))
# counter = 0
# for tt in range(len(time)-1):
#     #within = [date for date in tWave if time[tt] < date < time[tt+1]]
#     withinIndex = np.where((tWave > np.datetime64(time[tt])) & (tWave < np.datetime64(time[tt+1])))
#     if len(withinIndex[0]) < 2:
#         print('index is not working after time {}'.format(time[tt]))
#         print('which is index = {}'.format(tt))
#     else:
#
#         diff = time[tt+1]-time[tt]
#         timeDiff[counter] = diff.days
#         HsSurvey[counter] = np.nanmean(Hs[withinIndex[0]])
#         HsMaxSurvey[counter] = np.nanmax(Hs[withinIndex[0]])
#         TmSurvey[counter] = np.nanmean(Tm[withinIndex[0]])
#         dirs = Dm[withinIndex[0]]-72
#         neg = np.where((dirs > 180))
#         dirs[neg[0]] = dirs[neg[0]]-360
#         waveNorm = dirs
#         offpos = np.where((waveNorm>90))
#         offneg = np.where((waveNorm<-90))
#         waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
#         waveNorm[offneg[0]] = waveNorm[offneg[0]]*0
#         DmSurvey[counter] = np.nanmean(dirs)
#
#         c0 = 9.81*Tm[withinIndex[0]]/(2 * np.pi)
#         HsSquare = np.square(Hs[withinIndex[0]])
#         EF = 1025*9.81*HsSquare*c0/16
#         EFlongshore = (1025*9.81*HsSquare*c0/16)*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
#         EFSurvey[counter] = np.sum(EF)
#         EFLSurvey[counter] = np.nanmean(np.abs(EFlongshore))
#         WPSurvey[counter] = np.nanmean(HsSquare*Tm[withinIndex[0]])
#
#
#         # fig10 = plt.figure(figsize=(12, 12))
#         # ax100 = plt.subplot2grid((4, 3), (0, 0), colspan=3, rowspan=1)
#         # ax100.plot(tWave[withinIndex[0][0]:withinIndex[0][-1]], Hs[withinIndex[0][0:-1]])
#         # ax100.set_ylabel('Hs (m)')
#         # ax100.set_title('Duration = {} days'.format(diff.days))
#         # ax100.text(tWave[withinIndex[0][0]], HsSurvey[counter], 'Average Hs = {:.2f}'.format(HsSurvey[counter]))
#         # ax100.plot([tWave[withinIndex[0][0]], tWave[withinIndex[0][-1]]],[HsSurvey[counter],HsSurvey[counter]])
#         # #     ax7.text(320, 0, '{}'.format(timeSubset[timestep]))
#         # ax101 = plt.subplot2grid((4, 3), (1, 0), colspan=3, rowspan=1)
#         # ax101.plot(tWave[withinIndex[0][0]:withinIndex[0][-1]], Tm[withinIndex[0][0:-1]])
#         # ax101.text(tWave[withinIndex[0][0]], TmSurvey[counter], 'Average Tm = {:.2f}'.format(TmSurvey[counter]))
#         # ax101.plot([tWave[withinIndex[0][0]], tWave[withinIndex[0][-1]]],[TmSurvey[counter],TmSurvey[counter]])
#         # ax101.set_ylabel('Tm (s)')
#         # ax103 = plt.subplot2grid((4, 3), (2, 0), colspan=3, rowspan=1)
#         # ax103.plot(tWave[withinIndex[0][0]:withinIndex[0][-1]], dirs[0:-1])
#         # ax103.text(tWave[withinIndex[0][0]], DmSurvey[counter], 'Average Dm = {:.2f}'.format(DmSurvey[counter]))
#         # ax103.plot([tWave[withinIndex[0][0]], tWave[withinIndex[0][-1]]],[DmSurvey[counter],DmSurvey[counter]])
#         # ax103.set_ylabel('Dm (deg N)')
#         # ax102 = plt.subplot2grid((4, 3), (3, 0), colspan=3, rowspan=1)
#         # ax102.plot(xinterp,alllines[tt])
#         # ax102.plot(xinterp,alllines[tt+1])
#         #
#         #
#         #
#         #
#         # if tt < 10:
#         #     plt.savefig('Conditions00'+str(tt)+'.png')
#         #
#         # elif tt < 100:
#         #     plt.savefig('Conditions0'+str(tt)+'.png')
#         #
#         # else:
#         #     plt.savefig('Conditions'+str(tt)+'.png')
#         # plt.close()
#         counter = counter+1
#
#
#
# #HsReg = HsSurvey.copy()
# HsReg = HsSurvey[~np.isnan(HsSurvey)]
# HsRegMax = HsMaxSurvey[~np.isnan(HsSurvey)]
# TmReg = TmSurvey[~np.isnan(HsSurvey)]
# DmReg = DmSurvey[~np.isnan(HsSurvey)]
# WPReg =WPSurvey[~np.isnan(HsSurvey)]
# EFReg =EFSurvey[~np.isnan(HsSurvey)]
# EFLReg =EFLSurvey[~np.isnan(HsSurvey)]
#
# timeDiff = timeDiff[~np.isnan(HsSurvey)]
#
# wavetime = time[1:len(HsReg)+1]
#
#
#
#
#
#
# cpc1 = P2R(Rt[1:len(HsReg)+1, 0], phitSubset[1:len(HsReg)+1, 0])
# cpc2 = P2R(Rt[1:len(HsReg)+1, 1], phitSubset[1:len(HsReg)+1, 1])
# cpc3 = P2R(Rt[1:len(HsReg)+1, 2], phitSubset[1:len(HsReg)+1, 2])
#
# precpc1 = P2R(Rt[0:len(HsReg), 0], phitSubset[0:len(HsReg), 0])
# precpc2 = P2R(Rt[0:len(HsReg), 1], phitSubset[0:len(HsReg), 1])
# precpc3 = P2R(Rt[0:len(HsReg), 2], phitSubset[0:len(HsReg), 2])
#
# dcpc1 = P2R(Rt[0:len(HsReg)+1, 0], phitSubset[0:len(HsReg)+1, 0])
# dcpc2 = P2R(Rt[0:len(HsReg)+1, 1], phitSubset[0:len(HsReg)+1, 1])
# dcpc3 = P2R(Rt[0:len(HsReg)+1, 2], phitSubset[0:len(HsReg)+1, 2])
#
# realD1Dt = np.real(dcpc1[1:])-np.real(dcpc1[0:-1])
# imagD1Dt = np.imag(dcpc1[1:])-np.imag(dcpc1[0:-1])
# realD2Dt = np.real(dcpc2[1:])-np.real(dcpc2[0:-1])
# imagD2Dt = np.imag(dcpc2[1:])-np.imag(dcpc2[0:-1])
# realD3Dt = np.real(dcpc3[1:])-np.real(dcpc3[0:-1])
# imagD3Dt = np.imag(dcpc3[1:])-np.imag(dcpc3[0:-1])
# # fig4, ax4 = plt.subplots(2,1)
# # ax4[0].plot(time[1:-1],HsSurvey[0:-2])
# # ax4[1].plot(time, np.real(cpc1),'-')
# # ax4[1].plot(time, np.imag(cpc1),'-')
# #
# # ax4[0].set_xlim([time[t1], time[t2]])
# # ax4[1].set_xlim([time[t1], time[t2]])
#
# PC1 = Rt[0:len(HsReg)+1, 0]*np.sin(phit[0:len(HsReg)+1, 0]) + Rt[0:len(HsReg)+1, 0]*np.cos(phit[0:len(HsReg)+1, 0])
# PC2 = Rt[0:len(HsReg)+1, 1]*np.sin(phit[0:len(HsReg)+1, 1]) + Rt[0:len(HsReg)+1, 1]*np.cos(phit[0:len(HsReg)+1, 1])
# PC3 = Rt[0:len(HsReg)+1, 2]*np.sin(phit[0:len(HsReg)+1, 2]) + Rt[0:len(HsReg)+1, 2]*np.cos(phit[0:len(HsReg)+1, 2])
#
# dPC1 = PC1[1:] - PC1[0:-1]
# dPC2 = PC2[1:] - PC2[0:-1]
# dPC3 = PC3[1:] - PC3[0:-1]
#
# prePC1 = Rt[0:len(HsReg), 0]*np.sin(phit[0:len(HsReg), 0]) + Rt[0:len(HsReg), 0]*np.cos(phit[0:len(HsReg), 0])
# prePC2 = Rt[0:len(HsReg), 1]*np.sin(phit[0:len(HsReg), 1]) + Rt[0:len(HsReg), 1]*np.cos(phit[0:len(HsReg), 1])
# prePC3 = Rt[0:len(HsReg), 2]*np.sin(phit[0:len(HsReg), 2]) + Rt[0:len(HsReg), 2]*np.cos(phit[0:len(HsReg), 2])
#
# t1 = 300
# t2 = -1
#
# fig2, ax2 = plt.subplots(1,4)
# ax2[0].plot(EFReg/timeDiff,wavetime)
# ax2[1].pcolor(xg,tg,eofPred2.T, vmin=-1.25, vmax=1.25)
# ax2[2].plot(Rt[1:len(HsReg)+1,1],wavetime)
# ax2[3].plot(phit[1:len(HsReg)+1,1],wavetime)
#
# ax2[0].set_ylim([wavetime[t1], wavetime[t2]])
# ax2[1].set_ylim([wavetime[t1], wavetime[t2]])
# ax2[2].set_ylim([wavetime[t1], wavetime[t2]])
# ax2[3].set_ylim([wavetime[t1], wavetime[t2]])
#
#
#
# stackRI = np.vstack((realD1Dt, imagD1Dt, np.real(cpc1), np.imag(cpc1))).T
#
# stackRI2 = np.vstack((realD2Dt, imagD2Dt)).T
#
# from sklearn.linear_model import LinearRegression
#
# model1 = LinearRegression()
#
# #model1.fit(HsReg.reshape((-1, 1)),stackRI)
# model1.fit(stackRI, WPReg.reshape((-1, 1)))
#
#
# r_sq = model1.score(stackRI, WPReg.reshape((-1, 1)))
#
# X = np.vstack((realD1Dt, imagD1Dt, np.real(cpc1), np.imag(cpc1))).T
# Y = WPReg
#
# #X2 = np.vstack((np.real(precpc2), np.imag(precpc2), HsReg, TmReg, WPReg)).T
# X2 = np.vstack((EFReg/timeDiff, EFLReg/timeDiff, np.real(prePC3))).T
# Y2 = np.abs(np.real(dPC3))
#
# from statsmodels import api
#
# X2 = api.add_constant(X2)
#
# model2 = api.OLS(Y2, X2).fit() ## sm.OLS(output, input)
# predictions = model2.predict(X2)
# model2.summary()

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




