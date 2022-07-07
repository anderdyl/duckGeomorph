
from getdatatestbed import getDataFRF
import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from glob import glob
import os
from scipy.interpolate import interp1d
import sandBarTool.morphLib as mL
import scipy.io
from matplotlib import gridspec
import pickle
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
import datetime
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.metrics import silhouette_score
import scipy.cluster.hierarchy as shc
from scipy.interpolate import interp1d
import matplotlib.cm as cm
from scipy.stats import gaussian_kde



def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


#dbfile = open('10km_morpho_parameters_and_profiles.pickle', 'rb')
# dbfile = open('20km_16surveys_profiles.pickle', 'rb')
dbfile = open('20km_16surveys_fullprofiles_edited.pickle', 'rb')

profileData = pickle.load(dbfile)
dbfile.close()
#
# rawX = np.empty((1,210))
# rawZ = np.empty((1,210))

rawX = np.empty((1,307))
rawZ = np.empty((1,307))

allProfNumbers = []
allNumbers = []
counter = 1
for xx in profileData:
    print(xx)

    if len(xx) > 5:
        if xx == 'dist':
            print('did''t add the distance vector')
        else:
            rawZ = np.vstack((rawZ,profileData[xx]))
            surveyNum = counter*np.ones((1,len(profileData[xx])))
            allNumbers.append(surveyNum)
            profileNumbers = np.arange(0,1600)
            allProfNumbers.append(profileNumbers)
            counter = counter + 1


for xx in range(len(allNumbers)):
    if xx == 0:
        numberSurvs = allNumbers[xx]
        numberProfs = allProfNumbers[xx]
    else:
        numberSurvs = np.hstack((numberSurvs,allNumbers[xx]))
        numberProfs = np.hstack((numberProfs,allProfNumbers[xx]))


rawZ = rawZ[1:,:]
rawX = profileData['dist']
#
with open(r'case_0_499_profile_result.pickle', "rb") as input_file:
    outputProfiles = pickle.load(input_file)

#preProfiles = outputProfiles['pre']
#postProfiles = outputProfiles['post']
xbdist = outputProfiles['dist']

newX = np.arange(0,1199,1)
#newX = xbdist[8:208]

# plt.figure(figsize=(10,10))
allZs = []
allNums = []
allProfNums = []
for xx in range(len(rawZ)):
    tempXSurvey = rawX
    tempZSurvey = rawZ[xx,:]
    indexN = np.where((tempZSurvey < -40))
    if len(indexN[0]) > 0:
        tempZSurvey[indexN[0]] = 1.5*np.ones((np.shape(indexN[0])))
    tempNumber = numberSurvs[0][xx]
    tempProfNumber = numberProfs[xx]
    f2 = interp1d(tempXSurvey, tempZSurvey, kind='linear')
    #if xx == 0:
    zInterpOG = f2(newX)

    zInterp = moving_average(zInterpOG, 20)
    #else:
    #    zInterp = np.vstack((zInterp,f2(newX)))

    #indices1 = np.arange(10,499,1)
    #indices2 = np.arange(500,1000,1)
    #indices = np.hstack((indices1, indices2))
    #for zz in range(len(indices)):
    #    tempIndex = indices[zz]
    #
    #    if zz == 0:
    #        zInterp = f2(newX)
    #    else:
    #        zInterp = np.vstack((zInterp,f2(newX)))

    #zInterp2 = np.fliplr(zInterp)

    # zAligned = 1.5 * np.ones((np.shape(zInterp)))
    zAligned = 4.5 * np.ones((850,))

    numAligned = np.ones((len(zInterp),))*xx
    # lets go for the dune crest...
    indexI = np.where((zInterp == np.max(zInterp)))
    indexI = np.where((zInterp > 3.5) & (zInterp < 4.25))
    if len(indexI[0]) > 0:
        if np.max(zInterp) > 2.5:

            invertZ = (zInterp[(indexI[0][-1] - 70):])


            # zAligned[0:len(invertZ)] = invertZ
            zAligned[0:849] = invertZ[0:849]
            #plt.plot(newX[0:len(invertZ)], invertZ)#, 'k-')
            # tempX = np.arange(0,len(invertZ),1)
            tempX = np.arange(0, 849, 1)

            extendedX = np.arange(0,849,1)


            allZs.append(zAligned)
            allNums.append(tempNumber)
            allProfNums.append(tempProfNumber)


    # to get the front of the dune
    # if len(indexI[0]) > 0:
    #         # invertZ = tempZs[(indexI[0][-1]-100):]
    #     temp = np.diff(indexI)
    #     highInside = np.where((temp > 1))
    #     if len(highInside[0]) > 0:
    #         #invertZ = np.flipud(zInterp[(indexI[0][highInside[0][0]]-100):])
    #         #### invertZ = (zInterp[(indexI[0][highInside[-1][-1]]-150):])
    #         invertZ = (zInterp[(indexI[0][-1] - 40):])
    #
    #     else:
    #         invertZ = zInterp[(indexI[0][-1] - 40):]
    #
    #     # zAligned[0:len(invertZ)] = invertZ
    #     zAligned[0:849] = invertZ[0:849]
    #     #plt.plot(newX[0:len(invertZ)], invertZ)#, 'k-')
    #     # tempX = np.arange(0,len(invertZ),1)
    #     tempX = np.arange(0, 849, 1)
    #
    #     extendedX = np.arange(0,849,1)
    #
    #     allZs.append(zAligned)
    #     allNums.append(tempNumber)


        # if invertZ[0] > -5:
        #     if invertZ[-1] < -6.25:
        #             #     if invertZ[-20] < -7.8:
        #             #         if invertZ[670] < -5.9:
        #             #             if invertZ[625] < -5.82:
        #             #                 if invertZ[716] < -7:
        #             #                     if invertZ[645] < -6.18:
        #             #                         if invertZ[575] < -5.2:
        #             #                             if invertZ[490] < -4.3:
        #                     # if counter == 0:
        #                     #     zAligned = invertZ
        #                     #     plt.plot(newX,zAligned)
        #                     # else:
        #                     #     zAligned = np.vstack((zAligned,invertZ))
        #                     #
        #                     # plt.plot(newX[indexI[0]],zInterp2[zz,indexI[0]],'ro')
        #         extendedLowerBathyInds = np.where((extendedX > 600))
        #         lowerBathyInds = np.where((tempX > 600))
        #         upperBathyInds = np.where((tempX <= 600))
        #         if invertZ[lowerBathyInds[0][-1]] > -6.25:
        #            # xbot = np.array([tempX[lowerBathyInds[0][0]],tempX[lowerBathyInds[0][-1]]])
        #            xbot = np.array([tempX[lowerBathyInds[0][0]],extendedX[extendedLowerBathyInds[0][-1]]])
        #            zbot = np.array([invertZ[lowerBathyInds[0][0]],-8.5])
        #         else:
        #            xbot = np.array([tempX[lowerBathyInds[0][0]],tempX[lowerBathyInds[0][-1]]])
        #            zbot = np.array([invertZ[lowerBathyInds[0][0]],invertZ[lowerBathyInds[0][-1]]])
        #         mb = np.polyfit(xbot, zbot, 1)
        #         f = np.poly1d(mb)
        #         newZ = f(tempX[lowerBathyInds[0]])
        #         newZ2 = f(extendedX[extendedLowerBathyInds[0]])
        #         newInvert = np.hstack((invertZ[upperBathyInds[0]],newZ))
        #         newInvert2 = np.hstack((invertZ[upperBathyInds[0]], newZ2))
        #
        #             #moveAvg = moving_average(newInvert,10)
        #
        #
        #                                                 #moveMid = moving_average(invertZ,20)
        #                                                 #ending = moving_average(invertZ,40)
        #
        #                                                 #together = np.hstack((moveAvg[0:200],moveMid[190:350]))
        #                                                 #together = np.hstack((together,ending[330:]))
        #                                                 #together = moving_average(together, 10)
        #
        #         zAligned[0:len(newInvert2)] = newInvert2
        #         numAligned = xx
        #         plt.plot(extendedX,newInvert2)
        #         counter = counter + 1
        #
        #         # del indexI
        #
        #         tempDiff = np.diff(zAligned)
        #         tempDiffOrig = np.diff(tempZSurvey)
        #         if np.max(tempDiffOrig > 2):
        #             print('big jump in original profile')
        #         else:
        #             if np.max(tempDiff > 2):
        #                 # skip
        #                 print('big jump somewhere in the interpolated profile')
        #             else:
        #                 allZs.append(zAligned)
        #                 allNums.append(tempNumber)


# x = np.ma.MaskedArray.filled(xbdist[8:178])
# x = x-x[0]
# xEnd = np.ma.MaskedArray.filled(xOrig[100:])
# xInitial = np.arange(0,80,4)
# xIM = np.arange(80,100,2)
# xMiddle = np.arange(100,250,1)
# xME = np.arange(250,300,2)
# xME2 = np.arange(300,400,4)
# xME3 = np.arange(400,600,8)
# xME4 = np.arange(600,840,12)
xInitial = np.arange(0,70,7)
xInit2 = np.arange(70,90,2)
xIM = np.arange(90,100,1)
xIM2 = np.arange(100,180,0.5)
xIM3 = np.arange(180,230,1)
xMiddle = np.arange(230,260,2)
xMiddle2 = np.arange(260,440,10)
xME = np.arange(440,840,20)
# xME2 = np.arange(240,400,4)
# xME3 = np.arange(400,600,8)
# xME4 = np.arange(600,840,12)
# #
x = np.hstack((xInitial,xInit2))
x = np.hstack((x,xIM))
x = np.hstack((x,xIM2))
x = np.hstack((x,xIM3))
x = np.hstack((x,xMiddle))
x = np.hstack((x,xMiddle2))
x = np.hstack((x,xME))
# x = np.hstack((x,xME2))
# x = np.hstack((x,xME3))
# x = np.hstack((x,xME4))

# xbeachX = np.ma.filled(rawX[9:189]-rawX[9],0)



for xx in range(len(allZs)):
    if xx == 0:
        z = allZs[xx]
        f2 = interp1d(newX[0:len(z)], z, kind='linear')
        # if xx == 0:
        zDownscaled = f2(x)
        numberS = allNums[xx]
        numberP = allProfNums[xx]
    else:
        z = allZs[xx]
        f2 = interp1d(newX[0:len(z)], z, kind='linear')
        z2 = f2(x)
        zDownscaled = np.vstack((zDownscaled,z2))
        numberS = np.hstack((numberS,allNums[xx]))
        numberP = np.hstack((numberP,allProfNums[xx]))





# zShort = np.fliplr(z[:,0:751])
zShort = zDownscaled #(zDownscaled[:,0:751])
xShort = x#xbeachX#newX[0:751]


# xOrig = np.flipud(np.abs(rawX[0][0][14:]-858))
# xOrig = xOrig-xOrig[0]
# x = np.ma.MaskedArray.filled(xOrig)
# x = xShort






for zz in range(len(zShort)):
    f = interp1d(xShort, zShort[zz,:], kind='linear')
    if zz == 0:
        zInt = f(x)
    else:
        zInt = np.vstack((zInt, f(x)))

# zInt = zShort
# x = xShort





plt.figure(figsize=(10,10))
for xx in range(len(zInt)):
    plt.plot(x,zInt[xx,:])

meanZ = np.mean(zInt,axis=0)
plt.plot(x,meanZ,'k-',linewidth=3)
plt.xlabel('Cross-shore (m)')
plt.ylabel('Elevation (m)')
plt.title('All Profiles (black = mean)')
plt.savefig('allProfilesAligned.png')
plt.close()






plt.figure(figsize=(10,10))
# for xx in range(len(zInt)):
#     plt.plot(x,zInt[xx,:])

meanZ = np.mean(zInt,axis=0)
stdZ = np.std(zInt, axis=0)
sortedZ = np.sort(zInt,axis=0)
max5 = np.mean(sortedZ[-1200:-1,:],axis=0)
min5 = np.mean(sortedZ[0:1200,:],axis=0)
max500 = np.mean(sortedZ[-500:-1,:],axis=0)
min500 = np.mean(sortedZ[0:500,:],axis=0)
plt.plot(x,meanZ,'k-',linewidth=3)
plt.fill_between(x, max500,min500, facecolor='red', alpha=0.5)
# plt.fill_between(x, max5,min5, facecolor='blue', alpha=0.5)

# plt.fill_between(x, meanZ+2*stdZ, meanZ-2*stdZ, facecolor='red', alpha=0.5)
plt.fill_between(x, meanZ+stdZ, meanZ-stdZ, facecolor='blue', alpha=0.5)

plt.xlabel('Cross-shore (m)')
plt.ylabel('Elevation (m)')
plt.title('All Profiles (black = mean)')
plt.savefig('envelopeOfVariability.png')
plt.close()




# from sklearn.decomposition import PCA

demean = np.zeros((np.shape(zInt)))
normZ = np.zeros((np.shape(zInt)))
plt.figure()
for xx in range(len(demean)):
    demean[xx,:] = (zInt[xx, :] - meanZ)
    normZ[xx,:] = demean[xx,:] / stdZ
    plt.plot(x,demean[xx,:])


ipca = PCA(n_components=9)

# PCs = ipca.fit_transform(normZ)
PCs = ipca.fit_transform(zInt)

EOFs = ipca.components_
variance = ipca.explained_variance_
varianceRatio = ipca.explained_variance_ratio_
cumulativeVar = np.cumsum(varianceRatio)

colors = cm.rainbow(np.linspace(0, 1, 16))

plt.figure()
plt.plot(cumulativeVar[0:20],'o-')
plt.ylabel('Cumulative Variance Explained')
plt.xlabel('Mode #')
plt.savefig('varianceExplained.png')
plt.close()



uniqueNumbers = np.unique(numberS)
pcs1 = np.nan * np.zeros((np.shape(uniqueNumbers)))
pcs2 = np.nan * np.zeros((np.shape(uniqueNumbers)))
pcs3 = np.nan * np.zeros((np.shape(uniqueNumbers)))
pcs4 = np.nan * np.zeros((np.shape(uniqueNumbers)))
pcs5 = np.nan * np.zeros((np.shape(uniqueNumbers)))
pcs6 = np.nan * np.zeros((np.shape(uniqueNumbers)))
pcs7 = np.nan * np.zeros((np.shape(uniqueNumbers)))
pcs8 = np.nan * np.zeros((np.shape(uniqueNumbers)))
pcs9 = np.nan * np.zeros((np.shape(uniqueNumbers)))

for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = np.nanmean(PCs[indexN,0],axis=1)
    pcs1[int(tt)-1] = temp[0]
    temp2 = np.nanmean(PCs[indexN,1],axis=1)
    pcs2[int(tt)-1] = temp2[0]
    temp3 = np.nanmean(PCs[indexN,2],axis=1)
    pcs3[int(tt)-1] = temp3[0]
    temp4 = np.nanmean(PCs[indexN,3],axis=1)
    pcs4[int(tt)-1] = temp4[0]
    temp5 = np.nanmean(PCs[indexN,4],axis=1)
    pcs5[int(tt)-1] = temp5[0]
    temp6 = np.nanmean(PCs[indexN,5],axis=1)
    pcs6[int(tt)-1] = temp6[0]
    temp7 = np.nanmean(PCs[indexN,6],axis=1)
    pcs7[int(tt)-1] = temp7[0]
    temp8 = np.nanmean(PCs[indexN,7],axis=1)
    pcs8[int(tt)-1] = temp8[0]
    temp9 = np.nanmean(PCs[indexN,8],axis=1)
    pcs9[int(tt)-1] = temp9[0]


import datetime as dt
dates = np.asarray([dt.datetime(2010,11,1), dt.datetime(2011,6,1), dt.datetime(2011,11,1),dt.datetime(2012,6,1), dt.datetime(2012,11,1),
         dt.datetime(2013,6,1),dt.datetime(2014,6,1),dt.datetime(2015,6,1),
         dt.datetime(2016,6,1),dt.datetime(2016,10,1),dt.datetime(2017,7,1),dt.datetime(2017,11,1),
         dt.datetime(2018,5,1),dt.datetime(2019,4,1),
         dt.datetime(2019,8,1),dt.datetime(2019,11,1)])




formatter = mdates.DateFormatter("%Y-%m")

plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[0,:])
ax1.set_title('EOF 1')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)

data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,0]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]
x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]
ax2.boxplot(data, positions=pos, widths = 20)
ax2.set_xlim([0, (x2-x1).days ])
ax2.set_xticklabels(['{}-{}'.format(ff.year,ff.month) for ff in dates], rotation=45)
ax2.plot(pos,pcs1,'--',color=[0.2,0.2,0.2])
ax2.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
# plt.savefig('EOF1NagsHead.png')![](EOF6NagsHead.png)
# plt.close()
asdfg

import pickle
eofPickle = 'EOF1forPlotting.pickle'
plottingEOF1 = {}
plottingEOF1['x'] = x
plottingEOF1['EOFs'] = EOFs
plottingEOF1['PCs'] = PCs
plottingEOF1['uniqueNumbers'] = uniqueNumbers
plottingEOF1['numberS'] = numberS
plottingEOF1['numberP'] = numberP
plottingEOF1['dates'] = dates

with open(eofPickle,'wb') as f:
    pickle.dump(plottingEOF1, f)


def datetime2datevec(dtime):
    'Return matlab date vector from datetimes'
    return [dtime.year, dtime.month, dtime.day]

mdateVec = [datetime2datevec(x) for x in dates]

nourPCA = dict()
nourPCA['x'] = x
nourPCA['EOFs'] = EOFs
nourPCA['PCs'] = PCs
nourPCA['pcs1'] = pcs1
nourPCA['pcs2'] = pcs2

nourPCA['uniqueNumbers'] = uniqueNumbers
nourPCA['numberS'] = numberS
nourPCA['numberP'] = numberP

nourPCA['mdateVec'] = mdateVec

import scipy.io
scipy.io.savemat('nourishmentPCA.mat',nourPCA)



asdfg





plt.figure(figsize=(10,10))
axis1 = plt.subplot2grid((4,1),(0,0),rowspan=2,colspan=1)
for xx in range(len(zInt)):
    axis1.plot(x,zInt[xx,:], color=[0.5,0.5,0.5])

axis1.fill_between(x, meanZ+stdZ, meanZ-stdZ, facecolor=[0.75,0.75,0.75], alpha=0.5)

meanZ = np.mean(zInt,axis=0)
axis1.plot(x,meanZ,'k-',linewidth=3)
# axis1.set_xlabel('Cross-shore (m)')
# axis1.set_xlabel('')
axis1.set_ylabel('Elevation (m)')
axis1.text(300,6,'All Profiles (black = mean)')
axis1.set_xlim([0,800])
axis1.text(0.02,0.95,'a',transform=axis1.transAxes,fontsize=16,fontweight='bold',va='top')

axis2 = plt.subplot2grid((4,1),(2,0),rowspan=1,colspan=1)
# ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
axis2.plot(x,EOFs[0,:])
# axis2.set_title('EOF 1')
# axis2.set_title('')
axis2.text(360,0.07,'EOF #1')
axis2.set_xlabel('Cross-shore (m)')
axis2.set_xlim([0,800])
axis2.text(0.02,0.95,'b',transform=axis2.transAxes,fontsize=16,fontweight='bold',va='top')

axis3 = plt.subplot2grid((4,1),(3,0),rowspan=1,colspan=1)

data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,0]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]
x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]
axis3.boxplot(data, positions=pos, widths = 20)
axis3.set_xlim([0, (x2-x1).days ])
axis3.set_xticklabels(['{}/{}'.format(ff.month,ff.year) for ff in dates], rotation=45)
axis3.plot(pos,pcs1,'--',color=[0.2,0.2,0.2])
axis3.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
axis3.set_xlabel('Surveys')
axis3.set_ylabel('PCs')
axis3.text(0.02,0.95,'c',transform=axis3.transAxes,fontsize=16,fontweight='bold',va='top')
plt.tight_layout()



plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[1,:])
ax1.set_title('EOF 2')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.xaxis.set_major_formatter(formatter)
# ax2.plot(dates,pcs1,'k--')
#
data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,1]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]

x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]

ax2.boxplot(data, positions=pos, widths = 20)
ax2.set_xlim([0, (x2-x1).days ])
ax2.set_xticklabels(['{}-{}'.format(ff.year,ff.month) for ff in dates], rotation=45)
ax2.plot(pos,pcs2,'--',color=[0.2,0.2,0.2])
ax2.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF2NagsHead.png')
plt.close()






plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[2,:])
ax1.set_title('EOF 3')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.xaxis.set_major_formatter(formatter)
# ax2.plot(dates,pcs1,'k--')
#
data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,2]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]

x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]

ax2.boxplot(data, positions=pos, widths = 20)
ax2.set_xlim([0, (x2-x1).days ])
ax2.set_xticklabels(['{}-{}'.format(ff.year,ff.month) for ff in dates], rotation=45)
ax2.plot(pos,pcs3,'--',color=[0.2,0.2,0.2])
ax2.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF3NagsHead.png')
plt.close()



plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[3,:])
ax1.set_title('EOF 4')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.xaxis.set_major_formatter(formatter)
# ax2.plot(dates,pcs1,'k--')
#
data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,3]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]

x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]

ax2.boxplot(data, positions=pos, widths = 20)
ax2.set_xlim([0, (x2-x1).days ])
ax2.set_xticklabels(['{}-{}'.format(ff.year,ff.month) for ff in dates], rotation=45)
ax2.plot(pos,pcs4,'--',color=[0.2,0.2,0.2])
ax2.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF4NagsHead.png')
plt.close()






plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[4,:])
ax1.set_title('EOF 5')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.xaxis.set_major_formatter(formatter)
# ax2.plot(dates,pcs1,'k--')
#
data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,4]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]

x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]

ax2.boxplot(data, positions=pos, widths = 20)
ax2.set_xlim([0, (x2-x1).days ])
ax2.set_xticklabels(['{}-{}'.format(ff.year,ff.month) for ff in dates], rotation=45)
ax2.plot(pos,pcs5,'--',color=[0.2,0.2,0.2])
ax2.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF5NagsHead.png')
plt.close()






plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[5,:])
ax1.set_title('EOF 6')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.xaxis.set_major_formatter(formatter)
# ax2.plot(dates,pcs1,'k--')
#
data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,5]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]

x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]

ax2.boxplot(data, positions=pos, widths = 20)
ax2.set_xlim([0, (x2-x1).days ])
ax2.set_xticklabels(['{}-{}'.format(ff.year,ff.month) for ff in dates], rotation=45)
ax2.plot(pos,pcs6,'--',color=[0.2,0.2,0.2])
ax2.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF6NagsHead.png')
plt.close()



plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[6,:])
ax1.set_title('EOF 7')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.xaxis.set_major_formatter(formatter)
# ax2.plot(dates,pcs1,'k--')
#
data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,6]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]

x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]

ax2.boxplot(data, positions=pos, widths = 20)
ax2.set_xlim([0, (x2-x1).days ])
ax2.set_xticklabels(['{}-{}'.format(ff.year,ff.month) for ff in dates], rotation=45)
ax2.plot(pos,pcs7,'--',color=[0.2,0.2,0.2])
ax2.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF7NagsHead.png')
plt.close()



plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[7,:])
ax1.set_title('EOF 8')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.xaxis.set_major_formatter(formatter)
# ax2.plot(dates,pcs1,'k--')
#
data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,7]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]

x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]

ax2.boxplot(data, positions=pos, widths = 20)
ax2.set_xlim([0, (x2-x1).days ])
ax2.set_xticklabels(['{}-{}'.format(ff.year,ff.month) for ff in dates], rotation=45)
ax2.plot(pos,pcs8,'--',color=[0.2,0.2,0.2])
ax2.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF8NagsHead.png')
plt.close()





plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[8,:])
ax1.set_title('EOF 9')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.xaxis.set_major_formatter(formatter)
# ax2.plot(dates,pcs1,'k--')
#
data = []
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = PCs[indexN,8]
    data.append(temp[0])
    # temp = np.mean(PCs[indexN,0],axis=1)
    # pcs1[int(tt)] = temp[0]

x1 = dt.datetime(2010,1,1)
x2 = dt.datetime(2020,1,1)
pos = [(yy-x1).days for yy in dates]

ax2.boxplot(data, positions=pos, widths = 20)
ax2.set_xlim([0, (x2-x1).days ])
ax2.set_xticklabels(['{}-{}'.format(ff.year,ff.month) for ff in dates], rotation=45)
ax2.plot(pos,pcs9,'--',color=[0.2,0.2,0.2])
ax2.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF9NagsHead.png')
plt.close()






X_projects = ipca.inverse_transform(PCs)
zOrig = zInt

num = -2100
prof1 = meanZ + EOFs[0,:]*PCs[num,0]
prof2 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1]
prof3 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2]
prof4 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3]
prof5 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]
prof6 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]
prof7 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]
prof8 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]+ EOFs[7,:]*PCs[num,7]
prof9 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]+ EOFs[7,:]*PCs[num,7]+ EOFs[8,:]*PCs[num,8]

plt.figure()
plt.plot(x,meanZ,label='mean profile')
plt.plot(x,zOrig[num,:],'k-',linewidth=2,label='original profile')
# plt.plot(x,X_projects[num,:]+meanZ,'r-', label='Full EOF reconstruction')
plt.plot(x,prof1,label='1 EOF reconstruction')
plt.plot(x,prof2,label='2 EOFs reconstruction')
plt.plot(x,prof3,label='3 EOFs reconstruction')
plt.plot(x,prof4,label='4 EOFs reconstruction')
plt.plot(x,prof5,label='5 EOFs reconstruction')
plt.plot(x,prof6,label='6 EOFs reconstruction')
plt.plot(x,prof7,label='7 EOFs reconstruction')
plt.plot(x,prof8,label='8 EOFs reconstruction')
plt.plot(x,prof9,'m-',label='9 EOFs reconstruction')

plt.plot(x,X_projects[num,:],'k--',label='Inverse')
plt.legend()
plt.title('EOF Reconstruction Example.png')
plt.show()




dataPred = scipy.io.loadmat('repProfiles.mat')
repProfiles = dataPred['repProfiles']

pro1 = meanZ + EOFs[0,:]*repProfiles[0,0] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro2 = meanZ + EOFs[0,:]*repProfiles[0,1] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro3 = meanZ + EOFs[0,:]*repProfiles[0,2] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro4 = meanZ + EOFs[0,:]*repProfiles[0,3] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro5 = meanZ + EOFs[0,:]*repProfiles[0,4] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro6 = meanZ + EOFs[0,:]*repProfiles[0,5] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro7 = meanZ + EOFs[0,:]*repProfiles[0,6] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro8 = meanZ + EOFs[0,:]*repProfiles[0,7] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]

pro9 = meanZ + EOFs[0,:]*repProfiles[0,8] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro10 = meanZ + EOFs[0,:]*repProfiles[0,9] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro11 = meanZ + EOFs[0,:]*repProfiles[0,10] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro12 = meanZ + EOFs[0,:]*repProfiles[0,11] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro13 = meanZ + EOFs[0,:]*repProfiles[0,12] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro14 = meanZ + EOFs[0,:]*repProfiles[0,13] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
pro15 = meanZ + EOFs[0,:]*repProfiles[0,14] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]

# pro1 = meanZ + EOFs[0,:]*repProfiles[0,0] + EOFs[1,:]*repProfiles[1,0] + EOFs[2,:]*repProfiles[2,0] + EOFs[3,:]*repProfiles[3,0] + EOFs[4,:]*repProfiles[4,0]+ EOFs[5,:]*repProfiles[5,0]+ EOFs[6,:]*repProfiles[6,0]+ EOFs[7,:]*repProfiles[7,0]+ EOFs[8,:]*repProfiles[8,0]
# pro2 = meanZ + EOFs[0,:]*repProfiles[0,1] + EOFs[1,:]*repProfiles[1,1] + EOFs[2,:]*repProfiles[2,1] + EOFs[3,:]*repProfiles[3,1] + EOFs[4,:]*repProfiles[4,1]+ EOFs[5,:]*repProfiles[5,1]+ EOFs[6,:]*repProfiles[6,1]+ EOFs[7,:]*repProfiles[7,1]+ EOFs[8,:]*repProfiles[8,1]
# pro3 = meanZ + EOFs[0,:]*repProfiles[0,2] + EOFs[1,:]*repProfiles[1,2] + EOFs[2,:]*repProfiles[2,2] + EOFs[3,:]*repProfiles[3,2] + EOFs[4,:]*repProfiles[4,2]+ EOFs[5,:]*repProfiles[5,2]+ EOFs[6,:]*repProfiles[6,2]+ EOFs[7,:]*repProfiles[7,2]+ EOFs[8,:]*repProfiles[8,2]
# pro4 = meanZ + EOFs[0,:]*repProfiles[0,3] + EOFs[1,:]*repProfiles[1,3] + EOFs[2,:]*repProfiles[2,3] + EOFs[3,:]*repProfiles[3,3] + EOFs[4,:]*repProfiles[4,3]+ EOFs[5,:]*repProfiles[5,3]+ EOFs[6,:]*repProfiles[6,3]+ EOFs[7,:]*repProfiles[7,3]+ EOFs[8,:]*repProfiles[8,3]
# pro5 = meanZ + EOFs[0,:]*repProfiles[0,4] + EOFs[1,:]*repProfiles[1,4] + EOFs[2,:]*repProfiles[2,4] + EOFs[3,:]*repProfiles[3,4] + EOFs[4,:]*repProfiles[4,4]+ EOFs[5,:]*repProfiles[5,4]+ EOFs[6,:]*repProfiles[6,4]+ EOFs[7,:]*repProfiles[7,4]+ EOFs[8,:]*repProfiles[8,4]
# pro6 = meanZ + EOFs[0,:]*repProfiles[0,5] + EOFs[1,:]*repProfiles[1,5] + EOFs[2,:]*repProfiles[2,5] + EOFs[3,:]*repProfiles[3,5] + EOFs[4,:]*repProfiles[4,5]+ EOFs[5,:]*repProfiles[5,5]+ EOFs[6,:]*repProfiles[6,5]+ EOFs[7,:]*repProfiles[7,5]+ EOFs[8,:]*repProfiles[8,5]
# pro7 = meanZ + EOFs[0,:]*repProfiles[0,6] + EOFs[1,:]*repProfiles[1,6] + EOFs[2,:]*repProfiles[2,6] + EOFs[3,:]*repProfiles[3,6] + EOFs[4,:]*repProfiles[4,6]+ EOFs[5,:]*repProfiles[5,6]+ EOFs[6,:]*repProfiles[6,6]+ EOFs[7,:]*repProfiles[7,6]+ EOFs[8,:]*repProfiles[8,6]
# pro8 = meanZ + EOFs[0,:]*repProfiles[0,7] + EOFs[1,:]*repProfiles[1,7] + EOFs[2,:]*repProfiles[2,7] + EOFs[3,:]*repProfiles[3,7] + EOFs[4,:]*repProfiles[4,7]+ EOFs[5,:]*repProfiles[5,7]+ EOFs[6,:]*repProfiles[6,7]+ EOFs[7,:]*repProfiles[7,7]+ EOFs[8,:]*repProfiles[8,7]

import matplotlib.cm as cm
etcolors = cm.viridis(np.linspace(0, 1, 15))

plt.figure()
plt.plot(x,pro1,'k-',linewidth=2,label=dates[0])
plt.plot(x,pro2,label=dates[1],color=etcolors[0])
plt.plot(x,pro3,label=dates[2],color=etcolors[1])
plt.plot(x,pro4,label=dates[3],color=etcolors[2])
plt.plot(x,pro5,label=dates[4],color=etcolors[3])
plt.plot(x,pro6,label=dates[5],color=etcolors[4])
plt.plot(x,pro7,label=dates[6],color=etcolors[5])
plt.plot(x,pro8,label=dates[7],color=etcolors[6])
plt.plot(x,pro9,label=dates[8],color=etcolors[7])
plt.plot(x,pro10,label=dates[9],color=etcolors[8])
plt.plot(x,pro11,label=dates[10],color=etcolors[9])
plt.plot(x,pro12,label=dates[11],color=etcolors[10])
plt.plot(x,pro13,label=dates[12],color=etcolors[11])
plt.plot(x,pro14,label=dates[13],color=etcolors[12])
plt.plot(x,pro15,label=dates[14],color=etcolors[13])

plt.legend()


# num = 5100
# prof7 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]
# prof7l = meanZ + EOFs[0,:]*PCs[num-15,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]
# prof7r = meanZ + EOFs[0,:]*PCs[num+15,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]
#
# # prof8 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]+ EOFs[7,:]*PCs[num,7]
#
# plt.figure()
# #plt.plot(x,meanZ,label='mean profile')
# #plt.plot(x,zOrig[num,:],'k-',linewidth=2,label='original profile')
# # plt.plot(x,X_projects[num,:]+meanZ,'r-', label='Full EOF reconstruction')
# plt.plot(x,prof7,'k-',label='center')
# plt.plot(x,prof7l,label='left')
# plt.plot(x,prof7r,label='right')
#
# # plt.plot(x,prof8,'m-',label='8 EOFs reconstruction')
# #plt.plot(x,X_projects[num,:],'k--',label='Inverse')
# plt.legend()
# plt.title('Profile Neighbors 150m away just changing EOF1')
# plt.show()
#








PCs1 = PCs[:,0]
PCs2 = PCs[:,1]
PCs3 = PCs[:,2]
PCs4 = PCs[:,3]
PCs5 = PCs[:,4]
PCs6 = PCs[:,5]
PCs7 = PCs[:,6]
PCs8 = PCs[:,7]
PCs9 = PCs[:,8]

#
#
# PCsall = []
# PCsallArray = np.empty((1,9))
# for yyy in range(17):
#     indexS = np.where((numberS == yyy))
#     PCs1temp = PCs1[indexS]
#     PCs1center = PCs1temp[16:-16]
#     PCs1left = PCs1temp[0:-32]
#     PCs1right = PCs1temp[32:]
#     PCs1diff= np.divide(np.subtract(PCs1left,PCs1right),2)
#
#     # PCs1diff= np.diff((PCs1left,PCs1right))
#
#     PCs2temp = PCs2[indexS]
#     PCs2center = PCs2temp[16:-16]
#     PCs2left = PCs2temp[0:-2]
#     PCs2right = PCs2temp[2:]
#
#     PCs3temp = PCs3[indexS]
#     PCs3center = PCs3temp[16:-16]
#     PCs4temp = PCs4[indexS]
#     PCs4center = PCs4temp[16:-16]
#     PCs5temp = PCs5[indexS]
#     PCs5center = PCs5temp[16:-16]
#     PCs6temp = PCs6[indexS]
#     PCs6center = PCs6temp[16:-16]
#     PCs7temp = PCs7[indexS]
#     PCs7center = PCs7temp[16:-16]
#     PCs8temp = PCs8[indexS]
#     PCs8center = PCs8temp[16:-16]
#
#     PCsTemp = np.empty((len(PCs1center),9))
#     PCsTemp[:,0] = PCs1center
#     PCsTemp[:,1] = PCs2center
#     PCsTemp[:,2] = PCs3center
#     PCsTemp[:,3] = PCs4center
#     PCsTemp[:,4] = PCs5center
#     PCsTemp[:,5] = PCs6center
#     PCsTemp[:,6] = PCs7center
#     PCsTemp[:,7] = PCs8center
#
#     PCsTemp[:,8] = PCs1diff
#     # PCsTemp[:,8] = PCs1right
#
#     if yyy > 0:
#         PCsall.append(PCsTemp)
#         PCsallArray = np.vstack((PCsallArray,PCsTemp))
#
#
# PCsallArray = PCsallArray[1:,:]
#
#
#
# PCs1c = PCsallArray[:,0]
# PCs2c = PCsallArray[:,1]
# PCs3c = PCsallArray[:,2]
# PCs4c = PCsallArray[:,3]
# PCs5c = PCsallArray[:,4]
# PCs6c = PCsallArray[:,5]
# PCs7c = PCsallArray[:,6]
# PCs8c = PCsallArray[:,7]
#
# PCs1diff = PCsallArray[:,8]
# #PCs1r = PCsallArray[:,8]

#
# plt.figure(figsize=(14,8))
# ax31 = plt.subplot2grid((5,5),(0,0),rowspan=1,colspan=1)
# ax31.hist(PCs1)
# ax31.set_ylabel('EOF1')
# ax31.set_xlabel('EOF1')
# xy = np.vstack([PCs2,PCs1])
# z = gaussian_kde(xy)(xy)
# ax32 = plt.subplot2grid((5,5),(0,1),rowspan=1,colspan=1)
# ax32.scatter(PCs2,PCs1,5,c=z)
#
# ax33 = plt.subplot2grid((5,5),(0,2),rowspan=1,colspan=1)
# xy = np.vstack([PCs3,PCs1])
# z = gaussian_kde(xy)(xy)
# ax33.scatter(PCs3,PCs1,5,c=z)
# ax34 = plt.subplot2grid((5,5),(0,3),rowspan=1,colspan=1)
# xy = np.vstack([PCs4,PCs1])
# z = gaussian_kde(xy)(xy)
# ax34.scatter(PCs4,PCs1,5,c=z)
# ax35 = plt.subplot2grid((5,5),(0,4),rowspan=1,colspan=1)
# xy = np.vstack([PCs5,PCs1])
# z = gaussian_kde(xy)(xy)
# ax35.scatter(PCs5,PCs1,5,c=z)
#
# ax36 = plt.subplot2grid((5,5),(1,1),rowspan=1,colspan=1)
# ax36.hist(PCs2)
# ax36.set_ylabel('EOF2')
# ax36.set_xlabel('EOF2')
# ax37 = plt.subplot2grid((5,5),(1,2),rowspan=1,colspan=1)
# xy = np.vstack([PCs3,PCs2])
# z = gaussian_kde(xy)(xy)
# ax37.scatter(PCs3,PCs2,5,c=z)
# ax38 = plt.subplot2grid((5,5),(1,3),rowspan=1,colspan=1)
# xy = np.vstack([PCs4,PCs2])
# z = gaussian_kde(xy)(xy)
# ax38.scatter(PCs4,PCs2,5,c=z)
# ax39 = plt.subplot2grid((5,5),(1,4),rowspan=1,colspan=1)
# xy = np.vstack([PCs5,PCs2])
# z = gaussian_kde(xy)(xy)
# ax39.scatter(PCs5,PCs2,5,c=z)
#
# ax40 = plt.subplot2grid((5,5),(2,2),rowspan=1,colspan=1)
# ax40.hist(PCs3)
# ax40.set_ylabel('EOF3')
# ax40.set_xlabel('EOF3')
# ax41 = plt.subplot2grid((5,5),(2,3),rowspan=1,colspan=1)
# xy = np.vstack([PCs4,PCs3])
# z = gaussian_kde(xy)(xy)
# ax41.scatter(PCs4,PCs3,5,c=z)
# ax42 = plt.subplot2grid((5,5),(2,4),rowspan=1,colspan=1)
# xy = np.vstack([PCs5,PCs3])
# z = gaussian_kde(xy)(xy)
# ax42.scatter(PCs5,PCs3,5,c=z)
#
# ax43 = plt.subplot2grid((5,5),(3,3),rowspan=1,colspan=1)
# ax43.hist(PCs4)
# ax43.set_ylabel('EOF4')
# ax43.set_xlabel('EOF4')
# ax44 = plt.subplot2grid((5,5),(3,4),rowspan=1,colspan=1)
# xy = np.vstack([PCs5,PCs4])
# z = gaussian_kde(xy)(xy)
# ax44.scatter(PCs5,PCs4,5,c=z)
#
# ax45 = plt.subplot2grid((5,5),(4,4),rowspan=1,colspan=1)
# ax45.hist(PCs5)
# ax45.set_ylabel('EOF5')
# ax45.set_xlabel('EOF5')
#
# plt.tight_layout()
# # plt.savefig('EOFJointProbabilities.png')
# # plt.close()
#
#



from numpy.random import seed
from numpy.random import rand
# seed random number generator
seed(1)
# generate some random numbers
profileNumbers = np.floor(rand(2000000)*len(PCs1))





#
#
#
# dataCop = []
# for xx in range(len(PCs1)):
#     dataCop.append(list([PCs1[xx],PCs2[xx],PCs3[xx],PCs4[xx],PCs5[xx],PCs6[xx]]))#,PCs7[xx],PCs8[xx]]))
#
# #
# # PCs1c = PCsallArray[:,0]
# # PCs2c = PCsallArray[:,1]
# # PCs3c = PCsallArray[:,2]
# # PCs4c = PCsallArray[:,3]
# # PCs5c = PCsallArray[:,4]
# # PCs6c = PCsallArray[:,5]
# # PCs7c = PCsallArray[:,6]
# # PCs8c = PCsallArray[:,7]
# # # PCs1r = PCsallArray[:,8]
#
# # dataCopWlr = []
# # for xx in range(len(PCs1c)):
# #     dataCopWlr.append(list([PCs1c[xx],PCs2c[xx],PCs3c[xx],PCs4c[xx],PCs5c[xx],PCs6c[xx],PCs7c[xx],PCs1l[xx]]))#,PCs1r[xx]])) #,PCs8[xx]]))
#
#
# from copula import pyCopula
# cop = pyCopula.Copula(dataCop)
# samples = cop.gendata(1000)
#
# samplePCs1 = [item[0] for item in samples]
# samplePCs2 = [item[1] for item in samples]
# samplePCs3 = [item[2] for item in samples]
# samplePCs4 = [item[3] for item in samples]
# samplePCs5 = [item[4] for item in samples]
# samplePCs6 = [item[5] for item in samples]
# samplePCs7 = [item[6] for item in samples]
# #samplePCs1l = [item[7] for item in samples]
# # samplePCs1r = [item[8] for item in samples]
#
# samplePCs8 = [item[7] for item in samples]







# for num in range(1000):
#         plt.figure()
#
#         ztemp = meanZ + EOFs[0, :] * samplePCs1[num] + EOFs[1, :] * samplePCs2[num] + EOFs[2, :] * samplePCs3[num] + EOFs[3, :] * \
#                 samplePCs4[num] + EOFs[4, :] * samplePCs5[num] + EOFs[5, :] * samplePCs6[num] + EOFs[6, :] * samplePCs7[num] + EOFs[7,:]*samplePCs8[num]
#
#         plt.plot(x, ztemp, '-')
#         plt.plot((x[0],x[-1]),(0,0),'.-')
#         plt.title('profile{}'.format(num))
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/hypoNourishments/profile{}.png'.format(num))
#         plt.close()










import pickle

with open(r"filteredStorms.pickle", "rb") as input_file:
   filteredStorms = pickle.load(input_file)

filteredHs = filteredStorms['filteredHs']
filteredTp = filteredStorms['filteredTp']
filteredDm = filteredStorms['filteredDm']
filteredNTR = filteredStorms['filteredNTR']
filteredDur = filteredStorms['filteredDur']


dataCop = []
for xx in range(len(filteredHs)):
    dataCop.append(list([filteredHs[xx],filteredTp[xx],filteredDm[xx],filteredNTR[xx],filteredDur[xx]]))


from copula import pyCopula
cop = pyCopula.Copula(dataCop)
samples = cop.gendata(2000000)

sampleHs = [item[0] for item in samples]
sampleTp = [item[1] for item in samples]
sampleDm = [item[2] for item in samples]
sampleNTR = [item[3] for item in samples]
sampleDur = [item[4] for item in samples]


sampleHs = np.array(sampleHs)
sampleTp = np.array(sampleTp)
sampleDm = np.array(sampleDm)
sampleNTR = np.array(sampleNTR)
sampleDur = np.array(sampleDur)

tooSteep = np.where((sampleHs > 3) & (sampleTp < 6))
sampleHs = np.delete(sampleHs,tooSteep)
sampleTp = np.delete(sampleTp,tooSteep)
sampleDm = np.delete(sampleDm,tooSteep)
sampleNTR = np.delete(sampleNTR,tooSteep)
sampleDur = np.delete(sampleDur,tooSteep)

tooSteep = np.where((sampleHs > 3.5) & (sampleTp < 6.6))
sampleHs = np.delete(sampleHs,tooSteep)
sampleTp = np.delete(sampleTp,tooSteep)
sampleDm = np.delete(sampleDm,tooSteep)
sampleNTR = np.delete(sampleNTR,tooSteep)
sampleDur = np.delete(sampleDur,tooSteep)

tooSteep = np.where((sampleHs > 4) & (sampleTp < 7.5))
sampleHs = np.delete(sampleHs,tooSteep)
sampleTp = np.delete(sampleTp,tooSteep)
sampleDm = np.delete(sampleDm,tooSteep)
sampleNTR = np.delete(sampleNTR,tooSteep)
sampleDur = np.delete(sampleDur,tooSteep)

tooSteep = np.where((sampleHs > 4.5) & (sampleTp < 8.1))
sampleHs = np.delete(sampleHs,tooSteep)
sampleTp = np.delete(sampleTp,tooSteep)
sampleDm = np.delete(sampleDm,tooSteep)
sampleNTR = np.delete(sampleNTR,tooSteep)
sampleDur = np.delete(sampleDur,tooSteep)

tooSteep = np.where((sampleHs > 5) & (sampleTp < 9))
sampleHs = np.delete(sampleHs,tooSteep)
sampleTp = np.delete(sampleTp,tooSteep)
sampleDm = np.delete(sampleDm,tooSteep)
sampleNTR = np.delete(sampleNTR,tooSteep)
sampleDur = np.delete(sampleDur,tooSteep)

dataS = np.empty((len(sampleHs),5))
dataS[:,0] = sampleHs
dataS[:,1] = sampleTp
dataS[:,2] = sampleNTR
dataS[:,3] = sampleDur
dataS[:,4] = sampleDm


from numpy.random import seed
from numpy.random import rand
# seed random number generator
seed(1)
# generate some random numbers
dataT = np.empty((len(sampleHs),4))
dataT[:,0] = rand(len(sampleHs))*360
dataT[:,1] = rand(len(sampleHs))*360
dataT[:,2] = rand(len(sampleHs))*360
dataT[:,3] = rand(len(sampleHs))*360
# dataT[:,4] = rand(len(sampleHs))*360
# dataT[:,5] = rand(len(sampleHs))*360
# dataT[:,6] = rand(len(sampleHs))*360
# dataT[:,7] = rand(len(sampleHs))*360
#


M2amp = 0.49
N2amp = 0.114
K1amp = 0.087
S2amp = 0.088
O1amp = 0.059
SAamp = 0.059
SSAamp = 0.037
P1amp = 0.028
K2amp = 0.023
NU2amp = 0.023


M2speed = 28.984104
N2speed = 28.43973
K1speed = 15.041069
S2speed = 30.0
O1speed = 13.943035
SAspeed = 0.0410686
SSAspeed = 0.0821373
P1speed = 14.958931
K2speed = 30.082138
NU2speed = 28.512583

M2period = 360/M2speed
N2period = 360/N2speed
K1period = 360/K1speed
S2period = 360/S2speed
O1period = 360/O1speed
SAperiod = 360/SAspeed
SSAperiod = 360/SSAspeed
P1period = 360/P1speed
K2period = 360/K2speed
NU2period = 360/NU2speed


# example figure
hours = np.arange(0,24*5,1)
randNum = 115
M2tide = M2amp*np.cos((hours*(2*np.pi/M2period)+dataT[randNum,0]*(2*np.pi/360)))
N2tide = N2amp*np.cos((hours*(2*np.pi/N2period)+dataT[randNum,1]*(2*np.pi/360)))
K1tide = K1amp*np.cos((hours*(2*np.pi/K1period)+dataT[randNum,2]*(2*np.pi/360)))
S2tide = S2amp*np.cos((hours*(2*np.pi/S2period)+dataT[randNum,3]*(2*np.pi/360)))
# O1tide = O1amp*np.cos((hours*(2*np.pi/O1period)+dataT[randNum,4]*(2*np.pi/360)))
# P1tide = P1amp*np.cos((hours*(2*np.pi/P1period)+dataT[randNum,5]*(2*np.pi/360)))
# K2tide = K2amp*np.cos((hours*(2*np.pi/K2period)+dataT[randNum,6]*(2*np.pi/360)))
# NU2tide = NU2amp*np.cos((hours*(2*np.pi/NU2period)+dataT[randNum,7]*(2*np.pi/360)))

plt.figure()
ax100 = plt.subplot2grid((2,2),(0,0),rowspan=1,colspan=2)
ax100.plot(hours,M2tide,label='M2')
ax100.plot(hours,N2tide,label='N2')
ax100.plot(hours,K1tide,label='K1')
ax100.plot(hours,S2tide,label='S2')
# ax100.plot(hours,O1tide,label='O1')
# ax100.plot(hours,P1tide,label='P1')
# ax100.plot(hours,K2tide,label='K2')
# ax100.plot(hours,NU2tide,label='NU2')
ax100.set_xlabel('Hourly values (5 days worth)')
ax100.set_title('Random Phase tide construction')
ax100.legend()

ax101 = plt.subplot2grid((2,2),(1,0),rowspan=1,colspan=2)
ax101.plot(hours,M2tide+N2tide+K1tide+S2tide,label='Top 4')
# ax101.plot(hours,M2tide+N2tide+K1tide+S2tide+O1tide+P1tide+K2tide+NU2tide,label='Top 8 Daily')
ax101.legend()

from scipy.optimize import curve_fit

def objective(t,a,b,c):
    return a*t*t + b*t + c

def objective2(t,a,b):
    return a*t + b

def objective(t,a,b,c):
    return a*t*t + b*t + c


lt = np.array([2015,2045])
ncSLR26low = np.array([0,0.12])
ncSLR26mean = np.array([0,0.18])
ncSLR26high = np.array([0,0.238])
ncSLR85low = np.array([0,0.139])
ncSLR85mean = np.array([0,0.205])
ncSLR85high = np.array([0,0.269])


ht1 = np.arange(2000,2045,1/8766)
ht2 = np.arange(2045,2070,2/8766)
ht = np.hstack((ht1,ht2))
time = np.array([2020,2025,2030,2035,2040,2045,2050,2055,2060,2065,2070,2075,2080,2085,2090,2095,2100])
USACE = np.array([0.025,0.057,0.09,0.125,0.161,0.199,0.238,0.278,0.32,0.363,0.407,0.453,0.5,0.548,0.598,0.649,0.701])
USACElow = np.array([0.004,0.027,0.051,0.075,0.099,0.123,0.147,0.171,0.194,0.218,0.242,0.266,0.29,0.314,0.337,0.361,0.385])
linearSLR = (ht-2020)*0.00479 + 0.0
linearSLRlow = (time-2020)*(0.00479-0.00062) + 0.0
linearSLRhigh = (time-2020)*(0.00479+0.00062) + 0.0

popt, _ = curve_fit(objective,time,USACE)
popt2, _ = curve_fit(objective,time,USACElow)
a,b,c = popt
a2,b2,c2 = popt2
x_new = np.arange(2000,2070,1/8766)
usaceLow = objective(ht,a,b,c)
usaceMed = objective(ht,a2,b2,c2)

popt26low, _ = curve_fit(objective2,lt,ncSLR26low)
popt26mean, _ = curve_fit(objective2,lt,ncSLR26mean)
popt26high, _ = curve_fit(objective2,lt,ncSLR26high)
popt85low, _ = curve_fit(objective2,lt,ncSLR85low)
popt85mean, _ = curve_fit(objective2,lt,ncSLR85mean)
popt85high, _ = curve_fit(objective2,lt,ncSLR85high)

y26low = objective2(ht,popt26low[0],popt26low[1])
y26mean = objective2(ht,popt26mean[0],popt26mean[1])
y26high = objective2(ht,popt26high[0],popt26high[1])
y85low = objective2(ht,popt85low[0],popt85low[1])
y85mean = objective2(ht,popt85mean[0],popt85mean[1])
y85high = objective2(ht,popt85high[0],popt85high[1])


plt.figure()
#plt.plot(time,USACE)

plt.plot(time,linearSLRlow,'k--')
plt.plot(time,linearSLRhigh,'k--',label='tide fit uncertainty')
plt.plot(ht,linearSLR,'k-',label='linear tide gauge obs. (Duck)')
plt.plot(ht,y26low,'r-.',label='2.6 uncertainty')
plt.plot(lt,ncSLR26mean,'r-o',label='NC SLR RCP2.6mean')
plt.plot(ht,y26mean,'r--')
plt.plot(ht,y26high,'r-.')
plt.plot(ht,y85low,'m-.')
plt.plot(lt,ncSLR85mean,'m-o',label='NC SLR RCP8.5mean')
plt.plot(ht,y85mean,'m--')
plt.plot(ht,y85high,'m-.',label='8.5 uncertainty')
plt.plot(ht,usaceLow,label='USACE Mid',linewidth=2)
plt.plot(ht,usaceMed,label='USACE Low',linewidth=2)
plt.xlabel('Year')
plt.ylabel('NAVD88 (m)')
plt.legend()
plt.title('Future Sea Level Rise')
plt.ylim([-0.14,0.75])






slrNumbers = np.floor(rand(len(sampleHs))*len(usaceMed))

dataSLR = np.empty((len(sampleHs),1))
for ii in range(len(sampleHs)):
    top = usaceMed[int(slrNumbers[ii])]
    bot = usaceLow[int(slrNumbers[ii])]

    slrTemp = np.random.uniform(bot,top)
    dataSLR[ii] = slrTemp




dataN = np.empty((len(sampleHs),9))
for ii in range(len(sampleHs)):

    # if PCs1[int(profileNumbers[ii])] == np.nan:
    #     profileNumbers[ii] = profileNumbers[ii]+1
    # else:
    dataN[ii,0] = PCs1[int(profileNumbers[ii])]
    dataN[ii,1] = PCs2[int(profileNumbers[ii])]
    dataN[ii,2] = PCs3[int(profileNumbers[ii])]
    dataN[ii,3] = PCs4[int(profileNumbers[ii])]
    dataN[ii,4] = PCs5[int(profileNumbers[ii])]
    dataN[ii,5] = PCs6[int(profileNumbers[ii])]
    dataN[ii,6] = PCs7[int(profileNumbers[ii])]
    dataN[ii,7] = PCs8[int(profileNumbers[ii])]
    dataN[ii,8] = PCs9[int(profileNumbers[ii])]



#dataN[:,7] = samplePCs1l[0:len(sampleHs)]
#dataN[:,8] = samplePCs1r[0:len(sampleHs)]



# Should we normalize dataN by the variance before MDA?


def Normalize(data, ix_scalar, ix_directional, minis=[], maxis=[]):
    '''
    Normalize data subset - norm = val - min) / (max - min)

    data - data to normalize, data variables at columns.
    ix_scalar - scalar columns indexes
    ix_directional - directional columns indexes
    '''

    data_norm = np.zeros(data.shape) * np.nan

    # calculate maxs and mins
    if minis==[] or maxis==[]:

        # scalar data
        for ix in ix_scalar:
            v = data[:, ix]
            mi = np.amin(v)
            ma = np.amax(v)
            data_norm[:, ix] = (v - mi) / (ma - mi)
            minis.append(mi)
            maxis.append(ma)

        minis = np.array(minis)
        maxis = np.array(maxis)

    # # max and mins given
    # else:
    #
    #     # scalar data
    #     for c, ix in enumerate(ix_scalar):
    #         v = data[:, ix]
    #         mi = minis[c]
    #         ma = maxis[c]
    #         data_norm[:,ix] = (v - mi) / (ma - mi)

    # directional data
    for ix in ix_directional:
        v = data[:,ix]
        data_norm[:,ix] = v * np.pi / 180.0


    return data_norm, minis, maxis

def DeNormalize(data_norm, ix_scalar, ix_directional, minis, maxis):
    '''
    DeNormalize data subset for MaxDiss algorithm

    data - data to normalize, data variables at columns.
    ix_scalar - scalar columns indexes
    ix_directional - directional columns indexes
    '''

    data = np.zeros(data_norm.shape) * np.nan

    # scalar data
    for c, ix in enumerate(ix_scalar):
        v = data_norm[:,ix]
        mi = minis[c]
        ma = maxis[c]
        data[:, ix] = v * (ma - mi) + mi

    # directional data
    for ix in ix_directional:
        v = data_norm[:,ix]
        data[:, ix] = v * 180 / np.pi

    return data

def Normalized_Distance(M, D, ix_scalar, ix_directional):
    '''
    Normalized distance

    M -
    D -
    ix_scalar - scalar columns indexes
    ix_directional - directional columns indexes
    '''

    dif = np.zeros(M.shape)

    # scalar
    for ix in ix_scalar:
        dif[:,ix] = D[:,ix] - M[:,ix]

    # directional
    for ix in ix_directional:
        ab = np.absolute(D[:,ix] - M[:,ix])
        dif[:,ix] = np.minimum(ab, 2*np.pi - ab)/np.pi

    dist = np.sum(dif**2,1)
    return dist

def MaxDiss_Simplified_NoThreshold(data, num_centers, ix_scalar, ix_directional):
    '''
    Normalize data and calculate centers using
    maxdiss simplified no-threshold algorithm

    data - data to apply maxdiss algorithm, data variables at columns
    num_centers - number of centers to calculate
    ix_scalar - scalar columns indexes
    ix_directional - directional columns indexes
    '''

    # TODO: REFACTOR / OPTIMIZE

    print('\nMaxDiss waves parameters: {0} --> {1}\n'.format(
        data.shape[0], num_centers))

    # normalize scalar and directional data
    data_norm, minis, maxis = Normalize(data, ix_scalar, ix_directional)

    # mda seed
    seed = np.where(data_norm[:,0] == np.amax(data_norm[:,0]))[0][0]

    # initialize centroids subset
    subset = np.array([data_norm[seed]])
    train = np.delete(data_norm, seed, axis=0)

    # repeat till we have desired num_centers
    n_c = 1
    while n_c < num_centers:
        m = np.ones((train.shape[0],1))
        m2 = subset.shape[0]

        if m2 == 1:
            xx2 = np.repeat(subset, train.shape[0], axis=0)
            d_last = Normalized_Distance(train, xx2, ix_scalar, ix_directional)

        else:
            xx = np.array([subset[-1,:]])
            xx2 = np.repeat(xx, train.shape[0], axis=0)
            d_prev = Normalized_Distance(train, xx2, ix_scalar, ix_directional)
            d_last = np.minimum(d_prev, d_last)

        qerr, bmu = np.amax(d_last), np.argmax(d_last)

        if not np.isnan(qerr):
            subset = np.append(subset, np.array([train[bmu,:]]), axis=0)
            train = np.delete(train, bmu, axis=0)
            d_last = np.delete(d_last, bmu, axis=0)

            # log
            fmt = '0{0}d'.format(len(str(num_centers)))
            print('   MDA centroids: {1:{0}}/{2:{0}}'.format(
                fmt, subset.shape[0], num_centers), end='\r')

        n_c = subset.shape[0]
    print('\n')

    # normalize scalar and directional data
    centroids = DeNormalize(subset, ix_scalar, ix_directional, minis, maxis)

    return centroids



# dataAll = np.empty((len(sampleHs),18))
# dataAll[:,0:5] = dataS
# dataAll[:,5:9] = dataT
# dataAll[:,9:18] = dataN
# dataAll2 = dataAll[:,0:17]
dataNplusSLR = np.hstack((dataN,dataSLR))
dataNplusSLRplusS = np.hstack((dataNplusSLR,dataS))
dataNplusSLRplusSplusT = np.hstack((dataNplusSLRplusS,dataT))

scalars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
directionals = [14,15,16,17,18]
num_centers = 2100


# normalize scalar and directional data
data_norm, minis, maxis = Normalize(data=dataNplusSLRplusSplusT, ix_scalar=scalars, ix_directional=directionals, minis=[], maxis=[])
# now adjust for variance of the EOFs...
varDivide = np.divide(varianceRatio,varianceRatio[0])
onesDivide = np.array((1,1,1,1,1,1,1,1,1,1))
allDivide = np.hstack((varDivide,onesDivide))

data_normDivide = np.multiply(data_norm,allDivide)

# mda seed
seed = np.where(data_normDivide[:, 0] == np.amax(data_normDivide[:, 0]))[0][0]

# initialize centroids subset
subset = np.array([data_normDivide[seed]])
train = np.delete(data_normDivide, seed, axis=0)

# repeat till we have desired num_centers
n_c = 1
while n_c < num_centers:
    m = np.ones((train.shape[0], 1))
    m2 = subset.shape[0]

    if m2 == 1:
        xx2 = np.repeat(subset, train.shape[0], axis=0)
        d_last = Normalized_Distance(train, xx2, ix_scalar=scalars, ix_directional=directionals)

    else:
        xx = np.array([subset[-1, :]])
        xx2 = np.repeat(xx, train.shape[0], axis=0)
        d_prev = Normalized_Distance(train, xx2, ix_scalar=scalars, ix_directional=directionals)
        d_last = np.minimum(d_prev, d_last)

    qerr, bmu = np.amax(d_last), np.argmax(d_last)

    if not np.isnan(qerr):
        subset = np.append(subset, np.array([train[bmu, :]]), axis=0)
        train = np.delete(train, bmu, axis=0)
        d_last = np.delete(d_last, bmu, axis=0)

        # log
        fmt = '0{0}d'.format(len(str(num_centers)))
        print('   MDA centroids: {1:{0}}/{2:{0}}'.format(
            fmt, subset.shape[0], num_centers), end='\r')

    n_c = subset.shape[0]
print('\n')

subsetDivide = np.divide(subset,allDivide)
# normalize scalar and directional data
centroids = DeNormalize(subsetDivide, ix_scalar=scalars, ix_directional=directionals, minis=minis, maxis=maxis)
mdaSub = centroids
# print('\nMaxDiss waves parameters: {0} --> {1}\n'.format(data_norm.shape[0], num_centers))


# mdaSub = MaxDiss_Simplified_NoThreshold(dataNplusSplusT,2100,scalars,directionals)


# Python function to print common elements in three sorted arrays
def findCommon(ar1, ar2, ar3, n1, n2, n3):
    # Initialize starting indexes for ar1[], ar2[] and ar3[]
    i, j, k = 0, 0, 0

    # Iterate through three arrays while all arrays have elements
    while (i < n1 and j < n2 and k < n3):

        # If x = y and y = z, print any of them and move ahead
        # in all arrays
        if (ar1[i] == ar2[j] and ar2[j] == ar3[k]):
            print(ar1[i])
            output = ar1[i]
            i += 1
            j += 1
            k += 1

        # x < y
        elif ar1[i] < ar2[j]:
            i += 1

        # y < z
        elif ar2[j] < ar3[k]:
            j += 1

        # We reach here when x > y and z < y, i.e., z is smallest
        else:
            k += 1
    return output

# mdaSubtest = mda.MaxDiss_Simplified_NoThreshold(dataN,100,[0,1,2,3,4,5,6,7,8],[])
# mdaSubtest2 = mda.MaxDiss_Simplified_NoThreshold(dataNplusS,100,[0,1,2,3,4,5,6,7,8,9,10,11,12],[13])#,minis=np.amin(dataNplusS,axis=0),maxis=np.amax(dataNplusS,axis=0))
# mdaSubtest3 = mda.MaxDiss_Simplified_NoThreshold(dataS,100,[0,1,2,3],[4])
# mdaSubtest4 = mda.MaxDiss_Simplified_NoThreshold(dataSplusT,100,[0,1,2,3],[4,5,6,7,8])
# mdaSubtest5 = mda.MaxDiss_Simplified_NoThreshold(dataSplusTplusN,100,[0,1,2,3,9,10,11,12,13,14,15,16,17],[4,5,6,7,8])
# # mdaSubtest6 = mda.MaxDiss_Simplified_NoThreshold(dataAll,100,[0,1,2,3,9,10,11,12,13,14,15,16,17],[4,5,6,7,8])


plt.figure()
nr = 0
nc = 0
for num in range(100):

        ax = plt.subplot2grid((5,5),(nr,nc),rowspan=1,colspan=1)
        finder = np.where(((mdaSub[num,0]-.001)<PCs1) & ((mdaSub[num,0]+.001)>PCs1))
        finder2 = np.where(((mdaSub[num,1]-.001)<PCs2) & ((mdaSub[num,1]+.001)>PCs2))
        finder3 = np.where(((mdaSub[num,2]-.001)<PCs3) & ((mdaSub[num,2]+.001)>PCs3))

        ztemp = meanZ + EOFs[0, :]*mdaSub[num,0] + EOFs[1, :]*mdaSub[num,1] + EOFs[2, :] * mdaSub[num,2] + EOFs[3, :] * \
                mdaSub[num,3] + EOFs[4, :]*mdaSub[num,4] + EOFs[5, :]*mdaSub[num,5] + EOFs[6, :]*mdaSub[num,6] + EOFs[7,:]*mdaSub[num,7] + EOFs[8,:]*mdaSub[num,8]

        output = findCommon(finder[0],finder2[0],finder3[0],len(finder[0]),len(finder2[0]),len(finder3[0]))
        ax.plot(x, ztemp, '-',label='Xbeach Profile')
        ax.plot((x[0],x[-1]),(0,0),'.-',label='zero water level')
        ax.plot((x[0],x[-1]),(mdaSub[num,9],mdaSub[num,9]),label='MSL')
        # ax.plot(xShort, zInt[output, :],'--')
        # ax.set_title('profile{}'.format(num))
        ax.text(250,2,'profile{}'.format(num))
        ax.text(250,-1,'MSL={}m'.format(np.round(mdaSub[num,9]*1000)/1000))
        ax.text(250,-2,'HS={}m'.format(np.round(mdaSub[num,10]*10)/10))
        ax.text(250,-3,'Surge={}m'.format(np.round(mdaSub[num,12]*10)/10))

        ax.set_xlim([0, 400])
        ax.set_ylim([-4, 4])
        if nc < 4:
            nc = nc+1
        else:
            nc = 0
            nr = nr+1
        # plt.savefig('/home/dylananderson/projects/duckGeomorph/hypoNourishments/profile{}.png'.format(num))
        # plt.close()


plt.figure()
plt.hist(mdaSub[0:500,9],20)
plt.xlabel('MSL (m)')
plt.title('First 500 scenarios')


clusterPickle = 'hypoNourishStormTides2.pickle'
output = {}
output['hypo'] = mdaSub

import pickle
with open(clusterPickle,'wb') as f:
    pickle.dump(output, f)



eofPickle = 'nourishmentProfileEOFs2.pickle'
outputEOFs = {}
outputEOFs['meanZ'] = meanZ
outputEOFs['EOFs'] = EOFs
outputEOFs['PCs'] = PCs
outputEOFs['variance'] = variance
outputEOFs['varianceRatio'] = varianceRatio
outputEOFs['cumulativeVar'] = cumulativeVar
outputEOFs['zInt'] = zInt
outputEOFs['x'] = x

with open(eofPickle,'wb') as f:
    pickle.dump(outputEOFs, f)


hypoPickle = 'hypotheticalTrainingConditions2.pickle'
outputHypo = {}
outputHypo['dataN'] = dataN
outputHypo['dataT'] = dataT
outputHypo['dataS'] = dataS
outputHypo['M2amp'] = M2amp
outputHypo['N2amp'] = N2amp
outputHypo['K1amp'] = K1amp
outputHypo['S2amp'] = S2amp
outputHypo['M2period'] = M2period
outputHypo['N2period'] = N2period
outputHypo['K1period'] = K1period
outputHypo['S2period'] = S2period
outputHypo['M2speed'] = M2speed
outputHypo['N2speed'] = N2speed
outputHypo['K1speed'] = K1speed
outputHypo['S2speed'] = S2speed

outputHypo['PCs1'] = PCs1
outputHypo['PCs2'] = PCs2
outputHypo['PCs3'] = PCs3
outputHypo['PCs4'] = PCs4
outputHypo['PCs5'] = PCs5
outputHypo['PCs6'] = PCs6
outputHypo['PCs7'] = PCs7
outputHypo['PCs8'] = PCs8
outputHypo['PCs9'] = PCs9

with open(hypoPickle,'wb') as f:
    pickle.dump(outputHypo, f)



#
#
# ax21 = plt.subplot2grid((5,5),(1,0),rowspan=1,colspan=1)
# ax21.plot(sampleNTR,sampleHs,'.',color=[0.5,0.5,0.5])
# ax22 = plt.subplot2grid((5,5),(2,0),rowspan=1,colspan=1)
# ax22.plot(sampleNTR,sampleTp,'.',color=[0.5,0.5,0.5])
# ax23 = plt.subplot2grid((5,5),(3,0),rowspan=1,colspan=1)
# ax23.plot(sampleNTR,sampleDm,'.',color=[0.5,0.5,0.5])
# ax24 = plt.subplot2grid((5,5),(4,0),rowspan=1,colspan=1)
# ax24.plot(sampleNTR,sampleDur,'.',color=[0.5,0.5,0.5])
# ax24.set_xlabel('NTR (m)')
# ax25 = plt.subplot2grid((5,5),(2,1),rowspan=1,colspan=1)
# ax25.plot(sampleHs,sampleTp,'.',color=[0.5,0.5,0.5])
# ax26 = plt.subplot2grid((5,5),(3,1),rowspan=1,colspan=1)
# ax26.plot(sampleHs,sampleDm,'.',color=[0.5,0.5,0.5])
# ax27 = plt.subplot2grid((5,5),(4,1),rowspan=1,colspan=1)
# ax27.plot(sampleHs,sampleDur,'.',color=[0.5,0.5,0.5])
# ax27.set_xlabel('Hs (m)')
# ax28 = plt.subplot2grid((5,5),(3,2),rowspan=1,colspan=1)
# ax28.plot(sampleTp,sampleDm,'.',color=[0.5,0.5,0.5])
# ax29 = plt.subplot2grid((5,5),(4,2),rowspan=1,colspan=1)
# ax29.plot(sampleTp,sampleDur,'.',color=[0.5,0.5,0.5])
# ax29.set_xlabel('Tp (s)')
# ax30 = plt.subplot2grid((5,5),(4,3),rowspan=1,colspan=1)
# ax30.plot(sampleDm,sampleDur,'.',color=[0.5,0.5,0.5])
# ax30.set_xlabel('Dm (s)')
#
#
#
# mdaSub = mda.MaxDiss_Simplified_NoThreshold(data,200,[0,1,3,4],[2])
#
# ax21.plot(mdaSub[0:100,3],mdaSub[0:100,0],'.',color='red')
# ax22.plot(mdaSub[0:100,3],mdaSub[0:100,1],'.',color='red')
# ax23.plot(mdaSub[0:100,3],mdaSub[0:100,2],'.',color='red')
# ax24.plot(mdaSub[0:100,3],mdaSub[0:100,4],'.',color='red')
# ax25.plot(mdaSub[0:100,0],mdaSub[0:100,1],'.',color='red')
# ax26.plot(mdaSub[0:100,0],mdaSub[0:100,2],'.',color='red')
# ax27.plot(mdaSub[0:100,0],mdaSub[0:100,4],'.',color='red')
# ax28.plot(mdaSub[0:100,1],mdaSub[0:100,2],'.',color='red')
# ax29.plot(mdaSub[0:100,1],mdaSub[0:100,4],'.',color='red')
# ax30.plot(mdaSub[0:100,2],mdaSub[0:100,4],'.',color='red')
#
# ax21.plot(mdaSub[100:200,3],mdaSub[100:200,0],'.',color='orange')
# ax22.plot(mdaSub[100:200,3],mdaSub[100:200,1],'.',color='orange')
# ax23.plot(mdaSub[100:200,3],mdaSub[100:200,2],'.',color='orange')
# ax24.plot(mdaSub[100:200,3],mdaSub[100:200,4],'.',color='orange')
# ax25.plot(mdaSub[100:200,0],mdaSub[100:200,1],'.',color='orange')
# ax26.plot(mdaSub[100:200,0],mdaSub[100:200,2],'.',color='orange')
# ax27.plot(mdaSub[100:200,0],mdaSub[100:200,4],'.',color='orange')
# ax28.plot(mdaSub[100:200,1],mdaSub[100:200,2],'.',color='orange')
# ax29.plot(mdaSub[100:2000,1],mdaSub[100:200,4],'.',color='orange')
# ax30.plot(mdaSub[100:200,2],mdaSub[100:200,4],'.',color='orange')
#
#
# clusterPickle = 'hypotheticalStorms.pickle'
# output = {}
# output['hypotheticalStorms'] = mdaSub


















#
# plt.figure(figsize=(12,8))
# ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
# ax1.plot(x,EOFs[1,:])
# ax1.set_title('EOF 2')
# ax1.set_xlabel('Cross-shore (m)')
# ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,1],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
# ax2.set_ylabel('PCs')
# plt.tight_layout()
# plt.savefig('EOF2NagsHead.png')












import numpy as np
import scipy.optimize as sciopt
T = 15
d = 10
g = 9.81
kguess_deep = np.square((2*np.pi/T))/g
kguess_shallow = 2*np.pi/(T*np.sqrt(np.multiply(g,d)))

# Use average of these two limiting values
kguess = 0.5*(kguess_shallow + kguess_deep)

def ksolver(x, d,T,g):
    y = np.square((2 * np.pi / T)) - g * x * np.tanh(x * d)
    return y
# Use scipy.optimize.fsolve to find the true value of k
k = sciopt.fsolve(func=ksolver,x0=kguess,args=(d,T,g))
# Print out L
L = 2*np.pi/k
print(L)


















