

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
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.metrics import silhouette_score
import scipy.cluster.hierarchy as shc
from scipy.interpolate import interp1d
import matplotlib.cm as cm
from scipy.stats import gaussian_kde
#
#
# dbfile = open('5km_morpho_parameters.pickle', 'rb')
# data = pickle.load(dbfile)
# dbfile.close()
# # alllines = data['alllines']
# # xinterp = data['xinterp']
# # time = data['time']
# duneCrestElevation = []
# duneCrestHeight = []
# duneToeElevation = []
# beachWidth = []
# shoreline = []
# duneSlope = []
# beachSlope = []
#
# for xx in data:
#     print(xx)
#     duneCrestElevation = np.append(duneCrestElevation, data[xx]['DuneCrest_Elevation'])
#     duneCrestHeight = np.append(duneCrestHeight, data[xx]['DuneCrest_Height'])
#     duneToeElevation = np.append(duneToeElevation, data[xx]['DuneToe_Elevation'])
#     beachWidth = np.append(beachWidth, data[xx]['Beach_Width'])
#     shoreline = np.append(shoreline, data[xx]['Shoreline'])
#     duneSlope = np.append(duneSlope, data[xx]['Dune_Slope'])
#     beachSlope = np.append(beachSlope, data[xx]['Beach_Slope'])
#
#
# goodPicks = np.where((duneToeElevation>0.25))
#
# duneCrestHeight = duneCrestHeight[goodPicks[0]]
# duneCrestElevation = duneCrestElevation[goodPicks[0]]
# duneToeElevation = duneToeElevation[goodPicks[0]]
# beachSlope = beachSlope[goodPicks[0]]
# beachWidth = beachWidth[goodPicks[0]]
# shoreline = shoreline[goodPicks[0]]
# duneSlope = duneSlope[goodPicks[0]]
#
# goodPicks2 = np.where((beachWidth>10))
#
# duneCrestHeight = duneCrestHeight[goodPicks2[0]]
# duneCrestElevation = duneCrestElevation[goodPicks2[0]]
# duneToeElevation = duneToeElevation[goodPicks2[0]]
# beachSlope = beachSlope[goodPicks2[0]]
# beachWidth = beachWidth[goodPicks2[0]]
# shoreline = shoreline[goodPicks2[0]]
# duneSlope = duneSlope[goodPicks2[0]]
#
#
# avgDuneCrestHeight = np.mean(duneCrestHeight)
# avgDuneCrestElevation = np.mean(duneCrestElevation)
# avgDuneToeElevation = np.mean(duneToeElevation)
# avgDuneSlope = np.mean(duneSlope)
# avgBeachWidth = np.mean(beachWidth)
#
#
# x = np.arange(0,150,1)
# duneY = -avgDuneSlope*x + avgDuneCrestElevation
# duneToeX = (avgDuneToeElevation-avgDuneCrestHeight)/-avgDuneSlope
# shorelineX = duneToeX + avgBeachWidth
# shorelineXPlusSTD = shorelineX+np.std(beachWidth)
# shorelineXMinusSTD = shorelineX-np.std(beachWidth)
# shorelineXPlus2STD = shorelineX+2*np.std(beachWidth)
# shorelineXMinus2STD = shorelineX-2*np.std(beachWidth)
#
# plt.figure()
# ax1 = plt.subplot2grid((1,1),(0,0),colspan=1,rowspan=1)
# # ax1.plot(x[0:20],duneY[0:20],label='Avg Dune Slope')
# ax1.plot(0,avgDuneCrestElevation,'ro',label='Avg Dune Crest Elevation')
# ax1.plot(duneToeX,avgDuneToeElevation,'bo',label='Avg Dune Toe Elevation')
# ax1.plot([0, duneToeX],[avgDuneCrestElevation,avgDuneToeElevation],'k--')
# ax1.plot(shorelineX,0,'ko',label='Avg Beach Width')
# ax1.plot(shorelineXPlusSTD,0,'go',label='SD of Beach Width')
# ax1.plot(shorelineXMinusSTD,0,'go')
# ax1.plot(shorelineXPlus2STD,0,'yo',label='Two SD of Beach Width')
# ax1.plot(shorelineXMinus2STD,0,'yo')
# ax1.plot([duneToeX,shorelineX],[avgDuneToeElevation,0],'k--')
# ax1.plot([duneToeX,shorelineXMinus2STD],[avgDuneToeElevation,0],'--',color=[0.5, 0.5, 0.5])
# ax1.plot([duneToeX,shorelineXMinusSTD],[avgDuneToeElevation,0],'--',color=[0.5, 0.5, 0.5])
# ax1.plot([duneToeX,shorelineXPlus2STD],[avgDuneToeElevation,0],'--',color=[0.5, 0.5, 0.5])
# ax1.plot([duneToeX,shorelineXPlusSTD],[avgDuneToeElevation,0],'--',color=[0.5, 0.5, 0.5])
#
# ax1.legend()
#
#
# plt.figure()
# ax31 = plt.subplot2grid((5,5),(0,0),rowspan=1,colspan=1)
# ax31.hist(duneCrestHeight)
# ax31.set_ylabel('Dune Crest Height')
# ax31.set_xlabel('Dune Crest Height')
# ax32 = plt.subplot2grid((5,5),(0,1),rowspan=1,colspan=1)
# ax32.scatter(duneCrestElevation,duneCrestHeight,10,color='black')
# ax33 = plt.subplot2grid((5,5),(0,2),rowspan=1,colspan=1)
# ax33.scatter(duneToeElevation,duneCrestHeight,10,color='black')
# ax34 = plt.subplot2grid((5,5),(0,3),rowspan=1,colspan=1)
# ax34.scatter(duneSlope,duneCrestHeight,10,color='black')
# ax35 = plt.subplot2grid((5,5),(0,4),rowspan=1,colspan=1)
# ax35.scatter(beachWidth,duneCrestHeight,10,color='black')
#
# ax36 = plt.subplot2grid((5,5),(1,1),rowspan=1,colspan=1)
# ax36.hist(duneCrestElevation)
# ax36.set_ylabel('Dune Crest Elevation')
# ax36.set_xlabel('Dune Crest Elevation')
# ax37 = plt.subplot2grid((5,5),(1,2),rowspan=1,colspan=1)
# ax37.scatter(duneToeElevation,duneCrestElevation,10,color='black')
# ax38 = plt.subplot2grid((5,5),(1,3),rowspan=1,colspan=1)
# ax38.scatter(duneSlope,duneCrestElevation,10,color='black')
# ax39 = plt.subplot2grid((5,5),(1,4),rowspan=1,colspan=1)
# ax39.scatter(beachWidth,duneCrestElevation,10,color='black')
#
# ax40 = plt.subplot2grid((5,5),(2,2),rowspan=1,colspan=1)
# ax40.hist(duneToeElevation)
# ax40.set_ylabel('Dune Toe Elevation')
# ax40.set_xlabel('Dune Toe Elevation')
# ax41 = plt.subplot2grid((5,5),(2,3),rowspan=1,colspan=1)
# ax41.scatter(duneSlope,duneToeElevation,10,color='black')
# ax42 = plt.subplot2grid((5,5),(2,4),rowspan=1,colspan=1)
# ax42.scatter(beachWidth,duneToeElevation,10,color='black')
#
# ax43 = plt.subplot2grid((5,5),(3,3),rowspan=1,colspan=1)
# ax43.hist(duneSlope)
# ax43.set_ylabel('Dune Slope')
# ax43.set_xlabel('Dune Slope')
# ax44 = plt.subplot2grid((5,5),(3,4),rowspan=1,colspan=1)
# ax44.scatter(beachWidth,duneSlope,10,color='black')
#
# ax45 = plt.subplot2grid((5,5),(4,4),rowspan=1,colspan=1)
# ax45.hist(beachWidth)
# ax45.set_ylabel('Beach Width')
# ax45.set_xlabel('Beach Width')
#
# data = np.empty((len(duneCrestHeight),5))
# data[:,0] = duneCrestHeight
# data[:,1] = duneCrestElevation
# data[:,2] = duneToeElevation
# data[:,3] = duneSlope
# data[:,4] = beachWidth
#
# import mda
# mdaSub = mda.MaxDiss_Simplified_NoThreshold(data,10,[0,1,2,3,4],[])
# import matplotlib.cm as cm
# colors = cm.rainbow(np.linspace(0, 1, 10))
# ax32.scatter(mdaSub[:,1],mdaSub[:,0],c=colors)#,color='red')
# ax33.scatter(mdaSub[:,2],mdaSub[:,0],c=colors)#,color='red')
# ax34.scatter(mdaSub[:,3],mdaSub[:,0],c=colors)#,color='red')
# ax35.scatter(mdaSub[:,4],mdaSub[:,0],c=colors)#,color='red')
# ax37.scatter(mdaSub[:,2],mdaSub[:,1],c=colors)#,color='red')
# ax38.scatter(mdaSub[:,3],mdaSub[:,1],c=colors)#,color='red')
# ax39.scatter(mdaSub[:,4],mdaSub[:,1],c=colors)#,color='red')
#
# ax41.scatter(mdaSub[:,3],mdaSub[:,2],c=colors)#,color='red')
# ax42.scatter(mdaSub[:,4],mdaSub[:,2],c=colors)#,color='red')
#
# ax44.scatter(mdaSub[:,4],mdaSub[:,3],c=colors)#,color='red')
#
#






# plt.figure()
# ax1 = plt.subplot2grid((1,1),(0,0),colspan=1,rowspan=1)
# # ax1.plot(x[0:20],duneY[0:20],label='Avg Dune Slope')
# ax1.plot(0,avgDuneCrestElevation,'ro',label='Avg Dune Crest Elevation')
# ax1.plot(duneToeX,avgDuneToeElevation,'bo',label='Avg Dune Toe Elevation')
# ax1.plot([0, duneToeX],[avgDuneCrestElevation,avgDuneToeElevation],'k--')
# ax1.plot(shorelineX,0,'ko',label='Avg Beach Width')
# ax1.plot(shorelineXPlusSTD,0,'go',label='SD of Beach Width')
# ax1.plot(shorelineXMinusSTD,0,'go')
# ax1.plot(shorelineXPlus2STD,0,'yo',label='Two SD of Beach Width')
# ax1.plot(shorelineXMinus2STD,0,'yo')
# ax1.plot([duneToeX,shorelineX],[avgDuneToeElevation,0],'k--')
# ax1.plot([duneToeX,shorelineXMinus2STD],[avgDuneToeElevation,0],'--',color=[0.5, 0.5, 0.5])
# ax1.plot([duneToeX,shorelineXMinusSTD],[avgDuneToeElevation,0],'--',color=[0.5, 0.5, 0.5])
# ax1.plot([duneToeX,shorelineXPlus2STD],[avgDuneToeElevation,0],'--',color=[0.5, 0.5, 0.5])
# ax1.plot([duneToeX,shorelineXPlusSTD],[avgDuneToeElevation,0],'--',color=[0.5, 0.5, 0.5])
#
# ax1.legend()


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n



dbfile = open('10km_morpho_parameters_and_profiles.pickle', 'rb')
profileData = pickle.load(dbfile)
dbfile.close()
#
# duneCrestElevation = []
# duneCrestHeight = []
# duneToeElevation = []
# beachWidth = []
# shoreline = []
# duneSlope = []
# beachSlope = []
# duneToeLoc = []
# bermLoc = []
rawX = []
rawZ = []
rawNum = []
counter = 1
for xx in profileData:
    print(xx)
    # duneCrestElevation = np.append(duneCrestElevation, profileData[xx]['DuneCrest_Elevation'])
    # duneCrestHeight = np.append(duneCrestHeight, profileData[xx]['DuneCrest_Height'])
    # duneToeElevation = np.append(duneToeElevation, profileData[xx]['DuneToe_Elevation'])
    # beachWidth = np.append(beachWidth, profileData[xx]['Beach_Width'])
    # shoreline = np.append(shoreline, profileData[xx]['Shoreline'])
    # duneSlope = np.append(duneSlope, profileData[xx]['Dune_Slope'])
    # beachSlope = np.append(beachSlope, profileData[xx]['Beach_Slope'])
    # duneToeLoc = np.append(duneToeLoc, profileData[xx]['DuneToe_loc'])
    # bermLoc = np.append(bermLoc, profileData[xx]['DuneToe_loc'])
    # rawX = np.append(rawX, profileData[xx]['dist_from_offshore'])
    # rawZ = np.append(rawZ, profileData[xx]['elevation'])

    rawX.append(profileData[xx]['dist_from_offshore'])
    rawZ.append(profileData[xx]['elevation'])
    rawNum.append(counter*np.ones(len(profileData[xx]['elevation'])))
    counter = counter + 1

# # zInterp = np.empty((np.shape(newX,)))
# for xx in range(len(xSubset)):

#
# plt.figure()
# for xx in range(len(rawX)):
#     plt.plot()


#
newX = np.arange(0,854,1)
plt.figure(figsize=(10,10))
allZs = []
allNums = []
for xx in range(len(rawX)):
    tempXSurvey = rawX[xx]
    tempZSurvey = rawZ[xx]
    indices1 = np.arange(10,499,1)
    indices2 = np.arange(500,1000,1)
    indices = np.hstack((indices1, indices2))
    for zz in range(len(indices)):
        tempIndex = indices[zz]
        f2 = interp1d(tempXSurvey[tempIndex], tempZSurvey[tempIndex], kind='linear')
        if zz == 0:
            zInterp = f2(newX)
        else:
            zInterp = np.vstack((zInterp,f2(newX)))

    zInterp2 = np.fliplr(zInterp)

    zAligned = -10 * np.ones((np.shape(zInterp2)))
    numAligned = np.ones((len(zInterp2),))*xx
    counter = 0
    for zz in range(len(zInterp2)):
        #tempI = zInterp2[zz,0]
        tempZs = zInterp2[zz,:]

        # indexI = np.where((np.abs(tempZs-2) == np.nanmin((np.abs(tempZs-2)))))
        indexI = np.where((tempZs>4.25) & (tempZs <5.25))


        if len(indexI[0]) > 0:
            # invertZ = tempZs[(indexI[0][-1]-100):]
            temp = np.diff(indexI)
            highInside = np.where((temp > 1))
            if len(highInside[0]) > 0:
                invertZ = tempZs[(indexI[0][highInside[0][0]]-100):]
            #else:
                #    invertZ = tempZs[(indexI[0][0] - 100):]
                plt.plot(newX, zInterp2[zz, :])#, 'k-')
                tempX = np.arange(0,len(invertZ),1)
                if invertZ[0] > -5:
                    if invertZ[-1] < -6.75:
                    #     if invertZ[-20] < -7.8:
                    #         if invertZ[670] < -5.9:
                    #             if invertZ[625] < -5.82:
                    #                 if invertZ[716] < -7:
                    #                     if invertZ[645] < -6.18:
                    #                         if invertZ[575] < -5.2:
                    #                             if invertZ[490] < -4.3:
                            # if counter == 0:
                            #     zAligned = invertZ
                            #     plt.plot(newX,zAligned)
                            # else:
                            #     zAligned = np.vstack((zAligned,invertZ))
                            #
                            # plt.plot(newX[indexI[0]],zInterp2[zz,indexI[0]],'ro')
                        lowerBathyInds = np.where((tempX > 475))
                        upperBathyInds = np.where((tempX <= 475))
                        if invertZ[lowerBathyInds[0][-1]] > -6.25:
                            xbot = np.array([tempX[lowerBathyInds[0][0]],tempX[lowerBathyInds[0][-1]]])
                            zbot = np.array([invertZ[lowerBathyInds[0][0]],-8])
                        else:
                            xbot = np.array([tempX[lowerBathyInds[0][0]],tempX[lowerBathyInds[0][-1]]])
                            zbot = np.array([invertZ[lowerBathyInds[0][0]],invertZ[lowerBathyInds[0][-1]]])
                        mb = np.polyfit(xbot, zbot, 1)
                        f = np.poly1d(mb)
                        newZ = f(tempX[lowerBathyInds[0]])
                        newInvert = np.hstack((invertZ[upperBathyInds[0]],newZ))

                        #moveAvg = moving_average(newInvert,10)


                                                        #moveMid = moving_average(invertZ,20)
                                                        #ending = moving_average(invertZ,40)

                                                        #together = np.hstack((moveAvg[0:200],moveMid[190:350]))
                                                        #together = np.hstack((together,ending[330:]))
                                                        #together = moving_average(together, 10)

                        zAligned[counter,0:len(newInvert)] = newInvert
                        numAligned[counter] = xx
                        plt.plot(newX,zAligned[counter,:])
                        counter = counter + 1

            del indexI

    allZs.append(zAligned[0:counter,:])
    allNums.append(numAligned[0:counter])


for xx in range(len(allZs)):
    if xx == 0:
        z = allZs[xx]
        numberS = allNums[xx]
    else:
        z = np.vstack((z,allZs[xx]))
        numberS = np.hstack((numberS,allNums[xx]))



zShort = z[:,0:726]
xShort = newX[0:726]


xOrig = np.flipud(np.abs(rawX[0][0][14:]-858))
xOrig = xOrig-xOrig[0]
x = np.ma.MaskedArray.filled(xOrig)
# xEnd = np.ma.MaskedArray.filled(xOrig[100:])
# xInitial = np.arange(0,50,4)
# xIM = np.arange(50,70,2)
# xMiddle = np.arange(70,150,1)
# xME = np.arange(150,200,2)
# xME2 = np.arange(200,270,4)
# xME3 = np.arange(275,400,8)
# xME4 = np.arange(410,720,12)
#
# x = np.hstack((xInitial,xIM))
# x = np.hstack((x,xMiddle))
# x = np.hstack((x,xME))
# x = np.hstack((x,xME2))
# x = np.hstack((x,xME3))
# x = np.hstack((x,xME4))

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
plt.savefig('allProfiles.png')
plt.close()
plt.close()
# from sklearn.decomposition import PCA

stdZ = np.std(zInt, axis=0)
demean = np.zeros((np.shape(zInt)))
normZ = np.zeros((np.shape(zInt)))
plt.figure()
for xx in range(len(demean)):
    demean[xx,:] = (zInt[xx, :] - meanZ)
    normZ[xx,:] = demean[xx,:] / stdZ
    plt.plot(x,demean[xx,:])


ipca = PCA(n_components=5)

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
for tt in uniqueNumbers:
    indexN = np.where((numberS == tt))
    temp = np.mean(PCs[indexN,0],axis=1)
    pcs1[int(tt)] = temp[0]
    temp2 = np.mean(PCs[indexN,1],axis=1)
    pcs2[int(tt)] = temp2[0]
    temp3 = np.mean(PCs[indexN,2],axis=1)
    pcs3[int(tt)] = temp3[0]
    temp4 = np.mean(PCs[indexN,3],axis=1)
    pcs4[int(tt)] = temp4[0]
    temp5 = np.mean(PCs[indexN,4],axis=1)
    pcs5[int(tt)] = temp5[0]

import datetime as dt
dates = np.asarray([dt.datetime(2010,11,1), dt.datetime(2011,6,1), dt.datetime(2011,11,1),dt.datetime(2012,6,1), dt.datetime(2012,11,1),
         dt.datetime(2013,6,1),dt.datetime(2014,6,1),dt.datetime(2015,6,1),
         dt.datetime(2016,6,1),dt.datetime(2016,10,1),dt.datetime(2017,7,1),dt.datetime(2018,5,1),dt.datetime(2019,4,1),
         dt.datetime(2019,8,1),dt.datetime(2019,11,1)])

import matplotlib.dates as mdates
formatter = mdates.DateFormatter("%Y-%m")

plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[0,:])
ax1.set_title('EOF 1')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax2.xaxis.set_major_formatter(formatter)
# ax2.plot(dates,pcs1,'k--')
#
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
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2.set_xlabel('Surveys')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF1NagsHead.png')
plt.close()

plt.figure(figsize=(12,8))
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
ax1.plot(x,EOFs[1,:])
ax1.set_title('EOF 2')
ax1.set_xlabel('Cross-shore (m)')
ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,1],10,numberS,cmap='tab20')
ax2.set_xlabel('Transect #')
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
ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,2],10,numberS,cmap='tab20')
ax2.set_xlabel('Transect #')
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
ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,3],10,numberS,cmap='tab20')
ax2.set_xlabel('Transect #')
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
ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,4],10,numberS,cmap='tab20')
ax2.set_xlabel('Transect #')
ax2.set_ylabel('PCs')
plt.tight_layout()
plt.savefig('EOF5NagsHead.png')
plt.close()





plt.figure(figsize=(12,8))
ax0 = plt.subplot2grid((3,6),(0,0),rowspan=1,colspan=3)
ax0.plot([0,725],[0.6,0.6],'b--',label='MHW')
ax0.set_xlim([0, 725])

for xx in range(len(zInt)):
    ax0.plot(x, zInt[xx,:], color=[0.5,0.5,0.5])
ax0.plot(x,meanZ,'k-',linewidth=2,label='Average Profile')
ax0.set_ylabel('Elevation (m)')
ax0.legend()
ax1 = plt.subplot2grid((3,6),(1,0),rowspan=1,colspan=3)
ax1.plot(x,EOFs[0,:])
# ax1.text(400,0.15,'EOF#1')
ax1.set_ylabel('EOF#1')
ax1.text(220,0.17,'Mode explains variability related to sub-aerial beach width')
ax1.set_xlim([0, 725])
# ax1.set_title('EOF 1')
# ax1.set_xlabel('Cross-shore (m)')
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
ax1b = plt.subplot2grid((3,6),(1,3),rowspan=1,colspan=3)

ax1b.boxplot(data, positions=pos, widths = 30)
ax1b.set_xlim([0, (x2-x1).days ])
# ax1b.set_xticklabels(['{}/{}'.format(ff.month,ff.year) for ff in dates], rotation=45)
ax1b.set_xticklabels([''])
ax1b.plot(pos,pcs1,'--',color=[0.2,0.2,0.2])
ax1b.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
# ax1b.set_xlabel('Surveys')
ax1b.set_ylabel('PCs')
ax1b.set_title('Beach width grows immediately after nourishment, and subsequently erodes')

ax2 = plt.subplot2grid((3,6),(2,0),rowspan=1,colspan=3)
ax2.plot(x,EOFs[1,:])
# ax2.text(400,0.17,'EOF#2')
ax2.set_xlabel('Cross-shore (m)')
ax2.set_ylabel('EOF#2')
ax2.text(140,0.18,'Mode explains variability of dune and surfzone trough/bar')
ax2.set_xlim([0, 725])



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
ax2b = plt.subplot2grid((3,6),(2,3),rowspan=1,colspan=3)

ax2b.boxplot(data, positions=pos, widths = 30)
ax2b.set_xlim([0, (x2-x1).days ])
ax2b.set_xticklabels(['{}/{}'.format(ff.month,ff.year) for ff in dates], rotation=45)
ax2b.plot(pos,pcs2,'--',color=[0.2,0.2,0.2])
ax2b.plot([0, (x2-x1).days],[0,0],'-',color=[0.5,0.5,0.5],linewidth=0.5)
# ax2.scatter(np.arange(0,len(PCs[:,0]),1),PCs[:,0],10,numberS,cmap='tab20')
# ax2.set_xlabel('Transect #')
ax2b.set_xlabel('Surveys')
ax2b.set_ylabel('PCs')
ax2b.set_title('Dunes have grown following the initial nourishment')

# ax3 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
# ax3.plot(x,EOFs[2,:])
# ax3.set_title('EOF 3')
# ax3.set_xlabel('Cross-shore (m)')





X_projects = ipca.inverse_transform(PCs)
zOrig = zInt



num = 1800
prof1 = meanZ + EOFs[0,:]*PCs[num,0]
prof2 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1]
prof3 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2]
prof4 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3]
prof5 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]
# prof6 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]
# prof7 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]
# prof8 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]+ EOFs[7,:]*PCs[num,7]


plt.figure()
plt.plot(x,meanZ,label='mean profile')
plt.plot(x,zOrig[num,:],'k-',linewidth=2,label='original profile')
# plt.plot(x,X_projects[num,:]+meanZ,'r-', label='Full EOF reconstruction')
plt.plot(x,prof1,label='1 EOF reconstruction')
plt.plot(x,prof2,label='2 EOFs reconstruction')
plt.plot(x,prof3,label='3 EOFs reconstruction')
plt.plot(x,prof4,label='4 EOFs reconstruction')
plt.plot(x,prof5,'r-',label='5 EOFs reconstruction')
# plt.plot(x,prof6,'b-',label='6 EOFs reconstruction')
# plt.plot(x,prof7,'g-',label='7 EOFs reconstruction')
# plt.plot(x,prof8,'m-',label='8 EOFs reconstruction')
plt.plot(x,X_projects[num,:],'k--',label='Inverse')
plt.legend()


PCs1 = PCs[:,0]
PCs2 = PCs[:,1]
PCs3 = PCs[:,2]
PCs4 = PCs[:,3]
PCs5 = PCs[:,4]
# PCs6 = PCs[:,5]
# PCs7 = PCs[:,6]
# PCs8 = PCs[:,7]

plt.figure(figsize=(14,8))
ax31 = plt.subplot2grid((5,5),(0,0),rowspan=1,colspan=1)
ax31.hist(PCs1)
ax31.set_ylabel('EOF1')
ax31.set_xlabel('EOF1')
xy = np.vstack([PCs2,PCs1])
z = gaussian_kde(xy)(xy)
ax32 = plt.subplot2grid((5,5),(0,1),rowspan=1,colspan=1)
ax32.scatter(PCs2,PCs1,5,c=z)

ax33 = plt.subplot2grid((5,5),(0,2),rowspan=1,colspan=1)
xy = np.vstack([PCs3,PCs1])
z = gaussian_kde(xy)(xy)
ax33.scatter(PCs3,PCs1,5,c=z)
ax34 = plt.subplot2grid((5,5),(0,3),rowspan=1,colspan=1)
xy = np.vstack([PCs4,PCs1])
z = gaussian_kde(xy)(xy)
ax34.scatter(PCs4,PCs1,5,c=z)
ax35 = plt.subplot2grid((5,5),(0,4),rowspan=1,colspan=1)
xy = np.vstack([PCs5,PCs1])
z = gaussian_kde(xy)(xy)
ax35.scatter(PCs5,PCs1,5,c=z)

ax36 = plt.subplot2grid((5,5),(1,1),rowspan=1,colspan=1)
ax36.hist(PCs2)
ax36.set_ylabel('EOF2')
ax36.set_xlabel('EOF2')
ax37 = plt.subplot2grid((5,5),(1,2),rowspan=1,colspan=1)
xy = np.vstack([PCs3,PCs2])
z = gaussian_kde(xy)(xy)
ax37.scatter(PCs3,PCs2,5,c=z)
ax38 = plt.subplot2grid((5,5),(1,3),rowspan=1,colspan=1)
xy = np.vstack([PCs4,PCs2])
z = gaussian_kde(xy)(xy)
ax38.scatter(PCs4,PCs2,5,c=z)
ax39 = plt.subplot2grid((5,5),(1,4),rowspan=1,colspan=1)
xy = np.vstack([PCs5,PCs2])
z = gaussian_kde(xy)(xy)
ax39.scatter(PCs5,PCs2,5,c=z)

ax40 = plt.subplot2grid((5,5),(2,2),rowspan=1,colspan=1)
ax40.hist(PCs3)
ax40.set_ylabel('EOF3')
ax40.set_xlabel('EOF3')
ax41 = plt.subplot2grid((5,5),(2,3),rowspan=1,colspan=1)
xy = np.vstack([PCs4,PCs3])
z = gaussian_kde(xy)(xy)
ax41.scatter(PCs4,PCs3,5,c=z)
ax42 = plt.subplot2grid((5,5),(2,4),rowspan=1,colspan=1)
xy = np.vstack([PCs5,PCs3])
z = gaussian_kde(xy)(xy)
ax42.scatter(PCs5,PCs3,5,c=z)

ax43 = plt.subplot2grid((5,5),(3,3),rowspan=1,colspan=1)
ax43.hist(PCs4)
ax43.set_ylabel('EOF4')
ax43.set_xlabel('EOF4')
ax44 = plt.subplot2grid((5,5),(3,4),rowspan=1,colspan=1)
xy = np.vstack([PCs5,PCs4])
z = gaussian_kde(xy)(xy)
ax44.scatter(PCs5,PCs4,5,c=z)

ax45 = plt.subplot2grid((5,5),(4,4),rowspan=1,colspan=1)
ax45.hist(PCs5)
ax45.set_ylabel('EOF5')
ax45.set_xlabel('EOF5')

plt.tight_layout()
# plt.savefig('EOFJointProbabilities.png')
# plt.close()




dataCop = []
for xx in range(len(PCs1)):
    dataCop.append(list([PCs1[xx],PCs2[xx],PCs3[xx],PCs4[xx],PCs5[xx]])) #,PCs6[xx],PCs7[xx],PCs8[xx]]))

from copula import pyCopula
cop = pyCopula.Copula(dataCop)
samples = cop.gendata(100)

samplePCs1 = [item[0] for item in samples]
samplePCs2 = [item[1] for item in samples]
samplePCs3 = [item[2] for item in samples]
samplePCs4 = [item[3] for item in samples]
samplePCs5 = [item[4] for item in samples]
#samplePCs6 = [item[4] for item in samples]
#samplePCs7 = [item[4] for item in samples]
# samplePCs8 = [item[4] for item in samples]



























plt.figure(figsize=(12,12))
hypoProfiles = []
for xx in range(len(samplePCs1)):

    # prof1 = meanZ + EOFs[0,:]*PCs[num,0]
    # prof2 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1]
    # prof3 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2]
    # prof4 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3]
    # prof5 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]
    # prof6 = meanZ + EOFs[0,:]*samplePCs1[xx] + EOFs[1,:]*samplePCs2[xx] + EOFs[2,:]*samplePCs3[xx] \
    #         + EOFs[3,:]*samplePCs4[xx] + EOFs[4,:]*samplePCs5[xx]+ EOFs[5,:]*samplePCs6[xx]
    # prof6 = meanZ + EOFs[0,:]*samplePCs1[xx] + EOFs[1,:]*samplePCs2[xx] + EOFs[2,:]*samplePCs3[xx] \
    #         + EOFs[3,:]*samplePCs4[xx] + EOFs[4,:]*samplePCs5[xx]+ EOFs[5,:]*samplePCs6[xx]
    prof5 = meanZ + EOFs[0,:]*samplePCs1[xx] + EOFs[1,:]*samplePCs2[xx] + EOFs[2,:]*samplePCs3[xx] \
            + EOFs[3,:]*samplePCs4[xx] + EOFs[4,:]*samplePCs5[xx] #+ EOFs[5,:]*samplePCs6[xx] + EOFs[6,:]*samplePCs7[xx] + EOFs[7,:]*samplePCs8[xx]
    # prof7 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]
    # prof8 = meanZ + EOFs[0,:]*PCs[num,0] + EOFs[1,:]*PCs[num,1] + EOFs[2,:]*PCs[num,2] + EOFs[3,:]*PCs[num,3] + EOFs[4,:]*PCs[num,4]+ EOFs[5,:]*PCs[num,5]+ EOFs[6,:]*PCs[num,6]+ EOFs[7,:]*PCs[num,7]
    hypoProfiles.append(prof5)
    plt.plot(x,prof5)

plt.xlabel('cross-shore (m)')
plt.title('100 Hypothetical Profiles')
# plt.savefig('hypotheticalProfiles.png')



# pickleName = 'hypotheticalProfiles1Oct2020.pickle'
# output = {}
# output['x'] = x
# output['profiles'] = hypoProfiles
# import pickle
# with open(pickleName, 'wb') as f:
#     pickle.dump(output,f)


del rawX
del rawZ
del rawNum
del allZs
del allNums
del tempXSurvey
del tempZSurvey

# dbfile = open('first_10storms_profile.pickle', 'rb')
# first10Profiles = pickle.load(dbfile)
# dbfile.close()

asfg

dbfile = open('150storms_profile.pickle', 'rb')
xbProfiles = pickle.load(dbfile)
dbfile.close()
xXbeach = xbProfiles['dist_from_backshore']



plt.figure()
c1 = 0
c2 = 0
for ii in xbProfiles:
    # if ii == 'storm_46':
    #     print('skipping 46')
    # elif ii == 'storm_97':
    #     print('skipping 97')
    # elif ii == 'dist_from_backshore':
    #     print('dont print x')
    # else:
    ax = plt.subplot2grid((15,10),(c2,c1),rowspan=1, colspan=1)
    ax.plot(xXbeach,xbProfiles[ii][0,:])
    ax.plot(xXbeach,xbProfiles[ii][1,:])
    ax.set_xlim([75, 400])
    ax.set_ylim([-5, 7.5])
    ax.text(190,4,'{}'.format(ii))
    if c1 > 0:
        ax.set_yticks([])
    if c2 < 14:
        ax.set_xticks([])

    if c1 < 9:
        c1 = c1 + 1
    elif c1 == 9:
        c1 = 0
        c2 = c2+1
plt.tight_layout()

oldProfiles = np.nan * np.ones((150,168))
newProfiles = np.nan * np.ones((150,168))

c = 0
for ii in xbProfiles:
    # if ii == 'storm_46':
    #     print('skipping 46')
    # elif ii == 'storm_97':
    #     print('skipping 97')
    # elif ii == 'storm_15':
    #     print('skipping 15')
    # elif ii == 'storm_44':
    #     print('skipping 44')
    # elif ii == 'storm_70':
    #     print('skipping 70')
    # # elif ii == 'storm_97':
    # #     print('skipping 97')
    if ii == 'dist_from_backshore':
        print('dont print x')
    else:
        oldProfiles[c,:] = xbProfiles[ii][0,0:168]
        newProfiles[c,:] = xbProfiles[ii][1,0:168]
        c = c +1




# dbfile = open('hypotheticalProfiles23Sept2020.pickle', 'rb')
# hypos = pickle.load(dbfile)
# dbfile.close()
#
#
#
# plt.figure()
#
# plt.plot(xXbeach[0:168],oldProfiles[0,:])
# plt.plot(x,hypos['profiles'][0])

oldEOFscores = ipca.transform(oldProfiles)
newEOFscores = ipca.transform(newProfiles)

# [0:34, 100, 35:46, 101, 47:87, 102, 88:100, 103:153]

pickleName = 'hypotheticalStorms.pickle'
dbfile = open(pickleName, 'rb')
storms = pickle.load(dbfile)
dbfile.close()
predictors = np.vstack((storms['hypotheticalStorms'][0:34,:],storms['hypotheticalStorms'][100,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][35:46,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][101,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][47:87,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][102,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][88:100,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][103:153,:]))
predictors = np.hstack((predictors,oldEOFscores))



plt.figure()
ax1 = plt.subplot2grid((1,5),(0,0),rowspan=1,colspan=1)
ax1.plot(oldEOFscores[:,0],newEOFscores[:,0],'.')#,8,predictors[:,0])
ax1.set_title('EOF 1')
ax1.plot([-12,12],[-12,12],'k--')
ax1.set_xlim([-12, 12])
ax1.set_ylabel('Post-Storm PCs')
ax1.set_xlabel('Pre-Storm PCs')

ax2 = plt.subplot2grid((1,5),(0,1),rowspan=1,colspan=1)
ax2.plot(oldEOFscores[:,1],newEOFscores[:,1],'.')
ax2.set_title('EOF 2')
ax2.plot([-9,9],[-9,9],'k--')
ax2.set_xlim([-9, 9])
ax2.set_xlabel('Pre-Storm PCs')

ax3 = plt.subplot2grid((1,5),(0,2),rowspan=1,colspan=1)
ax3.plot(oldEOFscores[:,2],newEOFscores[:,2],'.')
ax3.set_title('EOF 3')
ax3.plot([-6,6],[-6,6],'k--')
ax3.set_xlim([-6, 6])
ax3.set_xlabel('Pre-Storm PCs')

ax4 = plt.subplot2grid((1,5),(0,3),rowspan=1,colspan=1)
ax4.plot(oldEOFscores[:,3],newEOFscores[:,3],'.')
ax4.set_title('EOF 4')
ax4.plot([-6,6],[-6,6],'k--')
ax4.set_xlim([-6, 6])
ax4.set_xlabel('Pre-Storm PCs')

ax5 = plt.subplot2grid((1,5),(0,4),rowspan=1,colspan=1)
ax5.plot(oldEOFscores[:,4],newEOFscores[:,4],'.')
ax5.set_title('EOF 5')
ax5.plot([-4.5,4.5],[-4.5,4.5],'k--')
ax5.set_xlim([-4.5, 4.5])
ax5.set_xlabel('Pre-Storm PCs')





predictands = newEOFscores

output = dict()
output['predictors'] = predictors
output['predictands'] = predictands
import scipy.io
scipy.io.savemat('profileTrainingTemp2.mat',output)







# pickleName = 'eofsForHypotheticalProfiles23Sept2020.pickle'
# output = {}
# output['oldEOFscores'] = oldEOFscores
# output['oldProfiles'] = oldProfiles
# output['newEOFscores'] = newEOFscores
# output['newProfiles'] = newProfiles
# import pickle
# with open(pickleName, 'wb') as f:
#     pickle.dump(output,f)


# Last 5 ad hoc done from matlab

# predictions = np.asarray([[10.1612, -1.8977, -1.2696, -2.0518, 0.3854],
#                [4.4429, -1.7726, 3.8883, -2.8549, 0.0429],
#                [-9.8108, -6.1893, 4.4404, -0.0924, 0.0374],
#                [4.5846, 3.78909, -3.6755, 1.7636, 1.2592],
#                [3.8274, 0.7086, 4.2271, -2.7818, -1.1348]])
#

predictions1 = np.asarray([8.24051552882147,0.651288166601663,2.06792065670337,-0.00441327333890396,2.14055378197120,
                           -5.76508708091749,10.8516586914893,1.02964039135815,9.22626712043547,1.46434713659600,
                           8.14427096889788,1.26750871482987,1.23465059657112,0.167281545957887,3.42137740055490,
                           3.77665432354694,-0.444723362564683,-0.971675102120118,-4.12864488588518,-5.43697620252006,
                           -0.926732011185326,-1.02984620141484,4.30264012300377,-0.318420090547170,-0.414590571922267,
                           1.06596498488804,2.55819996047097,3.12353739338365,0.855173810334266,7.07598403626834])

predictions2 = np.asarray([1.82152636342608,1.79504352490753,-0.995002690251203,-4.18162723251791,-4.57456304947201,
                           1.20604789336448,1.97348065424107,6.40219203552991,3.87403126114057,3.2603876492145,
                           -3.87606308950695,-1.64772838720118,-3.67258475203456,-4.09667655303335,
                           2.66766204368773,0.628536557590298,-2.96730199524186,-2.03108192801391,-4.06408517705855,
                           4.70924401742090,-3.73185609919827,-0.0406094329869013,-2.63781967882023,-0.959135184166478,
                           -3.04127935216056,-0.604906507345651,-1.21654033722951,7.08052869684351,0.861388518883308,
                           0.431003179955089])

predictions3 = np.asarray([0.383803168191418,-0.998750164850993,-1.13915817661372,-1.67333686353148,1.55523146781139,
                           -2.62559192573100,-0.444883431359792,-1.52501368711512,0.216531216241704,-4.42021327884469,
                           0.455192783909792,0.570188502382415,1.45183716772740,2.07553075746929,-2.32857541603716,
                           -3.40909110027151,-0.0975201768365936,-1.46134680089794,-0.00311293203712193,
                           -2.44666781876326,-1.04899242209801,2.27298998698457,-1.36286148749277,0.386978918657746,
                           1.77873531930211,-4.29798819914894,-2.04276050120061,0.181972152972804,1.70275116291807,
                           -3.26675629961543])


predictions4 = np.asarray([-1.71246964883942,-2.36752974367223,2.79106745060750,-0.391727058179780,-1.12783492948230,
                           2.98241344050615,-2.89810143651943,4.69811785195051,0.864505323670072,0.692936494515054,
                           -1.75443224190770,0.000746427433099134,0.683456359047676,1.86506272489373,
                           -0.0452797586814782,3.25184787406419,-1.36842344295981,-0.944183026973077,-1.99909290039668,
                           -0.210684588819978,-3.28719881335771,-1.35932826321093,-2.50118756266531,-1.95140650243688,
                           1.68042338368387,-2.04760601096167,-1.62435142181038,-1.50439176087597,-2.56315566304863,-2.54925177479166])

predictions5 = np.asarray([1.56974194212655,2.54005289467676,-1.15111673041522,0.745811340515135,-0.354911689126626,
                           3.65301457253994,-0.568104685286717,4.03349792885840,0.790298142881911,-0.897332124779644,
                           -1.77684764152538,0.241039901656537,-0.0403385646118353,-0.371633432713014,0.785509826932261,
                           0.578685508328669,1.22948756479417,1.55844720671083,2.18881797758020,2.62810185524506,
                           0.478924617605723,2.52715319586365,-1.69165868607795,0.616283543187773,0.842365676790348,
                           1.30802968597117,-0.392473406146465,0.402179951644244,3.46141350154159,-0.0746910831616405])




predictions = np.asarray([[8.61478630154569,2.29554815436805,	0.118753113016612,	-1.71724347677536,	2.01333775766301],
[0.709109830176786,1.69103949204515,	-0.999308881672605,	-2.30651686263261,	2.64265037891111],
[2.21121018385551,-1.05875887465013,	-1.54324551573926,	2.64342041964814,	-1.14702234592147],
[0.108337507465734,-4.23274190196759,	-1.73193322117835,	-0.509022204932414,	0.888790109328465],
[2.20073086898406,-4.36306175483212,1.49119812321273,	-1.07479975732173,	-0.305801050823949],
[-5.38727634950308,1.18836602464844,	-2.34513397550795,	2.99519918888951,	3.37657086624662],
[11.2415997962125,2.52035804963716,	-0.585140115124694,	-2.77832602187350,	-0.530325996074424],
[1.01396973679870,6.43873624465687,	-2.22112652989338,	4.75874787376950,	4.10303840826472],
[9.72985526373801,4.06310087401961,	-0.0516506670672193,	1.00065069335779,	1.32492915944874],
[1.31567579341810,3.34199797273961,	-3.76002055346865,	0.652483576023726,	-1.35286589583143],
[8.50823391563771,-3.58087702166188,	0.351184821737250,	-1.56804942945799,	-2.14831359124264],
[1.15145937318416,-1.63839964195704,	0.556290619515988,	-0.217610115486571,	0.161404450470295],
[0.722370982884582,-3.83047632501440,	1.47686275324442,	0.577095713760294,	0.0928045012970151],
[-0.256757532141367,-3.88204060513738,	2.15497176269994,	1.64870246383254,	-0.220103098810994],
[3.21437713621113,2.45931968906742,	-2.57539533424476,	-0.0783114705977258,	0.928364247638447],
[3.85137304567536,0.367025699178320,	-3.79253467428765,	3.42695264686239,	0.726947552137049],
[-0.231120610943439,-2.74933487200510,	-0.231345865529368,	-1.34967370297286,	1.32453983446370],
[-0.913692040227585,-2.12144913339000,	-1.10976565227377,	-0.921611580146056,	1.66320221855209],
[-4.15128272948918,-4.46376276508577,	0.121938013325229,	-1.99381459944728,	2.20045367574750],
[-5.44321054715006,4.02039104975187,	-2.42451512938627,	-0.207423072889936,	2.56913829467659],
[-1.22052622386801,-3.73996146984148,	-0.806296415021845,	-3.26477553320633,	0.515583651393068],
[-0.900576028458174,-0.0453090832938398,	2.58732505095286,	-1.37574152135121,	2.58226928944121],
[4.06804366897089,-2.85857908900791,	-1.00871704531999,	-2.46473142891612,	-1.74647962684592],
[-0.113585366220498,-0.970942535905689,	0.632981898034370,	-1.91613676808130,	0.543113965558594],
[-0.865468974833355,-2.64974999906964,	1.76971952522414,	1.58525718873889,	1.05030002824940],
[0.989489865536071,-0.599970396835085,	-3.79773583314060,	-2.13278400556520,	1.34140994293427],
[2.65853413209522,-1.27847744602918,	-2.42411081536225,	-1.66169895248516,	-0.390557174614277],
[2.90810120852464,6.40921136534763,	0.428599520537487,	-1.61615450092294,	0.291409440093079],
[0.989669140776911,0.803185886208453,	2.09028999349598,	-2.39944290310130,	3.26115478442508],
[6.98256428409842,0.600888087935170,	-3.41329350071531,	-2.57875982776195,	-0.0518385657166410]])


plt.figure(figsize=(16,8))
calProfiles = []
c = 0
c2 = 0
for xx in range(len(predictions1)):
    prof5 = meanZ + EOFs[0,:]*predictions[xx,0] + EOFs[1,:]*predictions[xx,1] + EOFs[2,:]*predictions[xx,2] \
            + EOFs[3,:]*predictions[xx,3] + EOFs[4,:]*predictions[xx,4]
    # prof5 = meanZ + EOFs[0,:]*predictions1[xx] + EOFs[1,:]*predictions2[xx] + EOFs[2,:]*predictions3[xx] \
    #         + EOFs[3,:]*predictions4[xx] + EOFs[4,:]*predictions5[xx]
    calProfiles.append(prof5)
    ax = plt.subplot2grid((3,10),(c2,c),rowspan=1,colspan=1)
    ax.plot(x,oldProfiles[xx+120,:],label='Pre-storm ')
    ax.plot(x,newProfiles[xx+120,:],label='Post-storm Xbeach')
    ax.plot(x,prof5,label='Post-storm GPR')
    ax.set_xlim([70, 300])
    ax.set_ylim([-3.5,7])
    ax.set_title('Trial {}'.format(xx+120))
    if c == 9:
        c = 0
        c2 = c2+1
    else:
        c = c + 1
ax.legend()





plt.figure(figsize=(12,6))
showIndex = np.array([121,123,127,132,144])
c2 = 0
c = 0
for xx in range(len(predictions1)):
    getInd = showIndex[xx]-120
    prof5 = meanZ + EOFs[0,:]*predictions[getInd,0] + EOFs[1,:]*predictions[getInd,1] + EOFs[2,:]*predictions[getInd,2] \
            + EOFs[3,:]*predictions[getInd,3] + EOFs[4,:]*predictions[getInd,4]
    # prof5 = meanZ + EOFs[0,:]*predictions1[xx] + EOFs[1,:]*predictions2[xx] + EOFs[2,:]*predictions3[xx] \
    #         + EOFs[3,:]*predictions4[xx] + EOFs[4,:]*predictions5[xx]
    calProfiles.append(prof5)
    ax = plt.subplot2grid((1,5),(0,c),rowspan=1,colspan=1)
    ax.plot(x,oldProfiles[getInd+120,:],label='Pre-storm ',color='black')
    ax.plot(x,newProfiles[getInd+120,:],label='Post-storm Xbeach',color='orange')
    ax.plot(x,prof5,label='Post-storm GPR',color='green')
    ax.set_xlim([70, 300])
    ax.set_ylim([-3.5,7])
    ax.set_title('Trial {}'.format(xx+121))
    c = c + 1
    # if c == 9:
    #     c = 0
    #     c2 = c2+1
    # else:
    #     c = c + 1
ax.legend()

#
# pickleName = 'hypotheticalProfiles23Sept2020.pickle'
# output = {}
# output['x'] = x
# output['profiles'] = hypoProfiles
# import pickle
# with open(pickleName, 'wb') as f:
#     pickle.dump(output,f)
#
# for xx in range(10):#range(len(zAligned)):
#
#
#
#     invertZ = tempZs[indexI[0][-1]:]
#     zAligned[xx,0:len(invertZ)] = invertZ
#     #plt.plot(newX,zAligned[xx,:])


#
#
# del profileData
# xArray = np.empty((np.shape(rawX[0][0])))
# zArray = np.empty((np.shape(rawX[0][0])))
#
# for xx in rawX:
#     for yy in xx:
#         xArray = np.vstack((xArray,yy))
#
# for xx in rawZ:
#     for yy in xx:
#         zArray = np.vstack((zArray,yy))
#
#
# xSubset = xArray[1::50,:]
# zSubset = zArray[1::50,:]
#
# newX = np.arange(0,854,4)
# del xArray
# del zArray
# # zSubset = np.fliplr(zSubset)
#
# # zInterp = np.empty((np.shape(newX,)))
# for xx in range(len(xSubset)):
#     f2 = interp1d(xSubset[xx,:], zSubset[xx,:], kind='linear')
#     if xx == 0:
#         zInterp = f2(newX)
#     else:
#         zInterp = np.vstack((zInterp,f2(newX)))
#
#
# zInterp = np.fliplr(zInterp)
#
# plt.figure()
# zAligned = -10 * np.ones((np.shape(zInterp)))
# for xx in range(10):#range(len(zAligned)):
#     tempI = zInterp[xx,0]
#     tempZs = zInterp[xx,:]
#     indexI = np.where((tempZs == tempI))
#     invertZ = tempZs[indexI[0][-1]:]
#     zAligned[xx,0:len(invertZ)] = invertZ
#     #plt.plot(newX,zAligned[xx,:])
#     plt.plot(newX,zInterp[xx,:])