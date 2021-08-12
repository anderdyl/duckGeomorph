
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
from shapely.geometry import LineString



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

# rawX = np.empty((1,307))
# rawZ = np.empty((1,307))
#
# allNumbers = []
# counter = 1
# for xx in profileData:
#     print(xx)
#
#     if len(xx) > 5:
#         if xx == 'dist':
#             print('did''t add the distance vector')
#         else:
#             rawZ = np.vstack((rawZ,profileData[xx]))
#             surveyNum = counter*np.ones((1,len(profileData[xx])))
#             allNumbers.append(surveyNum)
#             counter = counter + 1
#
#
# for xx in range(len(allNumbers)):
#     if xx == 0:
#         numberSurvs = allNumbers[xx]
#     else:
#         numberSurvs = np.hstack((numberSurvs,allNumbers[xx]))
#
#
# rawZ = rawZ[1:,:]
# rawX = profileData['dist']

rawZ = profileData['2019_08_']
rawX = profileData['dist']


# newX = np.arange(0,1199,1)
# # plt.figure(figsize=(10,10))
# allZs = []
# #allNums = []
# allXs = []
# for xx in range(len(rawZ)):
#     tempXSurvey = rawX
#     tempZSurvey = rawZ[xx,:]
#     indexN = np.where((tempZSurvey < -50))
#     if len(indexN[0]) > 0:
#         tempZSurvey[indexN[0]] = 1.5*np.ones((np.shape(indexN[0])))
#     f2 = interp1d(tempXSurvey, tempZSurvey, kind='linear')
#     zInterp = f2(newX)
#     zAligned = 1.5 * np.ones((850,))
#     numAligned = np.ones((len(zInterp),))*xx
#     # lets go for the dune crest...
#     indexI = np.where((zInterp == np.max(zInterp)))
#     if np.max(zInterp) > 2:
#         invertZ = (zInterp[(indexI[0][0] - 30):])
#         backX = newX[(indexI[0][0] - 30)]
#         allXs.append(backX)
#         zAligned[0:849] = invertZ[0:849]
#         tempX = np.arange(0, 849, 1)
#         extendedX = np.arange(0,849,1)
#         allZs.append(zAligned)
#


newX = np.arange(0,1199,1)
allXs = []
allZs = []
allNums = []
for xx in range(len(rawZ)):
    tempXSurvey = rawX
    tempZSurvey = rawZ[xx,:]
    indexN = np.where((tempZSurvey < -40))
    if len(indexN[0]) > 0:
        tempZSurvey[indexN[0]] = 1.5*np.ones((np.shape(indexN[0])))
    #tempNumber = numberSurvs[0][xx]
    f2 = interp1d(tempXSurvey, tempZSurvey, kind='linear')
    #if xx == 0:
    zInterpOG = f2(newX)

    zInterp = moving_average(zInterpOG, 20)

    zAligned = 4.5 * np.ones((850,))

    numAligned = np.ones((len(zInterp),))*xx
    # lets go for the dune crest...
    indexI = np.where((zInterp == np.max(zInterp)))
    indexI = np.where((zInterp > 3.5) & (zInterp < 4.25))
    if xx > 1235 and xx < 1255:
        indexI = np.where((zInterp > 3.25) & (zInterp < 4.25))
    if xx > 1460 and xx < 1485:
        indexI = np.where((zInterp > 3.15) & (zInterp < 4.25))
    if xx > 1469 and xx < 1477:
        indexI = np.where((zInterp > 3) & (zInterp < 4.25))

    if len(indexI[0]) > 0:
        if np.max(zInterp) > 2.5:
            invertZ = (zInterp[(indexI[0][-1] - 70):])
            backX = newX[(indexI[0][-1] - 70)]
            allXs.append(backX)
            # zAligned[0:len(invertZ)] = invertZ
            zAligned[0:849] = invertZ[0:849]
            #plt.plot(newX[0:len(invertZ)], invertZ)#, 'k-')
            # tempX = np.arange(0,len(invertZ),1)
            tempX = np.arange(0, 849, 1)
            extendedX = np.arange(0,849,1)
            allZs.append(zAligned)
    else:
        print('skipped profile # {}'.format(xx))
            #allNums.append(tempNumber)





xInitial = np.arange(0,70,7)
xInit2 = np.arange(70,90,2)
xIM = np.arange(90,100,1)
xIM2 = np.arange(100,180,0.5)
xIM3 = np.arange(180,230,1)
xMiddle = np.arange(230,260,2)
xMiddle2 = np.arange(260,440,10)
xME = np.arange(440,840,20)

x = np.hstack((xInitial,xInit2))
x = np.hstack((x,xIM))
x = np.hstack((x,xIM2))
x = np.hstack((x,xIM3))
x = np.hstack((x,xMiddle))
x = np.hstack((x,xMiddle2))
x = np.hstack((x,xME))



for xx in range(len(allZs)):
    if xx == 0:
        z = allZs[xx]
        f2 = interp1d(newX[0:len(z)], z, kind='linear')
        # if xx == 0:
        zDownscaled = f2(x)
    else:
        z = allZs[xx]
        f2 = interp1d(newX[0:len(z)], z, kind='linear')
        z2 = f2(x)
        zDownscaled = np.vstack((zDownscaled,z2))



plt.figure(figsize=(10,10))
for xx in range(len(zDownscaled)):
    plt.plot(x,zDownscaled[xx,:])

meanZ = np.mean(zDownscaled,axis=0)
plt.plot(x,meanZ,'k-',linewidth=3)
plt.xlabel('Cross-shore (m)')
plt.ylabel('Elevation (m)')
plt.title('All Profiles (black = mean)')
plt.savefig('allProfilesAligned.png')
plt.close()



import pickle

dbfile = 'nourishmentProfileEOFs2.pickle'

with open(r'nourishmentProfileEOFs2.pickle', "rb") as input_file:
    outputEOFs = pickle.load(input_file)

meanZ = outputEOFs['meanZ']
EOFs = outputEOFs['EOFs']
PCs = outputEOFs['PCs']
variance = outputEOFs['variance']
varianceRatio = outputEOFs['varianceRatio']
cumulativeVar = outputEOFs['cumulativeVar']
zInt = outputEOFs['zInt']
x = outputEOFs['x']
dist = x
originalIPCA = PCA(n_components=9)

# PCs = ipca.fit_transform(normZ)
origPCs = originalIPCA.fit_transform(zInt)

newPCs = originalIPCA.transform(zDownscaled)
newProfiles = originalIPCA.inverse_transform(newPCs)


initPCs = dict()
initPCs['initPCs'] = newPCs

import scipy.io
scipy.io.savemat('init2019PCs.mat',initPCs)



# num = 0
# plt.figure()
# plt.plot(x,zDownscaled[num,:],'k-',linewidth=2,label='original profile')
# plt.plot(x,newProfiles[num,:],'k--',label='Inverse')
# plt.legend()
# plt.title('EOF Reconstruction Example Profile {}'.format(num))
# plt.show()

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

rmseProfiles = np.empty((1600,))
for hhh in range(len(zDownscaled)):
    rmseProfiles[hhh] = rmse(newProfiles[hhh],zDownscaled[hhh])


for hhh in range(len(rmseProfiles)):
    if hhh == 0:
        newXProfiles = x + allXs[hhh]+10
    else:
        newXs = x + allXs[hhh]+10
        newXProfiles = np.vstack((newXProfiles,newXs))




nourishmentFitPickle = 'nourishment201908FitToEOFs2.pickle'
outputNourishmentFit = {}
outputNourishmentFit['newXProfiles'] = newXProfiles
outputNourishmentFit['newProfiles'] = newProfiles

with open(nourishmentFitPickle,'wb') as f:
    pickle.dump(outputNourishmentFit, f)




#
# with open(r'nourishment201908FitToEOFs.pickle', "rb") as input_file:
#     inputProfs = pickle.load(input_file)
#
# oldXProfiles = inputProfs['newXProfiles']
# oldProfiles = inputProfs['newProfiles']
#
#
# i = 10
# plt.figure()
# plt.plot(oldXProfiles[i],oldProfiles[i])
# plt.plot(newXProfiles[i],newProfiles[i])
#


asdfg

profileNumber = 400


dataPred = scipy.io.loadmat('storm1Fits.mat')
dist1 = dataPred['dist1']
a1 = dataPred['a1']
b1 = dataPred['b1']

zeroLine = np.zeros((np.shape(newProfiles[0,:])))
surrogates1 = np.zeros((np.shape(newProfiles)))

surrogates = np.zeros((np.shape(newProfiles)))
fig2 = plt.figure()
c = 0
c2 = 0
for hh in range(len(surrogates)):
    profIndex = hh
    if hh == profileNumber:
        p5 = plt.subplot2grid((1, 1), (c, c2), rowspan=1, colspan=1)
        p5.plot(dist, newProfiles[profIndex, :], color=[0, 0, 0],label='August Survey')

    f = zeroLine
    g = newProfiles[hh, :]
    newx = np.arange(0, 100, 1)
    first_line = LineString(np.column_stack((dist, f)))
    second_line = LineString(np.column_stack((dist, g)))
    intersection = first_line.intersection(second_line)

    if intersection.geom_type == 'Multi3Point':
        print('shit, multipoint')

    elif intersection.geom_type == 'Point':
        x, y = intersection.xy
        zeroCrossingX = x[0]

    f2 = zeroLine
    g2 = zDownscaled[hh,:]
    first_line2 = LineString(np.column_stack((dist, f2)))
    second_line2 = LineString(np.column_stack((dist, g2)))
    intersection2 = first_line.intersection(second_line2)

    if intersection.geom_type == 'Point':
        x2, y2 = intersection2.xy
        if x2 > x:
            zeroCrossingX = x2[0]



    # Lets make a new line
    # everthing shoreward of the top of the scarp
    ogDuneInd = np.where((dist < (zeroCrossingX - dist1[0][profIndex])))
    scarpXDune = dist[ogDuneInd]
    scarpZDune = newProfiles[profIndex,ogDuneInd]

    # finding the points between the scarp and the zero point...
    # xIndex = np.where((dist <= zeroCrossingX) & (dist > (zeroCrossingX - dist1[0][profIndex])))

    # find the points between the scarp and the zero point + 100 m
    xIndex = np.where((dist <= zeroCrossingX+100) & (dist > (zeroCrossingX - dist1[0][profIndex])))

    z2 = a1[0][profIndex] * np.power(dist[xIndex]-dist[xIndex[0][0]], b1[0][profIndex])
    newx2 = dist[xIndex]-dist[xIndex[0][0]]
    scarpXBeach = newx2 + zeroCrossingX - dist1[0][profIndex]
    scarpZBeach = z2 + newProfiles[profIndex,xIndex[0][0]]

    # # finding negatives and making them zzero
    # findNegatives = np.where((scarpZBeach < 0))
    # scarpZBeach[findNegatives] = scarpZBeach[findNegatives]*0

    diff = scarpZBeach - newProfiles[profIndex, xIndex]
    keepers = np.where(diff[0] <= 0)

    scarpXBeach = scarpXBeach[keepers[0]]
    scarpZBeach = scarpZBeach[keepers[0]]

    if scarpZBeach[-1] < -2:
        alter = np.where((scarpZBeach < 0))
        below = np.where((newProfiles[profIndex,:] < -2))
        xlinear = [scarpXBeach[alter[0][0]],dist[below[0][0]]]
        zlinear = [0,-2]
        flinear = interp1d(xlinear,zlinear)
        offshore = np.where((dist > scarpXBeach[alter[0][0]]))
        xflat = dist[offshore[0][0]:below[0][0]]
        zflat = flinear(xflat)
        keepers2 = np.where((scarpZBeach>0))
        scarpXBeach = scarpXBeach[keepers2[0]]
        scarpZBeach = scarpZBeach[keepers2[0]]

        ogSubaqueous = np.where((dist > xflat[-1]))
        scarpXSubaqueous = dist[ogSubaqueous]
        scarpZSubaqueous = newProfiles[profIndex,ogSubaqueous]

        predScarp1X = np.hstack((scarpXDune,scarpXBeach))
        predScarp15X = np.hstack((predScarp1X,xflat))
        predScarp2X = np.hstack((predScarp15X,scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune[0], scarpZBeach))
        predScarp15Z = np.hstack((predScarp1Z,zflat))
        predScarp2Z = np.hstack((predScarp15Z, scarpZSubaqueous[0]))

    else:
        ogSubaqueous = np.where((dist > scarpXBeach[-1]))
        scarpXSubaqueous = dist[ogSubaqueous]
        scarpZSubaqueous = newProfiles[profIndex,ogSubaqueous]

        predScarp1X = np.hstack((scarpXDune,scarpXBeach))
        predScarp2X = np.hstack((predScarp1X,scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune[0], scarpZBeach))
        predScarp2Z = np.hstack((predScarp1Z, scarpZSubaqueous[0]))



    if hh == profileNumber:
        plt.plot(predScarp2X,predScarp2Z,label='Surrogate Storm 1')

    finterp = interp1d(predScarp2X,predScarp2Z,'cubic')

    surrogateProf = finterp(dist)
    smoothed = moving_average(surrogateProf,7)
    extended = np.hstack((surrogateProf[0:3], smoothed))
    extended2 = np.hstack((extended,surrogateProf[-4:-1]))
    surrogates[hh,:] = extended2 #
    surrogates1[hh,:] = surrogateProf #

    if hh == profileNumber:
        if c2 == 9:
            c = c+1
            c2 = 0
        else:
            c2 = c2+1



storm1newPCs = originalIPCA.transform(surrogates)
storm1newProfiles = originalIPCA.inverse_transform(storm1newPCs)
p5.plot(dist, storm1newProfiles[profileNumber, :], '--', color=[0.7, 0.7, 0.7])

storm1PCs = dict()
storm1PCs['storm1PCs'] = storm1newPCs
import scipy.io
scipy.io.savemat('storm1PCs.mat',storm1PCs)







dataPred2 = scipy.io.loadmat('storm2Fits.mat')
dist2 = dataPred2['dist2']
a2 = dataPred2['a2']
b2 = dataPred2['b2']


surrogates2 = np.zeros((np.shape(storm1newProfiles)))

c = 0
c2 = 0
counting = 0
for hh in range(len(surrogates2)):
    profIndex = hh


    f = zeroLine
    g = surrogates[hh, :]
    newx = np.arange(0, 100, 1)
    first_line = LineString(np.column_stack((dist, f)))
    second_line = LineString(np.column_stack((dist, g)))
    intersection = first_line.intersection(second_line)

    if intersection.geom_type == 'Multi3Point':
        print('shit, multipoint')
        # plt.plot(*LineString(intersection).xy, 'o')
        # x, y = LineString(intersection).xy
        # # zeroCrossingX[i] = x[0]
        # plt.xlim([x[0] - 100, x[0] + 100])
        # plt.ylim([-3, 5])
    elif intersection.geom_type == 'Point':
        x, y = intersection.xy
        zeroCrossingX = x[0]

    f2 = zeroLine
    g2 = storm1newProfiles[hh,:]
    first_line2 = LineString(np.column_stack((dist, f2)))
    second_line2 = LineString(np.column_stack((dist, g2)))
    intersection2 = first_line.intersection(second_line2)

    if intersection.geom_type == 'Point':
        x2, y2 = intersection2.xy
        if x2 > x:
            zeroCrossingX = x2[0]
            print('choosing EOF shoreline: {}'.format(counting))
            countering =+ 1
    # z = scarpPointsMid[i, 2] * np.power(newx, scarpPointsMid[i, 3])

    # Lets make a new line

    ogDuneInd = np.where((dist < (zeroCrossingX - dist2[0][profIndex])))

    scarpXDune = dist[ogDuneInd]
    scarpZDune = storm1newProfiles[profIndex,ogDuneInd]

    # finding the points between the scarp and the zero point...
    # xIndex = np.where((dist <= zeroCrossingX) & (dist > (zeroCrossingX - dist1[0][profIndex])))

    # find the points between the scarp and the zero point + 100 m
    xIndex = np.where((dist <= zeroCrossingX+100) & (dist > (zeroCrossingX - dist2[0][profIndex])))

    z2 = a2[0][profIndex] * np.power(dist[xIndex]-dist[xIndex[0][0]], b2[0][profIndex])
    newx2 = dist[xIndex]-dist[xIndex[0][0]]
    scarpXBeach = newx2 + zeroCrossingX - dist2[0][profIndex]
    scarpZBeach = z2 + storm1newProfiles[profIndex,xIndex[0][0]]

    # # finding negatives and making them zzero
    # findNegatives = np.where((scarpZBeach < 0))
    # scarpZBeach[findNegatives] = scarpZBeach[findNegatives]*0

    diff = scarpZBeach - storm1newProfiles[profIndex, xIndex]
    keepers = np.where(diff[0] <= 0)

    scarpXBeach = scarpXBeach[keepers[0]]
    scarpZBeach = scarpZBeach[keepers[0]]

    if scarpZBeach[-1] < -2:
        alter = np.where((scarpZBeach < 0))
        below = np.where((storm1newProfiles[profIndex,:] < -2))
        xlinear = [scarpXBeach[alter[0][0]],dist[below[0][0]]]
        zlinear = [0,-2]
        flinear = interp1d(xlinear,zlinear)
        offshore = np.where((dist > scarpXBeach[alter[0][0]]))
        xflat = dist[offshore[0][0]:below[0][0]]
        zflat = flinear(xflat)
        keepers2 = np.where((scarpZBeach>0))
        scarpXBeach = scarpXBeach[keepers2[0]]
        scarpZBeach = scarpZBeach[keepers2[0]]

        ogSubaqueous = np.where((dist > xflat[-1]))
        scarpXSubaqueous = dist[ogSubaqueous]
        scarpZSubaqueous = storm1newProfiles[profIndex,ogSubaqueous]

        predScarp1X = np.hstack((scarpXDune,scarpXBeach))
        predScarp15X = np.hstack((predScarp1X,xflat))
        predScarp2X = np.hstack((predScarp15X,scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune[0], scarpZBeach))
        predScarp15Z = np.hstack((predScarp1Z,zflat))
        predScarp2Z = np.hstack((predScarp15Z, scarpZSubaqueous[0]))

    else:
        ogSubaqueous = np.where((dist > scarpXBeach[-1]))
        scarpXSubaqueous = dist[ogSubaqueous]
        scarpZSubaqueous = storm1newProfiles[profIndex,ogSubaqueous]

        predScarp1X = np.hstack((scarpXDune,scarpXBeach))
        predScarp2X = np.hstack((predScarp1X,scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune[0], scarpZBeach))
        predScarp2Z = np.hstack((predScarp1Z, scarpZSubaqueous[0]))

    if hh == profileNumber:
        plt.plot(predScarp2X,predScarp2Z,label='Surrogate Storm 2')

    finterp = interp1d(predScarp2X,predScarp2Z,'cubic')

    surrogateProf = finterp(dist)

    smoothed = moving_average(surrogateProf,7)
    extended = np.hstack((surrogateProf[0:3], smoothed))
    extended2 = np.hstack((extended,surrogateProf[-4:-1]))
    #surrogates2[hh,:] = extended2 #
    surrogates2[hh,:] = surrogateProf





storm2newPCs = originalIPCA.transform(surrogates2)
storm2newProfiles = originalIPCA.inverse_transform(storm2newPCs)
p5.plot(dist, storm2newProfiles[profileNumber, :], '--', color=[0.7, 0.7, 0.7])

storm2PCs = dict()
storm2PCs['storm2PCs'] = storm2newPCs
import scipy.io
scipy.io.savemat('storm2PCs.mat',storm2PCs)






dataPred3 = scipy.io.loadmat('storm3Fits.mat')
dist3 = dataPred3['dist3']
a3 = dataPred3['a3']
b3 = dataPred3['b3']


surrogates3 = np.zeros((np.shape(storm2newProfiles)))

c = 0
c2 = 0
for hh in range(len(surrogates3)):
    profIndex = hh


    f = zeroLine
    g = surrogates2[hh, :]
    newx = np.arange(0, 100, 1)
    first_line = LineString(np.column_stack((dist, f)))
    second_line = LineString(np.column_stack((dist, g)))
    intersection = first_line.intersection(second_line)

    if intersection.geom_type == 'Multi3Point':
        print('shit, multipoint')
        # plt.plot(*LineString(intersection).xy, 'o')
        # x, y = LineString(intersection).xy
        # # zeroCrossingX[i] = x[0]
        # plt.xlim([x[0] - 100, x[0] + 100])
        # plt.ylim([-3, 5])
    elif intersection.geom_type == 'Point':
        x, y = intersection.xy
        zeroCrossingX = x[0]


    f2 = zeroLine
    g2 = storm2newProfiles[hh,:]
    first_line2 = LineString(np.column_stack((dist, f2)))
    second_line2 = LineString(np.column_stack((dist, g2)))
    intersection2 = first_line.intersection(second_line2)

    if intersection.geom_type == 'Point':
        x2, y2 = intersection2.xy
        if x2 > x:
            zeroCrossingX = x2[0]

    # z = scarpPointsMid[i, 2] * np.power(newx, scarpPointsMid[i, 3])

    # Lets make a new line

    ogDuneInd = np.where((dist < (zeroCrossingX - dist3[0][profIndex])))

    scarpXDune = dist[ogDuneInd]
    scarpZDune = storm2newProfiles[profIndex,ogDuneInd]

    # finding the points between the scarp and the zero point...
    # xIndex = np.where((dist <= zeroCrossingX) & (dist > (zeroCrossingX - dist1[0][profIndex])))

    # find the points between the scarp and the zero point + 100 m
    xIndex = np.where((dist <= zeroCrossingX+100) & (dist > (zeroCrossingX - dist3[0][profIndex])))

    z2 = a3[0][profIndex] * np.power(dist[xIndex]-dist[xIndex[0][0]], b3[0][profIndex])
    newx2 = dist[xIndex]-dist[xIndex[0][0]]
    scarpXBeach = newx2 + zeroCrossingX - dist3[0][profIndex]
    scarpZBeach = z2 + storm2newProfiles[profIndex,xIndex[0][0]]

    # # finding negatives and making them zzero
    # findNegatives = np.where((scarpZBeach < 0))
    # scarpZBeach[findNegatives] = scarpZBeach[findNegatives]*0

    diff = scarpZBeach - storm2newProfiles[profIndex, xIndex]
    keepers = np.where(diff[0] <= 0)

    scarpXBeach = scarpXBeach[keepers[0]]
    scarpZBeach = scarpZBeach[keepers[0]]

    if scarpZBeach[-1] < -2:
        alter = np.where((scarpZBeach < 0))
        below = np.where((storm2newProfiles[profIndex,:] < -2))
        xlinear = [scarpXBeach[alter[0][0]],dist[below[0][0]]]
        zlinear = [0,-2]
        flinear = interp1d(xlinear,zlinear)
        offshore = np.where((dist > scarpXBeach[alter[0][0]]))
        xflat = dist[offshore[0][0]:below[0][0]]
        zflat = flinear(xflat)
        keepers2 = np.where((scarpZBeach>0))
        scarpXBeach = scarpXBeach[keepers2[0]]
        scarpZBeach = scarpZBeach[keepers2[0]]

        ogSubaqueous = np.where((dist > xflat[-1]))
        scarpXSubaqueous = dist[ogSubaqueous]
        scarpZSubaqueous = storm2newProfiles[profIndex,ogSubaqueous]

        predScarp1X = np.hstack((scarpXDune,scarpXBeach))
        predScarp15X = np.hstack((predScarp1X,xflat))
        predScarp2X = np.hstack((predScarp15X,scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune[0], scarpZBeach))
        predScarp15Z = np.hstack((predScarp1Z,zflat))
        predScarp2Z = np.hstack((predScarp15Z, scarpZSubaqueous[0]))

    else:
        ogSubaqueous = np.where((dist > scarpXBeach[-1]))
        scarpXSubaqueous = dist[ogSubaqueous]
        scarpZSubaqueous = storm2newProfiles[profIndex,ogSubaqueous]

        predScarp1X = np.hstack((scarpXDune,scarpXBeach))
        predScarp2X = np.hstack((predScarp1X,scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune[0], scarpZBeach))
        predScarp2Z = np.hstack((predScarp1Z, scarpZSubaqueous[0]))

    if hh == profileNumber:
        plt.plot(predScarp2X,predScarp2Z,label='Surrogate Storm 3')

    finterp = interp1d(predScarp2X,predScarp2Z,'cubic')

    surrogateProf = finterp(dist)
    smoothed = moving_average(surrogateProf,7)
    extended = np.hstack((surrogateProf[0:3], smoothed))
    extended2 = np.hstack((extended,surrogateProf[-4:-1]))
    #surrogates3[hh,:] = extended2 #
    surrogates3[hh,:] = surrogateProf





storm3newPCs = originalIPCA.transform(surrogates3)
storm3newProfiles = originalIPCA.inverse_transform(storm3newPCs)
p5.plot(dist, storm3newProfiles[profileNumber, :], '--', color=[0.7, 0.7, 0.7])

storm3PCs = dict()
storm3PCs['storm3PCs'] = storm3newPCs
import scipy.io
scipy.io.savemat('storm3PCs.mat',storm3PCs)











dataPred4 = scipy.io.loadmat('storm4Fits.mat')
dist4 = dataPred4['dist4']
a4 = dataPred4['a4']
b4 = dataPred4['b4']


surrogates4 = np.zeros((np.shape(storm3newProfiles)))
surrogatesOG = np.zeros((np.shape(storm3newProfiles)))
c = 0
c2 = 0
for hh in range(len(surrogates3)):
    profIndex = hh


    f = zeroLine
    g = surrogates3[hh, :]
    newx = np.arange(0, 100, 1)
    first_line = LineString(np.column_stack((dist, f)))
    second_line = LineString(np.column_stack((dist, g)))
    intersection = first_line.intersection(second_line)

    if intersection.geom_type == 'Multi3Point':
        print('shit, multipoint')
        # plt.plot(*LineString(intersection).xy, 'o')
        # x, y = LineString(intersection).xy
        # # zeroCrossingX[i] = x[0]
        # plt.xlim([x[0] - 100, x[0] + 100])
        # plt.ylim([-3, 5])
    elif intersection.geom_type == 'Point':
        x, y = intersection.xy
        zeroCrossingX = x[0]


    f2 = zeroLine
    g2 = storm3newProfiles[hh,:]
    first_line2 = LineString(np.column_stack((dist, f2)))
    second_line2 = LineString(np.column_stack((dist, g2)))
    intersection2 = first_line.intersection(second_line2)

    if intersection.geom_type == 'Point':
        x2, y2 = intersection2.xy
        if x2 > x:
            zeroCrossingX = x2[0]

    # z = scarpPointsMid[i, 2] * np.power(newx, scarpPointsMid[i, 3])

    # Lets make a new line

    ogDuneInd = np.where((dist < (zeroCrossingX - dist4[0][profIndex])))

    scarpXDune = dist[ogDuneInd]
    scarpZDune = storm3newProfiles[profIndex,ogDuneInd]

    scarpXDuneOG = dist[ogDuneInd]
    scarpZDuneOG = zDownscaled[profIndex,ogDuneInd]
    # finding the points between the scarp and the zero point...
    # xIndex = np.where((dist <= zeroCrossingX) & (dist > (zeroCrossingX - dist1[0][profIndex])))

    # find the points between the scarp and the zero point + 100 m
    xIndex = np.where((dist <= zeroCrossingX+100) & (dist > (zeroCrossingX - dist4[0][profIndex])))

    z2 = a4[0][profIndex] * np.power(dist[xIndex]-dist[xIndex[0][0]], b4[0][profIndex])
    newx2 = dist[xIndex]-dist[xIndex[0][0]]
    scarpXBeach = newx2 + zeroCrossingX - dist4[0][profIndex]
    scarpZBeach = z2 + storm3newProfiles[profIndex,xIndex[0][0]]

   # z2 = a4[0][profIndex] * np.power(dist[xIndex]-dist[xIndex[0][0]], b4[0][profIndex])
   # newx2 = dist[xIndex]-dist[xIndex[0][0]]
    scarpXBeachOG = newx2 + zeroCrossingX - dist4[0][profIndex]
    scarpZBeachOG = z2 + zDownscaled[profIndex,xIndex[0][0]]

    # # finding negatives and making them zzero
    # findNegatives = np.where((scarpZBeach < 0))
    # scarpZBeach[findNegatives] = scarpZBeach[findNegatives]*0

    diff = scarpZBeach - storm3newProfiles[profIndex, xIndex]
    keepers = np.where(diff[0] <= 0)
    scarpXBeach = scarpXBeach[keepers[0]]
    scarpZBeach = scarpZBeach[keepers[0]]

    diffOG = scarpZBeachOG - zDownscaled[profIndex, xIndex]
    keepersOG = np.where(diffOG[0] <= 0)
    scarpXBeachOG = scarpXBeachOG[keepersOG[0]]
    scarpZBeachOG = scarpZBeachOG[keepersOG[0]]

    if scarpZBeach[-1] < -2:
        alter = np.where((scarpZBeach < 0))
        below = np.where((storm3newProfiles[profIndex,:] < -2))
        xlinear = [scarpXBeach[alter[0][0]],dist[below[0][0]]]
        zlinear = [0,-2]
        flinear = interp1d(xlinear,zlinear)
        offshore = np.where((dist > scarpXBeach[alter[0][0]]))
        xflat = dist[offshore[0][0]:below[0][0]]
        zflat = flinear(xflat)
        keepers2 = np.where((scarpZBeach>0))
        scarpXBeach = scarpXBeach[keepers2[0]]
        scarpZBeach = scarpZBeach[keepers2[0]]

        ogSubaqueous = np.where((dist > xflat[-1]))
        scarpXSubaqueous = dist[ogSubaqueous]
        scarpZSubaqueous = storm3newProfiles[profIndex,ogSubaqueous]

        predScarp1X = np.hstack((scarpXDune,scarpXBeach))
        predScarp15X = np.hstack((predScarp1X,xflat))
        predScarp2X = np.hstack((predScarp15X,scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune[0], scarpZBeach))
        predScarp15Z = np.hstack((predScarp1Z,zflat))
        predScarp2Z = np.hstack((predScarp15Z, scarpZSubaqueous[0]))

    else:
        ogSubaqueous = np.where((dist > scarpXBeach[-1]))
        scarpXSubaqueous = dist[ogSubaqueous]
        scarpZSubaqueous = storm3newProfiles[profIndex,ogSubaqueous]

        predScarp1X = np.hstack((scarpXDune,scarpXBeach))
        predScarp2X = np.hstack((predScarp1X,scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune[0], scarpZBeach))
        predScarp2Z = np.hstack((predScarp1Z, scarpZSubaqueous[0]))




    if scarpZBeachOG[-1] < -2:
        alter = np.where((scarpZBeachOG < 0))
        below = np.where((zDownscaled[profIndex,:] < -2))
        xlinear = [scarpXBeachOG[alter[0][0]],dist[below[0][0]]]
        zlinear = [0,-2]
        flinear = interp1d(xlinear,zlinear)
        offshore = np.where((dist > scarpXBeachOG[alter[0][0]]))
        xflat = dist[offshore[0][0]:below[0][0]]
        zflat = flinear(xflat)
        keepers2OG = np.where((scarpZBeachOG>0))
        scarpXBeachOG = scarpXBeachOG[keepers2OG[0]]
        scarpZBeachOG = scarpZBeachOG[keepers2OG[0]]

        ogSubaqueousOG = np.where((dist > xflat[-1]))
        scarpXSubaqueousOG = dist[ogSubaqueousOG]
        scarpZSubaqueousOG = zDownscaled[profIndex,ogSubaqueousOG]

        predScarp1XOG = np.hstack((scarpXDuneOG,scarpXBeachOG))
        predScarp15XOG = np.hstack((predScarp1XOG,xflat))
        predScarp2XOG = np.hstack((predScarp15XOG,scarpXSubaqueousOG))
        predScarp1ZOG = np.hstack((scarpZDuneOG[0], scarpZBeachOG))
        predScarp15ZOG = np.hstack((predScarp1ZOG,zflat))
        predScarp2ZOG = np.hstack((predScarp15ZOG, scarpZSubaqueousOG[0]))

    else:
        ogSubaqueousOG = np.where((dist > scarpXBeachOG[-1]))
        scarpXSubaqueousOG = dist[ogSubaqueousOG]
        scarpZSubaqueousOG = zDownscaled[profIndex,ogSubaqueousOG]

        predScarp1XOG = np.hstack((scarpXDuneOG,scarpXBeachOG))
        predScarp2XOG = np.hstack((predScarp1XOG,scarpXSubaqueousOG))
        predScarp1ZOG = np.hstack((scarpZDuneOG[0], scarpZBeachOG))
        predScarp2ZOG = np.hstack((predScarp1ZOG, scarpZSubaqueousOG[0]))




    if hh == profileNumber:
        plt.plot(predScarp2X,predScarp2Z,label='Surrogate Storm 4')
        plt.plot(predScarp2XOG,predScarp2ZOG,label='Surrogate on OG')
        plt.legend()

    finterp = interp1d(predScarp2X,predScarp2Z,'cubic')
    finterp2 = interp1d(predScarp2XOG,predScarp2ZOG,'cubic')

    surrogateProf = finterp(dist)
    surrogateProfOG = finterp2(dist)
    smoothed = moving_average(surrogateProf,7)
    extended = np.hstack((surrogateProf[0:3], smoothed))
    extended2 = np.hstack((extended,surrogateProf[-4:-1]))
    surrogates4[hh,:] = extended2 #
    surrogatesOG[hh,:] = surrogateProfOG
    #surrogates4[hh,:] = surrogateProf


storm4newPCs = originalIPCA.transform(surrogates4)
storm4newProfiles = originalIPCA.inverse_transform(storm4newPCs)
p5.plot(dist, storm4newProfiles[profileNumber, :], '--', color=[0.7, 0.7, 0.7])

storm4PCs = dict()
storm4PCs['storm4PCs'] = storm4newPCs
import scipy.io
scipy.io.savemat('storm4PCs.mat',storm4PCs)






















dbfile = open('20km_16surveys_fullprofiles_edited.pickle', 'rb')

profileData = pickle.load(dbfile)
dbfile.close()
rawZ = profileData['2019_11_']
rawX = profileData['dist']


newX = np.arange(0,1199,1)
allXs2 = []
allZs = []
allNums = []
for xx in range(len(rawZ)):
    tempXSurvey = rawX
    tempZSurvey = rawZ[xx,:]
    indexN = np.where((tempZSurvey < -40))
    if len(indexN[0]) > 0:
        tempZSurvey[indexN[0]] = 1.5*np.ones((np.shape(indexN[0])))
    #tempNumber = numberSurvs[0][xx]
    f2 = interp1d(tempXSurvey, tempZSurvey, kind='linear')
    #if xx == 0:
    zInterpOG = f2(newX)

    zInterp = moving_average(zInterpOG, 20)

    zAligned = 4.5 * np.ones((850,))

    numAligned = np.ones((len(zInterp),))*xx
    # lets go for the dune crest...
    indexI = np.where((zInterp == np.max(zInterp)))
    indexI = np.where((zInterp > 3.5) & (zInterp < 4.25))
    if xx > 1235 and xx < 1255:
        indexI = np.where((zInterp > 3.25) & (zInterp < 4.25))
    if xx > 1460 and xx < 1485:
        indexI = np.where((zInterp > 3.15) & (zInterp < 4.25))
    if xx > 1469 and xx < 1477:
        indexI = np.where((zInterp > 3) & (zInterp < 4.25))

    if len(indexI[0]) > 0:
        if np.max(zInterp) > 2.5:
            invertZ = (zInterp[(indexI[0][-1] - 70):])
            backX = newX[(indexI[0][-1] - 70)]
            allXs2.append(backX)
            # zAligned[0:len(invertZ)] = invertZ
            zAligned[0:849] = invertZ[0:849]
            #plt.plot(newX[0:len(invertZ)], invertZ)#, 'k-')
            # tempX = np.arange(0,len(invertZ),1)
            tempX = np.arange(0, 849, 1)
            extendedX = np.arange(0,849,1)
            allZs.append(zAligned)
    else:
        print('skipped profile # {}'.format(xx))
            #allNums.append(tempNumber)





xInitial = np.arange(0,70,7)
xInit2 = np.arange(70,90,2)
xIM = np.arange(90,100,1)
xIM2 = np.arange(100,180,0.5)
xIM3 = np.arange(180,230,1)
xMiddle = np.arange(230,260,2)
xMiddle2 = np.arange(260,440,10)
xME = np.arange(440,840,20)

x = np.hstack((xInitial,xInit2))
x = np.hstack((x,xIM))
x = np.hstack((x,xIM2))
x = np.hstack((x,xIM3))
x = np.hstack((x,xMiddle))
x = np.hstack((x,xMiddle2))
x = np.hstack((x,xME))



for xx in range(len(allZs)):
    if xx == 0:
        z = allZs[xx]
        f2 = interp1d(newX[0:len(z)], z, kind='linear')
        # if xx == 0:
        zDownscaled2 = f2(x)
    else:
        z = allZs[xx]
        f2 = interp1d(newX[0:len(z)], z, kind='linear')
        z2 = f2(x)
        zDownscaled2 = np.vstack((zDownscaled2,z2))


laterPCs = originalIPCA.transform(zDownscaled2)
laterProfiles = originalIPCA.inverse_transform(laterPCs)

p5.plot(dist, laterProfiles[profileNumber, :], '--', color=[0.5, 0.5, 0.5],label='November Survey')
plt.legend()






surroPickle = 'surrogateOf4Storms.pickle'
output = {}
output['x'] = x
output['surrogates4'] = surrogatesOG
output['surrogates3'] = surrogates3
output['surrogates2'] = surrogates2
output['surrogates1'] = surrogates1
output['eofProfileAfterStorm1'] = storm1newProfiles
output['eofProfileAfterStorm2'] = storm2newProfiles
output['eofProfileAfterStorm3'] = storm3newProfiles
output['surveyInitial'] = zDownscaled
output['surveyFinal'] = zDownscaled2
output['backXs'] = allXs

import pickle
with open(surroPickle,'wb') as f:
    pickle.dump(output, f)





plt.figure()
num =399
plt.plot(x,surrogatesOG[num,:],label='surrogate')
plt.plot(x,zDownscaled[num,:],label='pre')
plt.plot(x,zDownscaled2[num,:],label='post')

plt.xlim([0,350])
plt.ylim([-5,5])
plt.legend()

asdfg







from shapely.geometry import Polygon

actualArea = list()
predictedArea = list()

for hh in range(1599):

    # Difference between the EOF eroded volumes volumes
    x_y_curveTest1 = list()
    x_y_curveTest2 = list()

    z1 = zDownscaled2[hh, :]
    z2 = zDownscaled[hh, :]

    findAbove35 = np.where((z1 > 3.05))
    z1[0:findAbove35[0][-1]] = z1[0:findAbove35[0][-1]] * 0 + 3.5
    findAbove352 = np.where((z2 > 3.05))
    z2[0:findAbove352[0][-1]] = z2[0:findAbove352[0][-1]] * 0 + 3.5

    findNegatives1 = np.where((z1 < 0))
    z1[findNegatives1] = z1[findNegatives1] * 0
    findNegatives2 = np.where((z2 < 0))
    z2[findNegatives2] = z2[findNegatives2] * 0


    finder = np.where((z1 < z2))

    for n in range(len(finder[0])):
        x_y_curveTest1.append(list((dist[n], z1[n])))
        x_y_curveTest2.append(list((dist[n], z2[n])))

    polygon_points = []  # creates a empty list where we will append the points to create the polygon

    for xyvalue in x_y_curveTest1:
        polygon_points.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 1

    for xyvalue in x_y_curveTest2[::-1]:
        polygon_points.append([xyvalue[0], xyvalue[
            1]])  # append all xy points for curve 2 in the reverse order (from last point to first point)

    for xyvalue in x_y_curveTest1[0:1]:
        polygon_points.append(
            [xyvalue[0], xyvalue[1]])  # append the first point in curve 1 again, to it "closes" the polygon

    polygon = Polygon(polygon_points)
    area = polygon.area
    actualArea.append(area)


    # Difference between the EOF eroded volumes volumes
    x_y_curveTest1 = list()
    x_y_curveTest2 = list()

    z1 = surrogates4[hh, :]
    z2 = zDownscaled[hh, :]

    findAbove35 = np.where((z1 > 3.25))
    z1[0:findAbove35[0][-1]] = z1[0:findAbove35[0][-1]] * 0 + 3.5
    findAbove352 = np.where((z2 > 3.25))
    z2[0:findAbove352[0][-1]] = z2[0:findAbove352[0][-1]] * 0 + 3.5

    findNegatives1 = np.where((z1 < 0))
    z1[findNegatives1] = z1[findNegatives1] * 0
    findNegatives2 = np.where((z2 < 0))
    z2[findNegatives2] = z2[findNegatives2] * 0


    for n in range(len(finder[0])):
        x_y_curveTest1.append(list((dist[n], z1[n])))
        x_y_curveTest2.append(list((dist[n], z2[n])))

    polygon_points = []  # creates a empty list where we will append the points to create the polygon

    for xyvalue in x_y_curveTest1:
        polygon_points.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 1

    for xyvalue in x_y_curveTest2[::-1]:
        polygon_points.append([xyvalue[0], xyvalue[
            1]])  # append all xy points for curve 2 in the reverse order (from last point to first point)

    for xyvalue in x_y_curveTest1[0:1]:
        polygon_points.append(
            [xyvalue[0], xyvalue[1]])  # append the first point in curve 1 again, to it "closes" the polygon

    polygon = Polygon(polygon_points)
    area = polygon.area
    predictedArea.append(area)

# # (area)
# p5.text(zeroCrossing[0][profIndex] - 25, 3.5, '{} m^2'.format(np.round(area * 10) / 10))
# # fdu = integrate.simps()

plt.figure()
plt.plot(actualArea,predictedArea,'.')
plt.xlim([-25,175])
plt.ylim([-25,125])
plt.grid()



predArray = np.asarray(predictedArea)
actuArray = np.asarray(actualArea)
rootMeanSquared = rmse(predArray,actuArray)
errors = np.abs(predArray-actuArray)
medianError = np.median(errors)






differences = np.zeros((1599,))
differences2 = np.zeros((1599,))

c = 0
c2 = 0
for hh in range(1599):
    profIndex = hh

    f = zeroLine
    g = surrogates4[hh, :]
    newx = np.arange(0, 100, 1)
    first_line = LineString(np.column_stack((dist, f)))
    second_line = LineString(np.column_stack((dist, g)))
    intersection = first_line.intersection(second_line)

    if intersection.geom_type == 'Multi3Point':
        print('shit, multipoint')
    elif intersection.geom_type == 'Point':
        x, y = intersection.xy
        zeroCrossingX = x[0]

    f2 = zeroLine
    g2 = newProfiles[hh,:]
    first_line2 = LineString(np.column_stack((dist, f2)))
    second_line2 = LineString(np.column_stack((dist, g2)))
    intersection2 = first_line.intersection(second_line2)

    if intersection2.geom_type == 'Point':
        x2, y2 = intersection2.xy
        zeroCrossingX2 = x2[0]

    f3 = zeroLine
    g3 = laterProfiles[hh,:]
    first_line3 = LineString(np.column_stack((dist, f3)))
    second_line3 = LineString(np.column_stack((dist, g3)))
    intersection3 = first_line.intersection(second_line3)

    if intersection3.geom_type == 'Point':
        x3, y3 = intersection3.xy
        zeroCrossingX3 = x3[0]

    differences[hh] = zeroCrossingX2-zeroCrossingX
    differences2[hh] = zeroCrossingX2-zeroCrossingX3



plt.figure()
plt.hist(differences)
plt.hist(differences2)
plt.hist(differences-differences2)

np.mean(differences-differences2)

np.mean(np.divide(differences-differences2,differences2))



surroPickle = 'surrogateOf4Storms.pickle'
output = {}
output['x'] = x
output['surrogates4'] = surrogates4
output['backXs'] = allXs

import pickle
with open(surroPickle,'wb') as f:
    pickle.dump(output, f)





plt.figure()
num =399
plt.plot(x,surrogates4[num,:],label='surrogate')
plt.plot(x,zDownscaled[num,:],label='pre')
plt.plot(x,zDownscaled2[num,:],label='post')

plt.xlim([0,350])
plt.legend()