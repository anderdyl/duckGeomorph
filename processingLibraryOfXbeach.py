import pickle
import numpy as np
from scipy.interpolate import interp1d
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt



dbfile = open('revised_scarp_0_499.pkl', 'rb')
scarpPoints = pickle.load(dbfile)
dbfile.close()


scarpTopDist = scarpPoints[:,1,0]
scarpBotDist = scarpPoints[:,1,1]
scarpTopElev = scarpPoints[:,1,2]
scarpBotElev = scarpPoints[:,1,3]
profBotDist = scarpPoints[:,1,4]
profBotElev = scarpPoints[:,1,5]
outputBreaks = dict()
outputBreaks['scarpTopDist'] = scarpTopDist
outputBreaks['scarpBotDist'] = scarpBotDist
outputBreaks['scarpTopElev'] = scarpTopElev
outputBreaks['scarpBotElev'] = scarpBotElev
outputBreaks['profBotDist'] = profBotDist
outputBreaks['profBotElev'] = profBotElev

import scipy.io
scipy.io.savemat('scarpPoints.mat',outputBreaks)

nanIndex = np.where(np.isnan(scarpTopDist))

with open(r'hypoNourishStormTides.pickle', "rb") as input_file:
    outputTrials = pickle.load(input_file)
hypoTrials = outputTrials['hypo']

outputHypos = dict()
outputHypos['hypos'] = hypoTrials
scipy.io.savemat('hypoPoints.mat',outputHypos)


with open(r'case_0_499_profile_result.pickle', "rb") as input_file:
    outputProfiles = pickle.load(input_file)
preProfiles = outputProfiles['pre']
postProfiles = outputProfiles['post']
dist = outputProfiles['dist']

outProfiles = dict()
outProfiles['pre'] = preProfiles
outProfiles['post'] = postProfiles
outProfiles['dist'] = dist
scipy.io.savemat('trialProfiles.mat',outProfiles)


asdfg

# dbfile1 = open('profiles_break_allinfo.pkl','rb')
# bothBreaks = pickle.load(dbfile1)
# dbfile1.close()



dbfile = open('profiles_break_370.pkl', 'rb')
breaksAll = pickle.load(dbfile)
dbfile.close()




# dist = breaks['dist']
# pre = breaks['pre']
# post = breaks['post']
breaks = breaksAll['break1']
middleBreakTopDist = breaks[:,1,0]
middleBreakBotDist = breaks[:,1,1]
middleBreakTopElev = breaks[:,1,2]
middleBreakBotElev = breaks[:,1,3]
middleBreakSlope = breaks[:,1,4]

#middleBreakHeight = breaks[:,1,1]

outputBreaks = dict()
outputBreaks['middleBreakTopDist'] = middleBreakTopDist
outputBreaks['middleBreakBotDist'] = middleBreakBotDist
outputBreaks['middleBreakTopElev'] = middleBreakTopElev
outputBreaks['middleBreakBotElev'] = middleBreakBotElev
outputBreaks['middleBreakSlope'] = middleBreakSlope

import scipy.io
scipy.io.savemat('profileBreaks2.mat',outputBreaks)



needZeroCrossing = np.array((8,11,35,37,47,54,60,68,73,78,86,87,88,93,98,
                            113,116,129,143,149,155,158,162,166,169,176,
                            182,196,197,198,200,202,203,210,213,214,216,217,232,235,239))


# dbfile = open('case_0_104_profile_result.pickle', 'rb')
dbfile = open('case_0_369_profile_result.pickle', 'rb')
xbProfiles = pickle.load(dbfile)
dbfile.close()
dist = xbProfiles['dist']
pre = xbProfiles['pre']
post = xbProfiles['post']



# all had the same offset of 8...
startPos = 8
rawX = dist[startPos:]#-dist[startPos]
rawX[0] = 80
rawX = rawX-80
rawZpreN = pre[:,0,startPos:]
rawZpreM = pre[:,1,startPos:]
rawZpreS = pre[:,2,startPos:]
rawZpostN = post[:,0,startPos:]
rawZpostM = post[:,1,startPos:]
rawZpostS = post[:,2,startPos:]

zeroCrossing = np.empty((np.size(needZeroCrossing)))
xcounter = 0
for xxx in needZeroCrossing:
    tempProfile = np.abs(rawZpostM[xxx]-1)
    findMinimum = np.argmin(tempProfile)
    zeroCrossing[xcounter] = rawX[findMinimum]
    xcounter = xcounter+1


#
#
# newX = np.arange(0,1199,1)

xInitial = np.arange(0,30,5)
xIM = np.arange(30,180,2)
xMiddle = np.arange(180,440,5)
xME = np.arange(440,840,10)
x = np.hstack((xInitial,xIM))
x = np.hstack((x,xMiddle))
x = np.hstack((x,xME))

allZsPreN = []
allZsPreM = []
allZsPreS = []
allZsPostN = []
allZsPostM = []
allZsPostS = []

allNums = []

allXsPreN = []
allXsPreM = []
allXsPreS = []
allXsPostN = []
allXsPostM = []
allXsPostS = []
for xx in range(len(rawZpreN)):
    print(xx)
    # need a section in here to skip the scenarios that crashed 44,63,121,161,209
    if xx == 44 or xx == 63 or xx == 121 or xx == 161 or xx == 209 or xx == 308 or xx == 359:
        print('skipped {}'.format(xx))
    else:
    # all had the same offset of 8...
    # so we just need to interpolate on to the old EOF space?
        # starting with the North Profile
        tempZpreN = rawZpreN[xx, :]
        f2preN = interp1d(rawX[0:len(tempZpreN)], tempZpreN, kind='linear')
        zDownscaledPreN = f2preN(x)
        allZsPreN.append(zDownscaledPreN)
        # now for the pre Middle Profile
        tempZpreM = rawZpreM[xx, :]
        f2preM = interp1d(rawX[0:len(tempZpreM)], tempZpreM, kind='linear')
        zDownscaledPreM = f2preM(x)
        allZsPreM.append(zDownscaledPreM)
        # now for the pre South Profile
        tempZpreS = rawZpreS[xx, :]
        f2preS = interp1d(rawX[0:len(tempZpreS)], tempZpreS, kind='linear')
        zDownscaledPreS = f2preS(x)
        allZsPreS.append(zDownscaledPreS)

        # now for the post South Profile
        tempZpostS = rawZpostS[xx, :]
        f2postS = interp1d(rawX[0:len(tempZpostS)], tempZpostS, kind='linear')
        zDownscaledPostS = f2postS(x)
        allZsPostS.append(zDownscaledPostS)
        # now for the post Middle Profile
        tempZpostM = rawZpostM[xx, :]
        f2postM = interp1d(rawX[0:len(tempZpostM)], tempZpostM, kind='linear')
        zDownscaledPostM = f2postM(x)
        allZsPostM.append(zDownscaledPostM)
        # now for the post North Profile
        tempZpostN = rawZpostN[xx, :]
        f2postN = interp1d(rawX[0:len(tempZpostN)], tempZpostN, kind='linear')
        zDownscaledPostN = f2postN(x)
        allZsPostN.append(zDownscaledPostN)



allZsPostMarray = np.asarray(allZsPostM)
allZsPostNarray = np.asarray(allZsPostN)
allZsPostSarray = np.asarray(allZsPostS)
allZsPreMarray = np.asarray(allZsPreM)
allZsPreSarray = np.asarray(allZsPreS)
allZsPreNarray = np.asarray(allZsPreN)

dbfile = 'nourishmentProfileEOFs.pickle'

with open(r'nourishmentProfileEOFs.pickle', "rb") as input_file:
    outputEOFs = pickle.load(input_file)

OGmeanZ = outputEOFs['meanZ']
OGEOFs = outputEOFs['EOFs']
OGPCs = outputEOFs['PCs']
OGvariance = outputEOFs['variance']
OGvarianceRatio = outputEOFs['varianceRatio']
OGcumulativeVar = outputEOFs['cumulativeVar']
OGzInt = outputEOFs['zInt']
OGx = outputEOFs['x']



with open(r'hypoNourishStormTides.pickle', "rb") as input_file:
    outputTrials = pickle.load(input_file)

hypoTrials = outputTrials['hypo']

for xx in range(len(rawZpreN)):
    print(xx)

    # need a section in here to skip the scenarios that crashed
    if xx == 44 or xx == 63 or xx == 121 or xx == 161 or xx == 209 or xx == 308 or xx == 359:
        print('skipped {}'.format(xx))
    else:
        if xx == 0:
            runTrials = hypoTrials[xx, :]
        else:
            runTrials = np.vstack((runTrials,hypoTrials[xx,:]))

originalIPCA = PCA(n_components=7)

# PCs = ipca.fit_transform(normZ)
origPCs = originalIPCA.fit_transform(OGzInt)

# inputPCs = originalIPCA.transform()
shouldBeOrigPCs = originalIPCA.transform(allZsPreMarray)

# newPCs = originalIPCA.transform(zDownscaled)
# newProfiles = originalIPCA.inverse_transform(newPCs)

fig = plt.figure()
ax1 = plt.subplot2grid((1,4),(0,0),rowspan=1,colspan=1)
ax1.hist(np.subtract(runTrials[:,0],shouldBeOrigPCs[:,0]),10)
ax1.set_title('EOF1')
ax2 = plt.subplot2grid((1,4),(0,1),rowspan=1,colspan=1)
ax2.hist(np.subtract(runTrials[:,1],shouldBeOrigPCs[:,1]),10)
ax2.set_title('EOF2')
ax3 = plt.subplot2grid((1,4),(0,2),rowspan=1,colspan=1)
ax3.hist(np.subtract(runTrials[:,2],shouldBeOrigPCs[:,2]),10)
ax3.set_title('EOF3')
ax4 = plt.subplot2grid((1,4),(0,3),rowspan=1,colspan=1)
ax4.hist(np.subtract(runTrials[:,3],shouldBeOrigPCs[:,3]),10)
ax4.set_title('EOF4')


shouldBeOrigPCsM = originalIPCA.transform(allZsPreMarray)
shouldBeOrigPCsS = originalIPCA.transform(allZsPreSarray)
shouldBeOrigPCsN = originalIPCA.transform(allZsPreNarray)
postPCsM = originalIPCA.transform(allZsPostMarray)
postPCsS = originalIPCA.transform(allZsPostSarray)
postPCsN = originalIPCA.transform(allZsPostNarray)

diffpreNS = np.divide(np.subtract(shouldBeOrigPCsN[:,0],shouldBeOrigPCsS[:,0]),2)
diffpostNS = np.divide(np.subtract(postPCsN[:,0],postPCsS[:,0]),2)


# plt.figure()
# plt.hist(np.subtract(runTrials[:,7],diffpreNS))

newEOFscores = np.empty((363,8))
newEOFscores[:,0:7] = postPCsM
newEOFscores[:,7] = diffpostNS
# newEOFscores = np.hstack((postPCsM,diffpostNS))
# newEOFscores = np.hstack((newEOFscores1,runTrials))



predictands = newEOFscores

output = dict()
output['predictors'] = runTrials
output['predictands'] = predictands
# import scipy.io
scipy.io.savemat('profileTrainingLibrary2Updated3.mat',output)




####

dataPred = scipy.io.loadmat('predictions')
ypred = dataPred['ypred']
ypred2 = dataPred['ypred2']
ypred3 = dataPred['ypred3']
ypred4 = dataPred['ypred4']
ypred5 = dataPred['ypred5']
ypred6 = dataPred['ypred6']
ypred7 = dataPred['ypred7']
ypreddist = dataPred['ypreddist']
ypredheight = dataPred['ypredheight']
valTrials = dataPred['valTrials']

predictions = np.empty((24,7))
predictions[:,0] = ypred.flatten()
predictions[:,1] = ypred2.flatten()
predictions[:,2] = ypred3.flatten()
predictions[:,3] = ypred4.flatten()
predictions[:,4] = ypred5.flatten()
predictions[:,5] = ypred6.flatten()
predictions[:,6] = ypred7.flatten()

predictions = dataPred['final']



pyVals = valTrials-1



plt.figure(figsize=(16,8))
calProfiles = []
c = 0
c2 = 0
for xx in range(len(pyVals[0])):

    #prof5 = OGmeanZ + OGEOFs[0,:]*predictions[xx,0] + OGEOFs[1,:]*predictions[xx,1] + OGEOFs[2,:]*predictions[xx,2] \
    #        + OGEOFs[3,:]*predictions[xx,3] + OGEOFs[4,:]*predictions[xx,4] + OGEOFs[5,:]*predictions[xx,5] \
    #        + OGEOFs[6, :] * predictions[xx, 6]

    scarp = np.zeros((len(x)))
    len1 = ypredheight[xx]/.25
    len2 = ypredheight[xx]/.05

    findScarpStart = np.argmin(np.abs(x-ypreddist[xx]))
    findScarpEnd = np.argmin(np.abs(x-len1-ypreddist[xx]))
    findReturnEnd = np.argmin(np.abs(x-len1-len2-ypreddist[xx]))
    b1 = .25*x[findScarpStart]
    y1 = -.25*x+b1
    b2 = -.05*x[findReturnEnd]
    y2 = .05*x+b2
    scarp[findScarpStart:findScarpEnd] = y1[findScarpStart:findScarpEnd]
    scarp[findScarpEnd:findReturnEnd] = y2[findScarpEnd:findReturnEnd]
    # prof5 = meanZ + EOFs[0,:]*predictions1[xx] + EOFs[1,:]*predictions2[xx] + EOFs[2,:]*predictions3[xx] \
    #         + EOFs[3,:]*predictions4[xx] + EOFs[4,:]*predictions5[xx]
    justScarp = allZsPreMarray[pyVals[0][xx],:]+scarp
    #calProfiles.append(prof5)
    ax = plt.subplot2grid((5,5),(c2,c),rowspan=1,colspan=1)
    ax.plot(x,allZsPreMarray[pyVals[0][xx],:],label='Pre-storm ')
    ax.plot(x,allZsPostMarray[pyVals[0][xx],:],label='Post-storm Xbeach')
    ax.plot(x,justScarp,label='Just Scarp on Pre')
    #ax.plot(x,prof5,label='Post-storm GPR')
    ax.set_xlim([0, 300])
    ax.set_ylim([-3.5,7])
    ax.set_title('Trial {}'.format(pyVals[0][xx]))
    if c == 4:
        c = 0
        c2 = c2+1
    else:
        c = c + 1
ax.legend()