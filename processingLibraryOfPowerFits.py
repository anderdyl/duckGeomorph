import pickle
import numpy as np
from scipy.interpolate import interp1d
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os
from shapely.geometry import LineString


geomorphdir = '/home/dylananderson/projects/duckGeomorph/scarp_data_DAversion/'
files = os.listdir(geomorphdir)
files.sort()
files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]


for i in range(len(files)):
    dbfile = open(os.path.join(files_path[i], files[i]), 'rb')
    scarpPointsTemp = pickle.load(dbfile)
    dbfile.close()
    if i == 0:
        scarpPointsMid = scarpPointsTemp[:,0,:]
        #scarpPointsLeft = scarpPointsTemp[:,0,:]
        #scarpPointsRight = scarpPointsTemp[:,2,:]
    else:
        scarpPointsMid = np.vstack((scarpPointsMid, scarpPointsTemp[:,0,:]))
        #scarpPointsLeft = np.vstack((scarpPointsLeft, scarpPointsTemp[:,0,:]))
        #scarpPointsRight = np.vstack((scarpPointsRight, scarpPointsTemp[:,2,:]))


# plt.figure()
# ax1 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
# ax1.hist((np.subtract(scarpPointsMid[:,2],scarpPointsLeft[:,2])),20)
# ax1.set_title('Middle-Left')
# ax2 = plt.subplot2grid((1,2),(0,1),rowspan=1,colspan=1)
# ax2.hist((np.subtract(scarpPointsMid[:,2],scarpPointsRight[:,2])),30)
# ax2.set_title('Middle-Right')



with open(r'case_0_499_profile_result.pickle', "rb") as input_file:
    outputProfiles = pickle.load(input_file)

preProfiles = outputProfiles['pre']
postProfiles = outputProfiles['post']
dist = outputProfiles['dist']
zeroLine = np.zeros((np.shape(postProfiles[0,0,:])))
failedSims = np.where(np.isnan(scarpPointsMid[:,0]))
zeroCrossingX = np.nan * np.ones((np.shape(scarpPointsMid[:,0])))
for i in range(len(scarpPointsMid)):
    if i in failedSims[0]:
        print('found failed sim {}.'.format(i))
    else:
        plt.figure()
        f = zeroLine
        g = preProfiles[i,1,:]
        newx = np.arange(0,100,1)
        z = scarpPointsMid[i,2] * np.power(newx, scarpPointsMid[i,3])
        first_line = LineString(np.column_stack((dist, f)))
        second_line = LineString(np.column_stack((dist, g)))
        intersection = first_line.intersection(second_line)
        plt.plot(dist, f, '-')
        plt.plot(dist,postProfiles[i,1,:],'-')
        plt.plot(dist, g, '-')
        plt.plot(scarpPointsMid[i,0],scarpPointsMid[i,1],'ko')
        plt.plot(newx+scarpPointsMid[i,0],z+scarpPointsMid[i,1],'k-')
        if intersection.geom_type == 'MultiPoint':
            plt.plot(*LineString(intersection).xy, 'o')
            x, y = LineString(intersection).xy
            #zeroCrossingX[i] = x[0]
            plt.xlim([x[0] - 100, x[0] + 100])
            plt.ylim([-3, 5])
        elif intersection.geom_type == 'Point':
            plt.plot(*intersection.xy, 'o')
            x, y =intersection.xy
            zeroCrossingX[i] = x[0]
            plt.xlim([x[0] - 100, x[0] + 100])
            plt.ylim([-3, 5])

        plt.savefig('/home/dylananderson/projects/duckGeomorph/zeroCrossingFigs/profile{}.png'.format(i))
        plt.close()



weird = np.subtract(zeroCrossingX,scarpPointsMid[:,0])
weirdProfiles = np.where((weird < 0))
plt.figure()
ax1 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
ax1.hist(scarpPointsMid[:,0],20)
ax1.set_title('Distance from Origin')
ax2 = plt.subplot2grid((1,2),(0,1),rowspan=1,colspan=1)
ax2.hist((np.subtract(zeroCrossingX,scarpPointsMid[:,0])),30)
ax2.set_title('Erosion from original shoreline')


outPower = dict()
outPower['a'] = scarpPointsMid[:,2]
outPower['b'] = scarpPointsMid[:,3]
outPower['xLoc'] = np.subtract(zeroCrossingX,scarpPointsMid[:,0])
outPower['zLoc'] = scarpPointsMid[:,1]
outPower['zeroCrossing'] = zeroCrossingX
#outPower['weirdProfiles'] = weirdProfiles
import scipy.io
scipy.io.savemat('trialsAB2.mat',outPower)


asdfg
# dbfile = open('revised_scarp_0_499.pkl', 'rb')
# scarpPoints = pickle.load(dbfile)
# dbfile.close()
#
# dbfile = open('revised_scarp_0_499_(A-B-only).pkl', 'rb')
# powerPoints = pickle.load(dbfile)
# dbfile.close()


scarpTopDist = scarpPointsMid[:,0]
scarpTopElev = scarpPointsMid[:,1]
scarpA = scarpPointsMid[:,2]
scarpB = scarpPointsMid[:,3]

profIndex = 434
fig = plt.figure()
plt.plot(dist,preProfiles[profIndex,1,:])
plt.plot(dist,postProfiles[profIndex,1,:])
plt.plot(scarpTopDist[profIndex],scarpTopElev[profIndex],'o')

#xStartPoint = np.where((dist == scarpTopDist[profIndex]))
xStartPoint = np.where((dist >= scarpTopDist[profIndex]))
newx = dist[xStartPoint[0][0]:]-scarpTopDist[profIndex]
z = scarpA[profIndex]*np.power(newx,scarpB[profIndex])
plt.plot(newx+scarpTopDist[profIndex],z+scarpTopElev[profIndex])

outPower = dict()
outPower['a'] = scarpA
outPower['b'] = scarpB
outPower['xLoc'] = scarpTopDist
outPower['zLoc'] = scarpTopElev
import scipy.io
scipy.io.savemat('trialsAB.mat',outPower)


# scarpTopDist = scarpPoints[:,1,0]
# scarpBotDist = scarpPoints[:,1,1]
# scarpTopElev = scarpPoints[:,1,2]
# scarpBotElev = scarpPoints[:,1,3]
# profBotDist = scarpPoints[:,1,4]
# profBotElev = scarpPoints[:,1,5]

asdfg

dataPred = scipy.io.loadmat('predictions2.mat')
ypreddist = dataPred['ypreddist']
ypreda = dataPred['ypreda']
ypredb = dataPred['ypredb']
valTrials = dataPred['valTrials']


fig2 = plt.figure()
c = 0
c2 = 0
for hh in range(20):
    trial = hh+0
    profIndex = valTrials[0][hh] - 1 + 0
    p5 = plt.subplot2grid((4, 5), (c, c2), rowspan=1, colspan=1)
    p5.plot(dist, preProfiles[profIndex, 1, :], color=[0, 0, 0])
    p5.plot(dist, postProfiles[profIndex, 1, :], color=[0.5, 0.5, 0.5])
    p5.plot(scarpTopDist[profIndex], scarpTopElev[profIndex], 'go')

    xStartPoint = np.where((dist >= scarpTopDist[profIndex]))
    newx = dist[xStartPoint[0][0]:] - scarpTopDist[profIndex]
    z = scarpA[profIndex] * np.power(newx, scarpB[profIndex])
    p5.plot(newx + scarpTopDist[profIndex], z + scarpTopElev[profIndex], 'g--')
    newx2 = np.arange(0, 100, 2)
    z2 = ypreda[trial] * np.power(newx2, ypredb[trial])
    zIndex = np.where((dist>(zeroCrossingX[profIndex] -ypreddist[trial]-1.5)) & (dist<(zeroCrossingX[profIndex] -ypreddist[trial]+1.5)))

    p5.plot(newx2 +zeroCrossingX[profIndex] - ypreddist[trial], z2 + preProfiles[profIndex,1,zIndex[0][0]], 'r--')
    p5.set_xlim([(zeroCrossingX[profIndex] -ypreddist[trial]-25), (zeroCrossingX[profIndex] -ypreddist[trial]+75)])
    p5.set_ylim([-3, 5])
    if c2 == 4:
        c = c+1
        c2 = 0
    else:
        c2 = c2+1


