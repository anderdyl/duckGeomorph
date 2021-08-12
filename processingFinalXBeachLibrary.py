
import pickle
import numpy as np
from scipy.interpolate import interp1d
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os
from shapely.geometry import LineString


#
# geomorphdir = '/home/dylananderson/projects/duckGeomorph/scarp_data_final2/'
# files = os.listdir(geomorphdir)
# files.sort()
# files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]
#
#
#
# for i in range(10):
#     dbfile = open(os.path.join(files_path[i+7], files[i+7]), 'rb')
#     scarpPointsTemp = pickle.load(dbfile)
#     dbfile.close()
#     if i == 0:
#         scarpPointsMid = scarpPointsTemp[:,:]
#     else:
#         scarpPointsMid = np.vstack((scarpPointsMid, scarpPointsTemp[:,:]))
#
# scarpPointsMid = scarpPointsMid[:,0,:]
# for i in range(5):
#     dbfile = open(os.path.join(files_path[i+2], files[i+2]), 'rb')
#     scarpPointsTemp = pickle.load(dbfile)
#     dbfile.close()
#     scarpPointsMid = np.vstack((scarpPointsMid, scarpPointsTemp[:,:]))
#

geomorphdir = '/home/dylananderson/projects/duckGeomorph/scarp_data_full_DAversion/'
files = os.listdir(geomorphdir)
files.sort()
files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]



for i in range(25):
    dbfile = open(os.path.join(files_path[i+3], files[i+3]), 'rb')
    scarpPointsTemp = pickle.load(dbfile)
    dbfile.close()
    if i == 0:
        scarpPointsMid = scarpPointsTemp[:,:]
    else:
        scarpPointsMid = np.vstack((scarpPointsMid, scarpPointsTemp[:,:]))


#scarpPointsMid = scarpPointsMid[:,0,:]


with open(r'/home/dylananderson/projects/duckGeomorph/scarp_data_full_DAversion/newcase_profiles_0000_0500_revised.pkl', "rb") as input_file:
    outputProfiles = pickle.load(input_file)

with open(r'/home/dylananderson/projects/duckGeomorph/scarp_data_full_DAversion/newcase_profiles_0500_1000.pkl', "rb") as input_file:
    outputProfiles2 = pickle.load(input_file)

with open(r'/home/dylananderson/projects/duckGeomorph/scarp_data_full_DAversion/newcase_profiles_1000_1250.pkl', "rb") as input_file:
    outputProfiles3 = pickle.load(input_file)

preProfiles1 = np.vstack((outputProfiles['pre'],outputProfiles2['pre']))
postProfiles1 = np.vstack((outputProfiles['post'],outputProfiles2['post']))
preProfiles = np.vstack((preProfiles1,outputProfiles3['pre']))
postProfiles = np.vstack((postProfiles1,outputProfiles3['post']))


dist = outputProfiles['dist']
zeroLine = np.zeros((np.shape(postProfiles[0,:])))
failedSims = np.where(np.isnan(scarpPointsMid[:,0]))
zeroCrossingX = np.nan * np.ones((np.shape(scarpPointsMid[:,0])))
for i in range(len(scarpPointsMid)):
    if i in failedSims[0]:
        print('found failed sim {}.'.format(i))
    else:
        plt.figure()
        f = zeroLine
        g = preProfiles[i,:]
        newx = np.arange(0,100,1)
        z = scarpPointsMid[i,2] * np.power(newx, scarpPointsMid[i,3])
        first_line = LineString(np.column_stack((dist, f)))
        second_line = LineString(np.column_stack((dist, g)))
        intersection = first_line.intersection(second_line)
        plt.plot(dist, f, '-')
        plt.plot(dist,postProfiles[i,:],'-')
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

        plt.savefig('/home/dylananderson/projects/duckGeomorph/zeroCrossingFigsFinal/profile{}.png'.format(i))
        plt.close()



outPower = dict()
outPower['a'] = scarpPointsMid[:,2]
outPower['b'] = scarpPointsMid[:,3]
outPower['xLoc'] = np.subtract(zeroCrossingX,scarpPointsMid[:,0])
outPower['zLoc'] = scarpPointsMid[:,1]
outPower['zeroCrossing'] = zeroCrossingX
#outPower['weirdProfiles'] = weirdProfiles
import scipy.io
scipy.io.savemat('trialsABFinal1250.mat',outPower)

scarpTopDist = scarpPointsMid[:,0]
scarpTopElev = scarpPointsMid[:,1]
scarpBotDist = scarpPointsMid[:,4]
scarpBotElev = scarpPointsMid[:,5]
scarpA = scarpPointsMid[:,2]
scarpB = scarpPointsMid[:,3]

outScarpPoints = dict()
outScarpPoints['scarpTopDist'] = scarpPointsMid[:,0]
outScarpPoints['scarpTopElev'] = scarpPointsMid[:,1]
outScarpPoints['scarpBotDist'] = scarpPointsMid[:,4]
outScarpPoints['scarpBotElev'] = scarpPointsMid[:,5]
import scipy.io
scipy.io.savemat('scarpPointsFinal1250.mat',outScarpPoints)



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
    trial = hh+40
    profIndex = valTrials[0][trial] - 1 + 0
    p5 = plt.subplot2grid((4, 5), (c, c2), rowspan=1, colspan=1)
    p5.plot(dist, preProfiles[profIndex, :], color=[0, 0, 0])
    p5.plot(dist, postProfiles[profIndex, :], color=[0.5, 0.5, 0.5])
    p5.plot(scarpTopDist[profIndex], scarpTopElev[profIndex], 'go')

    xStartPoint = np.where((dist >= scarpTopDist[profIndex]))
    newx = dist[xStartPoint[0][0]:] - scarpTopDist[profIndex]
    z = scarpA[profIndex] * np.power(newx, scarpB[profIndex])
    p5.plot(newx + scarpTopDist[profIndex], z + scarpTopElev[profIndex], 'g--')
    newx2 = np.arange(0, 100, 2)
    z2 = ypreda[trial] * np.power(newx2, ypredb[trial])
    zIndex = np.where((dist>(zeroCrossingX[profIndex] -ypreddist[trial]-1.5)) & (dist<(zeroCrossingX[profIndex] -ypreddist[trial]+1.5)))

    p5.plot(newx2 +zeroCrossingX[profIndex] - ypreddist[trial], z2 + preProfiles[profIndex,zIndex[0][0]], 'r--')
    p5.set_xlim([(zeroCrossingX[profIndex] -ypreddist[trial]-25), (zeroCrossingX[profIndex] -ypreddist[trial]+75)])
    p5.set_ylim([-3, 5])
    if c2 == 4:
        c = c+1
        c2 = 0
    else:
        c2 = c2+1
