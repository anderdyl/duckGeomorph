import pickle
import numpy as np
from scipy.interpolate import interp1d
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os
from shapely.geometry import LineString
import scipy.io
import scipy.integrate as integrate
from shapely.geometry import Polygon


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

scarpPointsMid = scarpPointsMid[:,0,:]



import pickle
fitPickle = 'scarpFits.pickle'
scarpFits = {}
scarpFits['fits'] = scarpPointsMid
with open(fitPickle,'wb') as f:
    pickle.dump(scarpFits, f)


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


xBeachPred = scipy.io.loadmat('trialsABFinal1250.mat')
aParam = xBeachPred['a']
bParam = xBeachPred['b']
xLocParam = xBeachPred['xLoc']
zLocParam = xBeachPred['zLoc']
zeroCrossing = xBeachPred['zeroCrossing']

scarpTopDist = scarpPointsMid[:,0]
scarpTopElev = scarpPointsMid[:,1]
scarpBotDist = scarpPointsMid[:,4]
scarpBotElev = scarpPointsMid[:,5]
scarpA = scarpPointsMid[:,2]
scarpB = scarpPointsMid[:,3]

dataPred = scipy.io.loadmat('predictions3.mat')
ypreddist = dataPred['ypreddist']
ypreda = dataPred['ypreda']
ypredb = dataPred['ypredb']
valTrials = dataPred['valTrials']


fig2 = plt.figure()
c = 0
c2 = 0
actualArea = list()
predictedArea = list()


chosen = [26,4,66,35,58,49,80,19,94]

textScenario = ['i.','ii.','iii.','iv.','v.','vi.','vii.','viii.','ix.']

for hh in range(9):
    trial = chosen[hh]
    #trial = hh+0
    profIndex = valTrials[0][trial] - 1 + 0
    p5 = plt.subplot2grid((3, 3), (c, c2), rowspan=1, colspan=1)

    p5.plot(dist, preProfiles[profIndex, :], color=[0, 0, 0],label='Pre-Storm')
    p5.plot(dist, postProfiles[profIndex, :], color=[0.5, 0.5, 0.5],label='Post-Storm XBeach')

    #p5.plot(scarpTopDist[profIndex], scarpTopElev[profIndex], 'go')

    xStartPoint = np.where((dist >= scarpTopDist[profIndex]))
    newx = dist[xStartPoint[0][0]:] - scarpTopDist[profIndex]
    z = scarpA[profIndex] * np.power(newx, scarpB[profIndex])
    #p5.plot(newx + scarpTopDist[profIndex], z + scarpTopElev[profIndex], 'g--')
    newx2 = np.arange(0, 100, 1)
    z2 = ypreda[trial] * np.power(newx2, ypredb[trial])
    zIndex = np.where((dist>(zeroCrossing[0][profIndex] -ypreddist[trial]-1.5)) & (dist<(zeroCrossing[0][profIndex] -ypreddist[trial]+1.5)))

    p5.plot(newx2 +zeroCrossing[0][profIndex] - ypreddist[trial], z2 + preProfiles[profIndex,zIndex[0][0]], 'r--',label='Post-Storm Surrogate')
    p5.set_xlim([(zeroCrossing[0][profIndex] -ypreddist[trial]-25), (zeroCrossing[0][profIndex] -ypreddist[trial]+75)])
    p5.set_ylim([-3, 5])

    # Lets make a new line
    ogDuneInd = np.where((dist < (zeroCrossing[0][profIndex] - ypreddist[trial])))
    scarpXDune = dist[ogDuneInd]
    scarpZDune = preProfiles[profIndex,ogDuneInd]

    scarpXBeach = newx2 +zeroCrossing[0][profIndex] - ypreddist[trial]
    scarpZBeach = z2 + preProfiles[profIndex,zIndex[0][0]]

    ogSubaqueous = np.where((dist > scarpXBeach[-1]))
    scarpXSubaqueous = dist[ogSubaqueous]
    scarpZSubaqueous = preProfiles[profIndex,ogSubaqueous]

    predScarp1X = np.hstack((scarpXDune,scarpXBeach))
    predScarp2X = np.hstack((predScarp1X,scarpXSubaqueous))
    predScarp1Z = np.hstack((scarpZDune[0], scarpZBeach))
    predScarp2Z = np.hstack((predScarp1Z, scarpZSubaqueous[0]))

    findNegatives = np.where((predScarp2Z < 0))
    predScarp2Z[findNegatives] = predScarp2Z[findNegatives]*0
    #p5.plot(predScarp2X,predScarp2Z,'m--')

    # Difference between the XBeach volumes
    x_y_curveTest1 = list()
    x_y_curveTest2 = list()
    z1 = preProfiles[profIndex, :]
    z2 = postProfiles[profIndex, :]
    findNegatives2 = np.where((z2<0))
    z2[findNegatives2] = z2[findNegatives2]*0
    finder = np.where((z1<(z2-.15)))

    for n in range(finder[0][0]):
        x_y_curveTest1.append(list((dist[n], z1[n])))
        x_y_curveTest2.append(list((dist[n], z2[n])))

    polygon_points = []  # creates a empty list where we will append the points to create the polygon

    for xyvalue in x_y_curveTest1:
        polygon_points.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 1

    for xyvalue in x_y_curveTest2[::-1]:
        polygon_points.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 2 in the reverse order (from last point to first point)

    for xyvalue in x_y_curveTest1[0:1]:
        polygon_points.append([xyvalue[0], xyvalue[1]])  # append the first point in curve 1 again, to it "closes" the polygon

    polygon = Polygon(polygon_points)
    area = polygon.area
    actualArea.append(area)
    #(area)
    #p5.text(zeroCrossing[0][profIndex]-25,3.5,'{} m^2'.format(np.round(area*10)/10))
    #fdu = integrate.simps()


    # Difference between the XBeach volumes
    x_y_curveTest12 = list()
    x_y_curveTest22 = list()
    z12 = preProfiles[profIndex, :]
    f = interp1d(dist,z12,kind='cubic')
    z12new = f(np.ma.filled(predScarp2X))
    z22 = predScarp2Z
    finder = np.where((z12new<(z22-.1)))
    if finder[0][0]<(finder[0][1]-5):
        for n in range(finder[0][1]):
            x_y_curveTest12.append(list((predScarp2X[n], z12new[n])))
            x_y_curveTest22.append(list((predScarp2X[n], z22[n])))
    else:
        for n in range(finder[0][0]):
            x_y_curveTest12.append(list((predScarp2X[n], z12new[n])))
            x_y_curveTest22.append(list((predScarp2X[n], z22[n])))

    polygon_points2 = []  # creates a empty list where we will append the points to create the polygon

    for xyvalue in x_y_curveTest12:
        polygon_points2.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 1

    for xyvalue in x_y_curveTest22[::-1]:
        polygon_points2.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 2 in the reverse order (from last point to first point)

    for xyvalue in x_y_curveTest12[0:1]:
        polygon_points2.append([xyvalue[0], xyvalue[1]])  # append the first point in curve 1 again, to it "closes" the polygon

    polygon2 = Polygon(polygon_points2)
    area2 = polygon2.area
    predictedArea.append(area2)

    #(area)
    # p5.text(zeroCrossing[0][profIndex]-25,2,'{} m^2'.format(np.round(area2*10)/10),color='m')
    #p5.text()
    p5.text(0.05, 0.88, textScenario[hh], transform=p5.transAxes,
            size=14, weight='bold')
    if c2 == 0:
        p5.set_ylabel('elevation (m)')
    if c == 2:
        p5.set_xlabel('cross-shore (m)')
    if hh == 0:
        p5.legend()
    if c2 == 2:
        c = c+1
        c2 = 0
    else:
        c2 = c2+1


def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

predArray = np.asarray(predictedArea)
actuArray = np.asarray(actualArea)
rootMeanSquared = rmse(predArray,actuArray)
errors = np.abs(predArray-actuArray)
medianError = np.median(errors)
plt.figure()
plt.hist(np.asarray(predictedArea)-np.asarray(actualArea))
plt.title('RMSE = {} m^2, Median Absolute Error = {} m^2'.format(np.round(rootMeanSquared*100)/100,np.round(medianError*100)/100))
