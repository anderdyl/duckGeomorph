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

from netCDF4 import Dataset
import datetime as DT




waterLevelDir = '/media/dylananderson/Elements/frfWaterLevel'
files = os.listdir(waterLevelDir)
files.sort()
#files = files[0:228]
files_path = [os.path.join(os.path.abspath(waterLevelDir), x) for x in files]
wls = Dataset(files_path[0])

def getWaterLevel(file):
    wldata = Dataset(file)
    waterLevel = wldata.variables['waterLevel'][:]
    predictedWaterLevel = wldata.variables['predictedWaterLevel'][:]
    residualWaterLevel = wldata.variables['residualWaterLevel'][:]
    timeWl = wldata.variables['time'][:]
    output = dict()
    output['waterLevel'] = waterLevel
    output['predictedWaterLevel'] = predictedWaterLevel
    output['residualWaterLevel'] = residualWaterLevel
    output['time'] = timeWl
    return output

timeWaterLevelFRF = []
waterLevelFRF = []
predictedWaterLevelFRF = []
residualWaterLevelFRF = []
for i in files_path:
    print('getting file {}'.format(i))
    waterLevels = getWaterLevel(i)
    waterLevelFRF = np.append(waterLevelFRF,waterLevels['waterLevel'])
    predictedWaterLevelFRF = np.append(predictedWaterLevelFRF,waterLevels['predictedWaterLevel'])
    residualWaterLevelFRF = np.append(residualWaterLevelFRF,waterLevels['residualWaterLevel'])
    timeWaterLevelFRF = np.append(timeWaterLevelFRF,waterLevels['time'].flatten())


badWaterLevel = np.where((residualWaterLevelFRF < -99))
tWaterLevelFRF = [DT.datetime.fromtimestamp(x) for x in timeWaterLevelFRF]
tWaterLevelFRF = np.asarray(tWaterLevelFRF)
waterLevelFRF[badWaterLevel] = waterLevelFRF[badWaterLevel]*np.nan
# predictedWaterLevelFRF[badWaterLevel] = predictedWaterLevelFRF[badWaterLevel]*np.nan
# residualWaterLevelFRF[badWaterLevel] = residualWaterLevelFRF[badWaterLevel]*np.nan
#
# residualWaterLevelFRF = residualWaterLevelFRF-np.nanmean(residualWaterLevelFRF)
#

# timestampSeconds = [ff.timestamp() for ff in outputTide1[0]]
# timestampHours = np.array([ff/(60*60) for ff in timestampSeconds])
timeStampHours = timeWaterLevelFRF/(60*60)
plt.figure()
plt.plot(tWaterLevelFRF,waterLevelFRF)


outTides = dict()
outTides['predTide'] = waterLevelFRF
# outTides['datetime'] = tWaterLevelFRF
outTides['timeStampHours'] = timeStampHours
import scipy.io
scipy.io.savemat('frfTideFitting.mat',outTides)





with open(r'/home/dylananderson/projects/duckGeomorph/s4_tide.pkl', "rb") as input_file:
    outputTide1 = pickle.load(input_file)


#matlabTime = np.array([ff-DT.datetime(1,1,1)+366 for ff in outputTide1[0]])

timestampSeconds = [ff.timestamp() for ff in outputTide1[0]]
timestampHours = np.array([ff/(60*60) for ff in timestampSeconds])
timestampDays = np.array([ff/(24) for ff in timestampHours])
matlabTime = timestampDays+719529#-4/24

# m2 = 0.49 * np.cos(matlabTime*2*np.pi*1.9323 + - 161.1476*np.pi/180)
# s2 = 0.088 * np.cos(matlabTime*2*np.pi*2 + - 21.9045*np.pi/180)

m2 = 0.49 * np.cos(matlabTime*2*np.pi*1.932273615 - 159.4021*np.pi/180)
s2 = 0.088 * np.cos(matlabTime*2*np.pi*2 - 22.3066*np.pi/180)
n2 = 0.114 * np.cos(matlabTime*2*np.pi*1.8959819 - 73.214*np.pi/180)
k1 = 0.087 * np.cos(matlabTime*2*np.pi*1.0027379 - 168.0082*np.pi/180)
totalTide = m2+s2+n2+k1

plt.figure()
p1 = plt.subplot2grid((2,1),(0,0),rowspan=1,colspan=1)
p1.plot(outputTide1[0],outputTide1[1])
p1.plot(outputTide1[0],totalTide)

p2 = plt.subplot2grid((2,1),(1,0),rowspan=1,colspan=1)
p2.plot(outputTide1[0],m2,label='m2')
p2.plot(outputTide1[0],s2,label='s2')
p2.plot(outputTide1[0],n2,label='n2')
p2.plot(outputTide1[0],k1,label='k1')



