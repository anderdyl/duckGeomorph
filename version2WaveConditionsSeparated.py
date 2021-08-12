
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


dbfile = open('sandbarsSouthernTransect_referencedMHHW_reallylongWithNegatives.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()
alllines = data['alllines']
xinterp = data['xinterp']
time = data['time']


dt = time[1:]-time[0:-1]

days = [t.days for t in dt]#[item for sublist in m for item in sublist]
days2 = np.array(days)
# plt.figure()
# plt.plot(days)
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


def conditional_transition_matrix(transitions,num_of_states):
    n = num_of_states #1+ max(transitions) #number of states
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


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


#
# mat = scipy.io.loadmat('matlabsSouthernTransectClustersComplex.mat')
#
# numCenters = int(mat['clusters']['NumberCenters'][0][0].flatten())
# #group = mat['clusters']['Group'][0][0]
# bmus = mat['clusters']['PatternsGroup'][0][0].flatten()
# order = mat['order'].flatten()
# #sorted_bmus = np.zeros((len(bmus),))

clusterPickle = 'sandbarsSouthernTransectPythonClusters8.pickle'
dbfile = open(clusterPickle, 'rb')
dataClusters = pickle.load(dbfile)
dbfile.close()
numCenters = dataClusters['numCenters']
bmus = dataClusters['bmu']
order = dataClusters['order']
profiles = dataClusters['profiles']


sorted_bmus = np.tile(0,(len(bmus),), )

#sorted_time = np.tile(0,(len(kma.labels_),), )
numClusters = numCenters

for i in range(numCenters):
    posc = np.where(bmus == order[i])
    sorted_bmus[posc] = int(i)
    #sorted_time[posc] = time[posc]

# def transitionDictionary(bmus,dates):

bins = dict()
date = dict()
nextbin = dict()
nextdate = dict()
prevbin = dict()
prevdate = dict()
daysBetween = dict()
for xx in range(numCenters):
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
    daysBetween[xx] = days2[res2]






windsDir = '/media/dylananderson/Elements/frfWinds'
files = os.listdir(windsDir)
files.sort()
files_path = [os.path.join(os.path.abspath(windsDir), x) for x in files]
winds = Dataset(files_path[0])

def getWinds(file):
    winddata = Dataset(file)
    windSpeed = winddata.variables['windSpeed'][:]
    vectorSpeed = winddata.variables['vectorSpeed'][:]
    sustWindSpeed = winddata.variables['sustWindSpeed'][:]
    windGust = winddata.variables['windGust'][:]
    windDirection = winddata.variables['windDirection'][:]
    maxWindSpeed = winddata.variables['maxWindSpeed'][:]
    minWindSpeed = winddata.variables['minWindSpeed'][:]
    stdWindSpeed = winddata.variables['stdWindSpeed'][:]
    qcFlagS = winddata.variables['qcFlagS'][:]
    time = winddata.variables['time'][:]
    output = dict()
    output['windSpeed'] = windSpeed
    output['vectorSpeed'] = vectorSpeed
    output['sustWindSpeed'] = sustWindSpeed
    output['windGust'] = windGust
    output['windDirection'] = windDirection
    output['maxWindSpeed'] = maxWindSpeed
    output['minWindSpeed'] = minWindSpeed
    output['stdWindSpeed'] = stdWindSpeed
    output['qcFlagS'] = qcFlagS
    output['time'] = time

    return output


windSpeed = []
vectorSpeed = []
sustWindSpeed = []
windGust = []
windDirection = []
maxWindSpeed = []
minWindSpeed = []
stdWindSpeed = []
qcFlagS = []
timeWind = []

for i in files_path:
    windsFRF = getWinds(i)
    windSpeed = np.append(windSpeed,windsFRF['windSpeed'])
    vectorSpeed = np.append(vectorSpeed,windsFRF['sustWindSpeed'])
    sustWindSpeed = np.append(sustWindSpeed,windsFRF['sustWindSpeed'])
    windGust = np.append(windGust,windsFRF['windGust'])
    windDirection = np.append(windDirection,windsFRF['windDirection'])
    maxWindSpeed = np.append(maxWindSpeed,windsFRF['maxWindSpeed'])
    minWindSpeed = np.append(minWindSpeed,windsFRF['minWindSpeed'])
    qcFlagS = np.append(qcFlagS,windsFRF['qcFlagS'])
    timeWind = np.append(timeWind,windsFRF['time'])

tWinds = [DT.datetime.fromtimestamp(x) for x in timeWind]
tWinds = np.asarray(tWinds)
badWinds = np.where((sustWindSpeed < -1))
sustWindSpeed[badWinds] = sustWindSpeed[badWinds]*np.nan
windDirection[badWinds] = windDirection[badWinds]*np.nan
windSpeed[badWinds] = windSpeed[badWinds]*np.nan
windGust[badWinds] = windGust[badWinds]*np.nan

badWinds = np.where((qcFlagS == 3))
sustWindSpeed[badWinds] = sustWindSpeed[badWinds]*np.nan
windDirection[badWinds] = windDirection[badWinds]*np.nan
windSpeed[badWinds] = windSpeed[badWinds]*np.nan
windGust[badWinds] = windGust[badWinds]*np.nan

badWinds = np.where((qcFlagS == 2))
sustWindSpeed[badWinds] = sustWindSpeed[badWinds]*np.nan
windDirection[badWinds] = windDirection[badWinds]*np.nan
windSpeed[badWinds] = windSpeed[badWinds]*np.nan
windGust[badWinds] = windGust[badWinds]*np.nan


#
#
# tWaterLevelFRF = [DT.datetime.fromtimestamp(x) for x in timeWaterLevelFRF]
#
# waterLevelFRF[badWaterLevel] = waterLevelFRF[badWaterLevel]*np.nan
# predictedWaterLevelFRF[badWaterLevel] = predictedWaterLevelFRF[badWaterLevel]*np.nan
# residualWaterLevelFRF[badWaterLevel] = residualWaterLevelFRF[badWaterLevel]*np.nan



waterLevelDir = '/media/dylananderson/Elements/frfWaterLevel'
files = os.listdir(waterLevelDir)
files.sort()
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
    waterLevels = getWaterLevel(i)
    waterLevelFRF = np.append(waterLevelFRF,waterLevels['waterLevel'])
    predictedWaterLevelFRF = np.append(predictedWaterLevelFRF,waterLevels['predictedWaterLevel'])
    residualWaterLevelFRF = np.append(residualWaterLevelFRF,waterLevels['residualWaterLevel'])
    timeWaterLevelFRF = np.append(timeWaterLevelFRF,waterLevels['time'].flatten())


badWaterLevel = np.where((residualWaterLevelFRF < -99))
tWaterLevelFRF = [DT.datetime.fromtimestamp(x) for x in timeWaterLevelFRF]
tWaterLevelFRF = np.asarray(tWaterLevelFRF)
waterLevelFRF[badWaterLevel] = waterLevelFRF[badWaterLevel]*np.nan
predictedWaterLevelFRF[badWaterLevel] = predictedWaterLevelFRF[badWaterLevel]*np.nan
residualWaterLevelFRF[badWaterLevel] = residualWaterLevelFRF[badWaterLevel]*np.nan


wavedir = '/media/dylananderson/Elements1/WIS_ST63218/'

# Need to sort the files to ensure correct temporal order...
files = os.listdir(wavedir)
files.sort()
files_path = [os.path.join(os.path.abspath(wavedir), x) for x in files]

wis = Dataset(files_path[0])

def getWIS(file):
    waves = Dataset(file)

    waveHs = waves.variables['waveHs'][:]
    waveTp = waves.variables['waveTp'][:]
    waveMeanDirection = waves.variables['waveMeanDirection'][:]

    waveTm = waves.variables['waveTm'][:]
    waveTm1 = waves.variables['waveTm1'][:]
    waveTm2 = waves.variables['waveTm2'][:]

    waveHsWindsea = waves.variables['waveHsWindsea'][:]
    waveTmWindsea = waves.variables['waveTmWindsea'][:]
    waveMeanDirectionWindsea = waves.variables['waveMeanDirectionWindsea'][:]
    waveSpreadWindsea = waves.variables['waveSpreadWindsea'][:]

    timeW = waves.variables['time'][:]

    waveTpSwell = waves.variables['waveTpSwell'][:]
    waveHsSwell = waves.variables['waveHsSwell'][:]
    waveMeanDirectionSwell = waves.variables['waveMeanDirectionSwell'][:]
    waveSpreadSwell = waves.variables['waveSpreadSwell'][:]


    output = dict()
    output['waveHs'] = waveHs
    output['waveTp'] = waveTp
    output['waveMeanDirection'] = waveMeanDirection
    output['waveTm'] = waveTm
    output['waveTm1'] = waveTm1
    output['waveTm2'] = waveTm2
    output['waveTpSwell'] = waveTpSwell
    output['waveHsSwell'] = waveHsSwell
    output['waveMeanDirectionSwell'] = waveMeanDirectionSwell
    output['waveSpreadSwell'] = waveSpreadSwell
    output['waveHsWindsea'] = waveHsWindsea
    output['waveTpWindsea'] = waveTmWindsea
    output['waveMeanDirectionWindsea'] = waveMeanDirectionWindsea
    output['waveSpreadWindsea'] = waveSpreadWindsea

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
hsSwell = []
tpSwell = []
dmSwell = []
hsWindsea = []
tpWindsea = []
dmWindsea = []

timeWave = []
for i in files_path:
    waves = getWIS(i)
    Hs = np.append(Hs,waves['waveHs'])
    Tp = np.append(Tp,waves['waveTp'])
    Dm = np.append(Dm,waves['waveMeanDirection'])
    hsSwell = np.append(hsSwell,waves['waveHsSwell'])
    tpSwell = np.append(tpSwell,waves['waveTpSwell'])
    dmSwell = np.append(dmSwell,waves['waveMeanDirectionSwell'])
    hsWindsea = np.append(hsWindsea,waves['waveHsWindsea'])
    tpWindsea = np.append(tpWindsea,waves['waveTpWindsea'])
    dmWindsea = np.append(dmWindsea,waves['waveMeanDirectionWindsea'])
    timeWave = np.append(timeWave,waves['t'].flatten())


def getArray(file):
    waves = Dataset(file)
    waveHs = waves.variables['waveHs'][:]
    waveTp = waves.variables['waveTp'][:]
    waveMeanDirection = waves.variables['waveMeanDirection'][:]
    timeW = waves.variables['time'][:]
    output = dict()
    output['waveHs'] = waveHs
    output['waveTp'] = waveTp
    output['waveMeanDirection'] = waveMeanDirection
    output['t'] = timeW
    return output

wavedir26 = '/media/dylananderson/Elements/26mArray/'
# Need to sort the files to ensure correct temporal order...
files = os.listdir(wavedir26)
files.sort()
files_path = [os.path.join(os.path.abspath(wavedir26), x) for x in files]
array26m = Dataset(files_path[0])
Hs26m = []
Tp26m = []
Dm26m = []
timeWave26m = []
for i in files_path:
    waves26m = getArray(i)
    Hs26m = np.append(Hs26m,waves26m['waveHs'][0:-1:2])
    Tp26m = np.append(Tp26m,waves26m['waveTp'][0:-1:2])
    Dm26m = np.append(Dm26m,waves26m['waveMeanDirection'][0:-1:2])
    timeWave26m = np.append(timeWave26m,waves26m['t'][0:-1:2])

ind = np.where((Hs26m > 0))
hs26m = Hs26m[ind]
tp26m = Tp26m[ind]
dm26m = Dm26m[ind]
t26m = timeWave26m[ind]
tWave26m = [DT.datetime.fromtimestamp(x) for x in t26m]


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

hsCombined = np.append(Hs,hs26m)
hsSmooth = moving_average(hsCombined,3)
tpCombined = np.append(Tp,tp26m)
dmCombined = np.append(Dm,dm26m)
tWave = [DT.datetime.fromtimestamp(x) for x in timeWave]
tC = np.append(np.array(tWave),tWave26m)




# hsCombined = Hs
# tpCombined = Tp
# dmCombined = Dm
# hsSwellCombined = hsSwell
# tpSwellCombined = tpSwell
# dmSwellCombined = dmSwell
# hsWindseaCombined = hsWindsea
# tpWindseaCombined = tpWindsea
# dmWindseaCombined = dmWindsea

# tSwell = np.array(tWave)
# tWindsea = np.array(tWave)

#
# plt.figure(figsize=(10,10))
# ax = plt.subplot2grid((3,3,),(0,0),rowspan=1,colspan=3)
# ax.plot(tWave,Hs,label='reconstruction')
# ax.plot(tWave17m,hs17m,label='17mArray')
# # ax.plot(tWave8m,hs8m,label='8mArray')
# plt.legend()
# ax2 = plt.subplot2grid((3,3,),(1,0),rowspan=1,colspan=3)
# ax2.plot(tWave,Tp)
# ax2.plot(tWave17m,tp17m)
# # ax2.plot(tWave8m,tp8m)
# ax3 = plt.subplot2grid((3,3,),(2,0),rowspan=1,colspan=3)
# ax3.plot(tWave,Dm)
# ax3.plot(tWave17m,dm17m)
# # ax3.plot(tWave8m,dm8m)
#
# plt.show()


badDirs = np.where((dmCombined > 360))
dmCombined[badDirs] = dmCombined[badDirs]*np.nan

badtp = np.where((tpCombined < 1))

tpCombined[badtp] = tpCombined[badtp]*np.nan
# badDirsSwell = np.where((dmSwellCombined > 360))
# dmSwellCombined[badDirsSwell] = dmSwellCombined[badDirsSwell]*np.nan
# badDirsWindsea = np.where((dmWindseaCombined > 360))
# dmWindseaCombined[badDirsWindsea] = dmWindseaCombined[badDirsWindsea]*np.nan

waveNorm = dmCombined - 72
neg = np.where((waveNorm > 180))
waveNorm[neg[0]] = waveNorm[neg[0]]-360
offpos = np.where((waveNorm>90))
offneg = np.where((waveNorm<-90))
waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
waveNorm[offneg[0]] = waveNorm[offneg[0]]*0
#
# waveNormSwell = dmSwellCombined - 72
# negSwell = np.where((waveNormSwell > 180))
# waveNormSwell[negSwell[0]] = waveNormSwell[negSwell[0]]-360
# offposSwell = np.where((waveNormSwell>90))
# offnegSwell = np.where((waveNormSwell<-90))
# waveNormSwell[offposSwell[0]] = waveNormSwell[offposSwell[0]]*0
# waveNormSwell[offnegSwell[0]] = waveNormSwell[offnegSwell[0]]*0
#
# waveNormWindsea = dmWindseaCombined - 72
# negWindsea = np.where((waveNormWindsea > 180))
# waveNormWindsea[negWindsea[0]] = waveNormWindsea[negWindsea[0]]-360
# offposWindsea = np.where((waveNormWindsea>90))
# offnegWindsea = np.where((waveNormWindsea<-90))
# waveNormWindsea[offposWindsea[0]] = waveNormWindsea[offposWindsea[0]]*0
# waveNormWindsea[offnegWindsea[0]] = waveNormWindsea[offnegWindsea[0]]*0



Lo = (9.81/(2*np.pi)) * np.square(tpCombined)
Ir = 0.122/(np.sqrt((hsCombined/Lo)))
HoverL = hsCombined/Lo
lwpC = 1025*np.square(hsCombined)*tpCombined*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
weC = np.square(hsCombined)*tpCombined
ws = (2.65-1)*9.81*np.square(0.00015)/(18*0.000001)
fV = hsCombined/(ws*tpCombined)

# lwpSwell = 1025*np.square(hsSwellCombined)*tpSwellCombined*(9.81/(64*np.pi))*np.cos(waveNormSwell*(np.pi/180))*np.sin(waveNormSwell*(np.pi/180))
# weSwell = np.square(hsSwellCombined)*tpSwellCombined
#
# lwpWindsea = 1025*np.square(hsWindseaCombined)*tpWindseaCombined*(9.81/(64*np.pi))*np.cos(waveNormWindsea*(np.pi/180))*np.sin(waveNormWindsea*(np.pi/180))
# weWindsea = np.square(hsWindseaCombined)*tpWindseaCombined


hsList = []
tpList = []
dmList = []
lwpList = []
weList = []
fvList = []
irList = []
hlList = []
for xx in range(len(time)-1):
    t1 = time[xx]
    t2 = time[xx+1]
    tempWave = np.where((tC < t2) & (tC > t1))

    # hsCombinedTemp = hsCombined[tempWave].astype(int)
    # array_of_tuples = map(tuple, hsCombinedTemp)
    # tuple_of_tuples = tuple(array_of_tuples)
    #
    # hsList.append(tuple_of_tuples)
    hsList.append(hsCombined[tempWave])
    tpList.append(tpCombined[tempWave])
    dmList.append(dmCombined[tempWave])
    lwpList.append(hsSmooth[tempWave])
    weList.append(weC[tempWave])
    fvList.append(fV[tempWave])
    irList.append(Ir[tempWave])
    hlList.append(HoverL[tempWave])


transitionType = np.empty((len(sorted_bmus)-1,))
# for ii in range(len(sorted_bmus)-1):
#     c = sorted_bmus[ii+1]
#     p = sorted_bmus[ii]
#     if c == 0:
#         if p == 0:
#             transitionType[ii] = 0
#         elif p > 0 and p <= 6:
#             transitionType[ii] = 2
#         elif p > 6:
#             transitionType[ii] = 1
#     elif c == 1:
#         if p == 1:
#             transitionType[ii] = 0
#         elif p > 1 and p <= 6:
#             transitionType[ii] = 2
#         elif p > 6:
#             transitionType[ii] = 1
#         elif p == 0:
#             transitionType[ii] = 1
#     elif c == 2:
#         if p == 2:
#             transitionType[ii] = 0
#         elif p > 2 and p <= 8:
#             transitionType[ii] = 2
#         elif p > 8:
#             transitionType[ii] = 1
#         elif p < 2:
#             transitionType[ii] = 1
#     elif c == 3:
#         if p == 3:
#             transitionType[ii] = 0
#         elif p > 3 and p <= 9:
#             transitionType[ii] = 2
#         elif p > 9:
#             transitionType[ii] = 1
#         elif p < 3:
#             transitionType[ii] = 1
#     elif c == 4:
#         if p == 4:
#             transitionType[ii] = 0
#         elif p > 4 and p <= 10:
#             transitionType[ii] = 2
#         elif p > 10:
#             transitionType[ii] = 1
#         elif p < 4:
#             transitionType[ii] = 1
#     elif c == 5:
#         if p == 5:
#             transitionType[ii] = 0
#         elif p > 5 and p <= 11:
#             transitionType[ii] = 2
#         elif p > 11:
#             transitionType[ii] = 1
#         elif p >= 0 and p <= 4:
#             transitionType[ii] = 1
#
#
#     elif c == 6:
#         if p == 6:
#             transitionType[ii] = 0
#         elif p > 6 and p <= 12:
#             transitionType[ii] = 2
#         elif p >= 0 and p <= 5:
#             transitionType[ii] = 1
#
#
#     elif c == 7:
#         if p == 7:
#             transitionType[ii] = 0
#         elif p > 7 and p <= 13:
#             transitionType[ii] = 2
#         elif p == 0:
#             transitionType[ii] = 2
#         elif p > 0 and p < 7:
#             transitionType[ii] = 1
#
#
#     elif c == 8:
#         if p == 8:
#             transitionType[ii] = 0
#         elif p > 8 and p <= 14:
#             transitionType[ii] = 2
#         elif p >= 0 and p <= 1:
#             transitionType[ii] = 2
#         elif p > 1 and p < 8:
#             transitionType[ii] = 1
#
#
#     elif c == 9:
#         if p == 9:
#             transitionType[ii] = 0
#         elif p > 9 and p <= 15:
#             transitionType[ii] = 2
#         elif p >= 0 and p <= 2:
#             transitionType[ii] = 2
#         elif p > 2 and p < 9:
#             transitionType[ii] = 1
#
#
#     elif c == 10:
#         if p == 10:
#             transitionType[ii] = 0
#         elif p > 10 and p <= 15:
#             transitionType[ii] = 2
#         elif p >= 0 and p <= 3:
#             transitionType[ii] = 2
#         elif p > 3 and p < 10:
#             transitionType[ii] = 1
#
#     elif c == 11:
#         if p == 11:
#             transitionType[ii] = 0
#         elif p > 11 and p <= 15:
#             transitionType[ii] = 2
#         elif p >= 0 and p <= 4:
#             transitionType[ii] = 2
#         elif p >4 and p < 11:
#             transitionType[ii] = 1
#
#     elif c == 12:
#         if p == 12:
#             transitionType[ii] = 0
#         elif p > 12 and p <= 15:
#             transitionType[ii] = 2
#         elif p >= 0 and p <= 5:
#             transitionType[ii] = 2
#         elif p > 5 and p < 12:
#             transitionType[ii] = 1
#
#     elif c == 13:
#         if p == 13:
#             transitionType[ii] = 0
#         elif p > 13 and p <= 15:
#             transitionType[ii] = 2
#         elif p >= 0 and p <= 6:
#             transitionType[ii] = 2
#         elif p > 6 and p < 13:
#             transitionType[ii] = 1




for ii in range(len(sorted_bmus)-1):
    c = sorted_bmus[ii+1]
    p = sorted_bmus[ii]
    if c == 0:
        if p == 0:
            transitionType[ii] = 0
        elif p > 0 and p <= 4:
            transitionType[ii] = 2
        elif p > 4:
            transitionType[ii] = 1
    elif c == 1:
        if p == 1:
            transitionType[ii] = 0
        elif p > 1 and p <= 5:
            transitionType[ii] = 2
        elif p > 5:
            transitionType[ii] = 1
        elif p < 1:
            transitionType[ii] = 2
    elif c == 2:
        if p == 2:
            transitionType[ii] = 0
        elif p > 2 and p < 6:
            transitionType[ii] = 2
        elif p >= 6:
            transitionType[ii] = 1
        elif p < 2:
            transitionType[ii] = 1
    elif c == 3:
        if p == 3:
            transitionType[ii] = 0
        elif p > 3 and p < 7:
            transitionType[ii] = 2
        elif p >= 7:
            transitionType[ii] = 1
        elif p < 3:
            transitionType[ii] = 1
    elif c == 4:
        if p == 4:
            transitionType[ii] = 0
        elif p > 4 and p < 8:
            transitionType[ii] = 2
        elif p == 0:
            transitionType[ii] = 2
        elif p < 4 and p > 0:
            transitionType[ii] = 1
    elif c == 5:
        if p == 5:
            transitionType[ii] = 0
        elif p > 5 and p < 9:
            transitionType[ii] = 2
        elif p >= 0 and p <= 1:
            transitionType[ii] = 2
        elif p > 1 and p < 5:
            transitionType[ii] = 1
    elif c == 6:
        if p == 6:
            transitionType[ii] = 0
        elif p > 6 and p < 9:
            transitionType[ii] = 2
        elif p >= 0 and p <= 2:
            transitionType[ii] = 2
        elif p > 2 and p < 6:
            transitionType[ii] = 1
    elif c == 7:
        if p == 7:
            transitionType[ii] = 0
        elif p >= 0 and p <= 3:
            transitionType[ii] = 2
        elif p > 3 and p < 7:
            transitionType[ii] = 1







sameIndices = np.where((transitionType==0))
upIndices = np.where((transitionType==1))
downIndices = np.where((transitionType==2))

upHsList = []
testing = [upHsList.append(hsList[xx]) for xx in upIndices[0]]
upHs = np.concatenate(upHsList,axis=0)
downHsList = []
testing = [downHsList.append(hsList[xx]) for xx in downIndices[0]]
downHs = np.concatenate(downHsList,axis=0)
sameHsList = []
testing = [sameHsList.append(hsList[xx]) for xx in sameIndices[0]]
sameHs = np.concatenate(sameHsList,axis=0)


upFvList = []
testing = [upFvList.append(fvList[xx]) for xx in upIndices[0]]
upFv = np.concatenate(upFvList,axis=0)
downFvList = []
testing = [downFvList.append(fvList[xx]) for xx in downIndices[0]]
downFv = np.concatenate(downFvList,axis=0)
sameFvList = []
testing = [sameFvList.append(fvList[xx]) for xx in sameIndices[0]]
sameFv = np.concatenate(sameFvList,axis=0)


import scipy
# function for calculating the t-test for two independent samples
def independent_ttest(data1, data2, alpha):
    # calculate means
    mean1, mean2 = np.nanmean(data1), np.nanmean(data2)
    # calculate standard errors
    se1, se2 = scipy.stats.sem(data1), scipy.stats.sem(data2)
    # standard error on the difference between the samples
    sed = scipy.sqrt(se1 ** 2.0 + se2 ** 2.0)
    # calculate the t statistic
    t_stat = (mean1 - mean2) / sed
    # degrees of freedom
    df = len(data1) + len(data2) - 2
    # calculate the critical value
    cv = scipy.stats.t.ppf(1.0 - alpha, df)
    # calculate the p-value
    p = (1.0 - scipy.stats.t.cdf(abs(t_stat), df)) * 2.0
    # return everything
    return t_stat, df, cv, p






#
# import matplotlib.cm as cm
# import matplotlib.colors as mcolors
# plt.style.use('default')
#
# from scipy.stats.kde import gaussian_kde
# dist_space = np.linspace(0, 4, 80)
# fig = plt.figure(figsize=(10,10))
# gs2 = gridspec.GridSpec(6, 5)
#
# import matplotlib.cm as cm
# import matplotlib.colors as mcolors
# from scipy.stats.kde import gaussian_kde
# dist_space = np.linspace(0, 5, 50)
# fig = plt.figure(figsize=(10,10))
# colorparam = np.zeros((numClusters*numClusters,))
# counter = 0
# for xx in range(numClusters):
#     for yy in range(numClusters):
#         ax = plt.subplot2grid((numClusters,numClusters), (yy,xx), rowspan=1, colspan=1)
#         #normalize = mcolors.Normalize(vmin=np.min(colorparams), vmax=np.max(colorparams))
#         normalize = mcolors.Normalize(vmin=.5, vmax=1.5)
#
#         ax.set_xlim([0,3])
#         ax.set_ylim([0,2])
#         data = wHs[xx][yy]
#         if len(data)>0:
#             kde = gaussian_kde(data)
#             colorparam[counter] = np.nanmean(data)
#             colormap = cm.Reds
#             color = colormap(normalize(colorparam[counter]))
#             ax.plot(dist_space, kde(dist_space),linewidth=1,color=color)
#             ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
#             ax.spines['top'].set_color([0.5, 0.5, 0.5])
#             ax.spines['right'].set_color([0.5, 0.5, 0.5])
#             ax.spines['left'].set_color([0.5, 0.5, 0.5])
#             #ax.text(1.8, 1, np.round(colorparam*100)/100, fontweight='bold')
#
#         else:
#             ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
#             ax.spines['top'].set_color([0.3, 0.3, 0.3])
#             ax.spines['right'].set_color([0.3, 0.3, 0.3])
#             ax.spines['left'].set_color([0.3, 0.3, 0.3])
#             if yy < 14:
#                 ax.xaxis.set_ticks([])
#         if yy < 14:
#             ax.xaxis.set_ticklabels([])
#         ax.yaxis.set_ticklabels([])
#         ax.yaxis.set_ticks([])
#         counter = counter+1
# plt.show()
# s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
# s_map.set_array(colorparam)
# fig.subplots_adjust(right=0.92)
# cbar_ax = fig.add_axes([0.94, 0.10, 0.02, 0.7])
# cbar = fig.colorbar(s_map, cax=cbar_ax)
# cbar.set_label('Mean Hs (m)')
# #cbar = ax.cax.colorbar(s_map, ticks=colorparam, format='%2.2g') # format='%2i' for integer
#
#






