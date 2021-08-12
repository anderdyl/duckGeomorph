
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
plt.figure()
plt.plot(days)
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



m,m2 = transition_matrix(sorted_bmus)
for row in m: print(' '.join('{0:.2f}'.format(x) for x in row))

flat_list = [item for sublist in m for item in sublist]
flatarray = np.asarray(flat_list)
flatarray.resize(numClusters, numClusters)

flat_list2 = [item for sublist in m2 for item in sublist]
flatarray2 = np.asarray(flat_list2)
flatarray2.resize(numClusters, numClusters)






plt.figure(figsize=(10,10))
plt.plot(time,sorted_bmus)
import datetime
years = np.arange(1979,2020)
summers = []
summerTransitions = []
for year in years:
    summerInd = np.where((time>datetime.datetime(year=year,month=6,day=1)) & (time < datetime.datetime(year=year,month=9,day=1)))
    summers.append(summerInd)
    if len(summerInd[0]>1):
        tran1, tran2 = conditional_transition_matrix(sorted_bmus[summerInd],numClusters)
        flat_listCond = [item for sublist in tran2 for item in sublist]
        flatarrayCond = np.asarray(flat_listCond)
        flatarrayCond.resize(numClusters, numClusters)

        summerTransitions.append(flatarrayCond)

summerAllTransitions = np.sum(summerTransitions,axis=0)

falls = []
fallTransitions = []

for year in years:
    fallInd = np.where((time>datetime.datetime(year=year,month=9,day=1)) & (time < datetime.datetime(year=year,month=12,day=1)))
    falls.append(fallInd)
    if len(fallInd[0]>1):
        tran1, tran2 = conditional_transition_matrix(sorted_bmus[fallInd],numClusters)
        flat_listCond = [item for sublist in tran2 for item in sublist]
        flatarrayCond = np.asarray(flat_listCond)
        flatarrayCond.resize(numClusters, numClusters)

        fallTransitions.append(flatarrayCond)
fallAllTransitions = np.sum(fallTransitions,axis=0)

winters = []
winterTransitions = []
for year in years:
    winterInd = np.where((time>datetime.datetime(year=year,month=12,day=1)) & (time < datetime.datetime(year=year+1,month=3,day=1)))
    winters.append(winterInd)
    if len(winterInd[0]>1):
        tran1, tran2 = conditional_transition_matrix(sorted_bmus[winterInd],numClusters)
        flat_listCond = [item for sublist in tran2 for item in sublist]
        flatarrayCond = np.asarray(flat_listCond)
        flatarrayCond.resize(numClusters, numClusters)

        winterTransitions.append(flatarrayCond)
winterAllTransitions = np.sum(winterTransitions,axis=0)


springs = []
springTransitions = []
for year in years:
    springInd = np.where((time>datetime.datetime(year=year,month=3,day=1)) & (time < datetime.datetime(year=year,month=6,day=1)))
    springs.append(springInd)
    if len(springInd[0]>1):
        tran1, tran2 = conditional_transition_matrix(sorted_bmus[springInd],numClusters)
        flat_listCond = [item for sublist in tran2 for item in sublist]
        flatarrayCond = np.asarray(flat_listCond)
        flatarrayCond.resize(numClusters, numClusters)

        springTransitions.append(flatarrayCond)
springAllTransitions = np.sum(springTransitions,axis=0)

def convertToProbabilities(totalTransitions):
    mm, nn = np.shape(totalTransitions)
    output = np.zeros((mm,nn))
    for qq in range(mm):
        s = sum(totalTransitions[qq])
        if s > 0:
            output[qq,0:len(totalTransitions[qq])] = [f / s for f in totalTransitions[qq]]
    return output


summerProbs = convertToProbabilities(summerAllTransitions)
springProbs = convertToProbabilities(springAllTransitions)
winterProbs = convertToProbabilities(winterAllTransitions)
fallProbs = convertToProbabilities(fallAllTransitions)



import pickle
dbfile = open('filteredStorms.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()

dataCop = data['dataCop']
filteredHs = data['filteredHs']
filteredTp = data['filteredTp']
filteredDm = data['filteredDm']
filteredNTR = data['filteredNTR']
filteredDur = data['filteredDur']
filteredTime = data['filteredTime']


upperThreshold = 10
lowerThreshold = 5
stormSubset = np.where((filteredHs > lowerThreshold) & (filteredHs < upperThreshold))
timeSubset = filteredTime[stormSubset]
waveSubset = filteredHs[stormSubset]
stormTransitions = []
for xx in range(len(stormSubset[0])):
    targetTime = timeSubset[xx]

    if targetTime > time[0]:
        earlyTime = targetTime - datetime.timedelta(days=240)
        earlierSurveys = np.where((time>earlyTime) & (time<targetTime))
        laterTime = targetTime + datetime.timedelta(days=240)
        laterSurveys = np.where((time<laterTime) & (time>targetTime))
        startSurvey = sorted_bmus[earlierSurveys[0][-1]]
        endSurvey = sorted_bmus[laterSurveys[0][0]]
        startSurveyTime = time[earlierSurveys[0][-1]]
        endSurveyTime = time[laterSurveys[0][0]]
        checkForStorms = np.where((filteredTime > startSurveyTime) & (filteredTime < endSurveyTime))
        wavesForStorms = filteredHs[checkForStorms]
        if np.max(wavesForStorms) == waveSubset[xx]:
            temp = np.hstack((startSurvey, endSurvey))
            tran1, tran2 = conditional_transition_matrix(temp,numClusters)
            flat_listCond = [item for sublist in tran2 for item in sublist]
            flatarrayCond = np.asarray(flat_listCond)
            flatarrayCond.resize(numClusters, numClusters)
            stormTransitions.append(flatarrayCond)

hugeStormTransitions = np.sum(stormTransitions,axis=0)



upperThreshold = 5
lowerThreshold = 4
stormSubset = np.where((filteredHs > lowerThreshold) & (filteredHs < upperThreshold))
timeSubset = filteredTime[stormSubset]
waveSubset = filteredHs[stormSubset]
stormTransitions = []
for xx in range(len(stormSubset[0])):
    targetTime = timeSubset[xx]

    if targetTime > time[0]:
        earlyTime = targetTime - datetime.timedelta(days=240)
        earlierSurveys = np.where((time>earlyTime) & (time<targetTime))
        laterTime = targetTime + datetime.timedelta(days=240)
        laterSurveys = np.where((time<laterTime) & (time>targetTime))
        startSurvey = sorted_bmus[earlierSurveys[0][-1]]
        endSurvey = sorted_bmus[laterSurveys[0][0]]
        startSurveyTime = time[earlierSurveys[0][-1]]
        endSurveyTime = time[laterSurveys[0][0]]
        checkForStorms = np.where((filteredTime > startSurveyTime) & (filteredTime < endSurveyTime))
        wavesForStorms = filteredHs[checkForStorms]
        if np.max(wavesForStorms) == waveSubset[xx]:
            temp = np.hstack((startSurvey, endSurvey))
            tran1, tran2 = conditional_transition_matrix(temp,numClusters)
            flat_listCond = [item for sublist in tran2 for item in sublist]
            flatarrayCond = np.asarray(flat_listCond)
            flatarrayCond.resize(numClusters, numClusters)
            stormTransitions.append(flatarrayCond)

highStormTransitions = np.sum(stormTransitions,axis=0)


upperThreshold = 4
lowerThreshold = 3
stormSubset = np.where((filteredHs > lowerThreshold) & (filteredHs < upperThreshold))
timeSubset = filteredTime[stormSubset]
waveSubset = filteredHs[stormSubset]
stormTransitions = []
for xx in range(len(stormSubset[0])):
    targetTime = timeSubset[xx]

    if targetTime > time[0]:
        earlyTime = targetTime - datetime.timedelta(days=240)
        earlierSurveys = np.where((time>earlyTime) & (time<targetTime))
        laterTime = targetTime + datetime.timedelta(days=240)
        laterSurveys = np.where((time<laterTime) & (time>targetTime))
        startSurvey = sorted_bmus[earlierSurveys[0][-1]]
        endSurvey = sorted_bmus[laterSurveys[0][0]]
        startSurveyTime = time[earlierSurveys[0][-1]]
        endSurveyTime = time[laterSurveys[0][0]]
        checkForStorms = np.where((filteredTime > startSurveyTime) & (filteredTime < endSurveyTime))
        wavesForStorms = filteredHs[checkForStorms]
        if np.max(wavesForStorms) == waveSubset[xx]:
            temp = np.hstack((startSurvey, endSurvey))
            tran1, tran2 = conditional_transition_matrix(temp,numClusters)
            flat_listCond = [item for sublist in tran2 for item in sublist]
            flatarrayCond = np.asarray(flat_listCond)
            flatarrayCond.resize(numClusters, numClusters)
            stormTransitions.append(flatarrayCond)

midStormTransitions = np.sum(stormTransitions,axis=0)



upperThreshold = 3
lowerThreshold = 2
stormSubset = np.where((filteredHs > lowerThreshold) & (filteredHs < upperThreshold))
timeSubset = filteredTime[stormSubset]
waveSubset = filteredHs[stormSubset]
stormTransitions = []
for xx in range(len(stormSubset[0])):
    targetTime = timeSubset[xx]

    if targetTime > time[0]:
        earlyTime = targetTime - datetime.timedelta(days=240)
        earlierSurveys = np.where((time>earlyTime) & (time<targetTime))
        laterTime = targetTime + datetime.timedelta(days=240)
        laterSurveys = np.where((time<laterTime) & (time>targetTime))
        startSurvey = sorted_bmus[earlierSurveys[0][-1]]
        endSurvey = sorted_bmus[laterSurveys[0][0]]
        startSurveyTime = time[earlierSurveys[0][-1]]
        endSurveyTime = time[laterSurveys[0][0]]
        checkForStorms = np.where((filteredTime > startSurveyTime) & (filteredTime < endSurveyTime))
        wavesForStorms = filteredHs[checkForStorms]
        if np.max(wavesForStorms) == waveSubset[xx]:
            temp = np.hstack((startSurvey, endSurvey))
            tran1, tran2 = conditional_transition_matrix(temp,numClusters)
            flat_listCond = [item for sublist in tran2 for item in sublist]
            flatarrayCond = np.asarray(flat_listCond)
            flatarrayCond.resize(numClusters, numClusters)
            stormTransitions.append(flatarrayCond)

lowStormTransitions = np.sum(stormTransitions,axis=0)


# lets find everything with no storms....
noStormTransitions = []
for xx in range(len(time)-1):
    earlyTime = time[xx]
    laterTime = time[xx+1]
    checkForStorms = np.where((filteredTime > earlyTime) & (filteredTime < laterTime))
    #print(checkForStorms)
    if len(checkForStorms[0]) == 0:
        s1 = sorted_bmus[xx]
        s2 = sorted_bmus[xx+1]
        temp = np.hstack((s1, s2))
        tran1, tran2 = conditional_transition_matrix(temp, numClusters)
        flat_listCond = [item for sublist in tran2 for item in sublist]
        flatarrayCond = np.asarray(flat_listCond)
        flatarrayCond.resize(numClusters, numClusters)
        noStormTransitions.append(flatarrayCond)
noStorms = np.sum(noStormTransitions,axis=0)







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

residualWaterLevelFRF = residualWaterLevelFRF-np.nanmean(residualWaterLevelFRF)






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
    Hs26m = np.append(Hs26m,waves26m['waveHs'])
    Tp26m = np.append(Tp26m,waves26m['waveTp'])
    Dm26m = np.append(Dm26m,waves26m['waveMeanDirection'])
    timeWave26m = np.append(timeWave26m,waves26m['t'])

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

plt.style.use('default')
plt.figure(figsize=(10,10))
ax = plt.subplot2grid((3,3,),(0,0),rowspan=1,colspan=3)
ax.plot(tC,hsCombined,color='k')
ax.set_ylabel('Hs (m)')

# ax.plot(tWave8m,hs8m,label='8mArray')
ax2 = plt.subplot2grid((3,3,),(1,0),rowspan=1,colspan=3)
ax2.plot(tC,tpCombined,color='k')
ax2.set_ylabel('Tp (s)')
ax2.set_ylim([0, 20])

# ax2.plot(tWave8m,tp8m)
ax3 = plt.subplot2grid((3,3,),(2,0),rowspan=1,colspan=3)
ax3.scatter(tC,dmCombined,c='k',s=1)
ax3.set_ylabel('Dir (deg)')
ax3.set_ylim([0, 360])
# ax3.plot(tWave8m,dm8m)

plt.show()




fig3 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(1, numClusters, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
import matplotlib.cm as cm
#colors = cm.rainbow(np.linspace(0, 1, numClusters))
colors = cm.hsv(np.linspace(0, 1, (numClusters+1)))

for i in range(numClusters):
    #getind = sorted[i]
    getind = order[i] #np.where(order==i+1)
    getProfs = np.where(bmus==getind)
    true = np.nanmean(alllines[getProfs,:],axis=1)

    #dev = deviation[sortedPeakInd[i],:]
    #true = profiles[sortedPeakInd[i],:]

    #peaks = sig.find_peaks(x=(true), prominence=0.05)

    #if len(peaks[0]) > 0:
    #    offshorePeaks[i] = np.max(peaks[0])

    ax = plt.subplot(gs[gc, gr])
    #ax.plot(xinterp,np.mean(alllines,axis=0),'w-')

    ax.plot(xinterp,true[0],color=colors[i])
    #ax.plot(xinterp[peaks[0]],true[peaks[0]],'ro')
    ax.set_xlim([0, 500])
    ax.set_ylim([-8, 0.5])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    #ax.text(400,0, group_size[sortedPeakInd[i]], fontweight='bold')
    #ax.text(400,0, sortedPhases[i], fontweight='bold')

    if gr > 0:
        ax.set_yticks([])

    # if gc < numClusters:
    #     ax.set_xticks([])
    #  counter
    gr += 1
    if gc >= numClusters:
        gc += 1
        gr = 0




fig3 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(numClusters, 1, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
import matplotlib.cm as cm
#colors = cm.rainbow(np.linspace(0, 1, numClusters))
colors = cm.hsv(np.linspace(0, 1, (numClusters+1)))

for i in range(numClusters):
    #getind = sorted[i]
    getind = order[i] #np.where(order==i+1)
    getProfs = np.where(bmus==getind)
    true = np.nanmean(alllines[getProfs,:],axis=1)

    #dev = deviation[sortedPeakInd[i],:]
    #true = profiles[sortedPeakInd[i],:]

    #peaks = sig.find_peaks(x=(true), prominence=0.05)

    #if len(peaks[0]) > 0:
    #    offshorePeaks[i] = np.max(peaks[0])

    ax = plt.subplot(gs[gc, gr])
    #ax.plot(xinterp,np.mean(alllines,axis=0),'w-')

    ax.plot(xinterp,true[0],color=colors[i])
    #ax.plot(xinterp[peaks[0]],true[peaks[0]],'ro')
    ax.set_xlim([0, 500])
    ax.set_ylim([-8, 0.5])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    #ax.text(400,0, group_size[sortedPeakInd[i]], fontweight='bold')
    #ax.text(400,0, sortedPhases[i], fontweight='bold')

    if gr > 0:
        ax.set_yticks([])

    # if gc < numClusters:
    #     ax.set_xticks([])
    #  counter
    gc += 1
    if gc >= numClusters:
        gc += 1
        gr = 0

def getDistributionPercentages(data,numClusters):

    #minData = np.nanmin(data)
    #maxData = np.nanmax(data)
    #dist_space = np.linspace(minData,maxData,numBins)
    param05 = np.zeros((numClusters*numClusters,))
    param25 = np.zeros((numClusters*numClusters,))
    param50 = np.zeros((numClusters*numClusters,))
    param75 = np.zeros((numClusters*numClusters,))
    param90 = np.zeros((numClusters*numClusters,))
    param95 = np.zeros((numClusters*numClusters,))

    counter = 0
    for xx in range(numClusters):
        for yy in range(numClusters):

            transitionData = data[xx][yy]
            if len(data)>0:
                #kdes = gaussian_kde(transitionData)
                #param50[counter] = np.nanmean(transitionData)
                param05[counter] = np.nanpercentile(transitionData,5)
                param25[counter] = np.nanpercentile(transitionData,25)
                param50[counter] = np.nanpercentile(transitionData,50)
                param75[counter] = np.nanpercentile(transitionData,75)
                param90[counter] = np.nanpercentile(transitionData,90)
                param95[counter] = np.nanpercentile(transitionData,95)

                #kde_ci = np.array(np.percentile(kdes, (5, 95), axis=0))
            counter = counter+1

    output = dict()
    output['param05'] = param05
    output['param25'] = param25
    output['param50'] = param50
    output['param75'] = param75
    output['param90'] = param90
    output['param95'] = param95

    return output



def getTransitionConditions(data,dataTime,numClusters,prevbin,date,waveTime):

    param = []
    param05 = []
    param25 = []
    param50 = []
    param75 = []
    param90 = []
    param95 = []

    for xx in range(numClusters):
        innerListParam = []
        innerListParam05 = []
        innerListParam25 = []
        innerListParam50 = []
        innerListParam75 = []
        innerListParam90 = []
        innerListParam95 = []

        for yy in range(numClusters):
            # if yy is different...
            wInd = np.where(prevbin[xx] == yy)
            if len(wInd[0]) > 0:
                tempParam = []
                tempParam05 = []
                tempParam25 = []
                tempParam50 = []
                tempParam75 = []
                tempParam90 = []
                tempParam95 = []

                for tt in range(len(wInd[0])):
                    if daysBetween[xx][wInd[0][tt]] < 92:
                        tempT = np.where((dataTime < date[xx][wInd[0][tt]]) & (dataTime > prevdate[xx][wInd[0][tt]]))
                        tempWave = np.where((waveTime < date[xx][wInd[0][tt]]) & (waveTime > prevdate[xx][wInd[0][tt]]))

                        if len(tempT[0]) > 0 and len(tempWave[0]) > 0:
                            tempParam = np.append(tempParam, data[tempT])
                            tempParam05 = np.append(tempParam05, np.nanpercentile(data[tempT],5))
                            tempParam25 = np.append(tempParam25, np.nanpercentile(data[tempT],25))
                            tempParam50 = np.append(tempParam50, np.nanpercentile(data[tempT],50))
                            tempParam75 = np.append(tempParam75, np.nanpercentile(data[tempT],75))
                            tempParam90 = np.append(tempParam90, np.nanpercentile(data[tempT],90))
                            tempParam95 = np.append(tempParam95, np.nanpercentile(data[tempT],95))
                        elif len(tempT[0]) == 0 and len(tempWave[0]) > 0:
                            tempParam = np.append(tempParam, -99)
                            tempParam05 = np.append(tempParam05, -99)
                            tempParam25 = np.append(tempParam25, -99)
                            tempParam50 = np.append(tempParam50, -99)
                            tempParam75 = np.append(tempParam75, -99)
                            tempParam90 = np.append(tempParam90, -99)
                            tempParam95 = np.append(tempParam95, -99)

            else:
                tempParam = []
                tempParam05 = []
                tempParam25 = []
                tempParam50 = []
                tempParam75 = []
                tempParam90 = []
                tempParam95 = []

            innerListParam.append(tempParam)
            innerListParam05.append(tempParam05)
            innerListParam25.append(tempParam25)
            innerListParam50.append(tempParam50)
            innerListParam75.append(tempParam75)
            innerListParam90.append(tempParam90)
            innerListParam95.append(tempParam95)

        param.append(innerListParam)
        param05.append(innerListParam05)
        param25.append(innerListParam25)
        param50.append(innerListParam50)
        param75.append(innerListParam75)
        param90.append(innerListParam90)
        param95.append(innerListParam95)

    output = dict()
    output['param'] = param
    output['param05'] = param05
    output['param25'] = param25
    output['param50'] = param50
    output['param75'] = param75
    output['param90'] = param90
    output['param95'] = param95
    return output





hsSplit = getTransitionConditions(hsSmooth, tC, numClusters, prevbin, date, tC)
tpSplit = getTransitionConditions(tpCombined, tC, numClusters, prevbin, date, tC)
dmSplit = getTransitionConditions(dmCombined, tC, numClusters, prevbin, date, tC)
lwpSplit = getTransitionConditions(np.abs(lwpC), tC, numClusters, prevbin, date, tC)
weSplit = getTransitionConditions(weC, tC, numClusters, prevbin, date, tC)
hoverlSplit = getTransitionConditions(HoverL, tC, numClusters, prevbin, date, tC)
loSplit = getTransitionConditions(Lo, tC, numClusters, prevbin, date, tC)
irSplit = getTransitionConditions(Ir, tC, numClusters, prevbin, date, tC)
fvSplit = getTransitionConditions(fV, tC, numClusters, prevbin, date, tC)

waterLevelSplit = getTransitionConditions(waterLevelFRF, tWaterLevelFRF, numClusters, prevbin, date, tC)
predWaterLevelSplit = getTransitionConditions(predictedWaterLevelFRF, tWaterLevelFRF, numClusters, prevbin, date, tC)
resWaterLevelSplit = getTransitionConditions(residualWaterLevelFRF, tWaterLevelFRF, numClusters, prevbin, date, tC)

sustWindSpeedSplit = getTransitionConditions(sustWindSpeed, tWinds, numClusters, prevbin, date, tC)
windSpeedSplit = getTransitionConditions(windSpeed, tWinds, numClusters, prevbin, date, tC)
windDirectionSplit = getTransitionConditions(windDirection, tWinds, numClusters, prevbin, date, tC)


from itertools import groupby
from operator import itemgetter
import more_itertools as mit

wCurrentBmu = []
wPrevBmu = []

wHs = []        # Offshore Wave Heights
wTp = []        # Offshore Wave Periods
wDm = []        # Offshore Wave Directions
wLo = []        # Deep Water Wave Length
wLWP = []       # Offshore Longshore Wave Power
wWE = []        # Offshore Wave Energy
wT = []
wIr = []        # Offshore Iribarren Number
wHoverL = []    # Offshore wave steepness

wCumuWE = []    # Cumulative Offshore Wave Energy
wCumuLWP = []   # Cumulative Offshore Longshore Wave Power


wStorm = []
wStormHs = []
wStormCumuWE = []
wStormCumuLWP = []

wCalm = []
wRatio = []
wNumStorms = []

wAvgHs = []
wAvgTp = []
wAvgWE = []
wAvgLWP = []
wAvgIr = []
wDays = []

# swellHs = []
# swellTp = []
# swellDm = []
# swellLWP = []
# windseaHs = []
# windseaTp = []
# windseaDm = []
# windseaLWP = []
for xx in range(numClusters):
    innerListCurrentBmu = []
    innerListPrevBmu = []
    innerListHs = []
    innerListAvgHs = []
    innerListTp = []
    innerListAvgTp = []
    innerListDm = []
    innerListHoverL = []
    innerListLo = []
    innerListLWP = []
    innerListCumuLWP = []
    innerListAvgLWP = []
    innerListWE = []
    innerListCumuWE = []
    innerListAvgWE = []
    innerListT = []
    innerListStorm = []
    innerListCalm = []
    innerListRatio = []
    innerListNumStorms = []
    innerListIr = []
    innerListAvgIr = []
    innerListDays = []
    innerListStormHs = []
    innerListStormCumuWE = []
    innerListStormCumuLWP = []

    # innerListSwellHs = []
    # innerListSwellTp = []
    # innerListSwellDm = []
    # innerListSwellLWP = []
    # innerListWindseaHs = []
    # innerListWindseaTp = []
    # innerListWindseaDm = []
    # innerListWindseaLWP = []
    for yy in range(numClusters):
        # if both are equal then the wave conditions didn't cause a transition and
        # everything between the last two dates should be added to a distribution

        # if yy is different...
        wInd = np.where(prevbin[xx]==yy)
        if len(wInd[0])>0:
            tempCurrentBmu = []
            tempPrevBmu = []
            tempHs = []
            tempAvgHs = []
            tempTp = []
            tempAvgTp = []
            tempDm = []
            tempHoverL = []
            tempLo = []
            tempLWP = []
            tempCumuLWP = []
            tempAvgLWP = []
            tempWE = []
            tempCumuWE = []
            tempAvgWE = []
            tempTi = []
            tempStorm = []
            tempCalm = []
            tempRatio =[]
            tempNumStorms =[]
            tempIr = []
            tempDays = []
            tempStormHs = []
            tempStormCumuWE = []
            tempStormCumuLWP = []
            tempAvgIr = []
            # tempSwellHs = []
            # tempSwellTp = []
            # tempSwellDm = []
            # tempSwellLWP = []
            # tempWindseaHs = []
            # tempWindseaTp = []
            # tempWindseaDm = []
            # tempWindseaLWP = []
            for tt in range(len(wInd[0])):
                if daysBetween[xx][wInd[0][tt]] < 92:
                    tempT = np.where((tC < date[xx][wInd[0][tt]]) & (tC > prevdate[xx][wInd[0][tt]]))

                    if len(tempT[0]) > 0:

                        # lets look for storms in the 26 m record
                        stormHsInd = np.where((hsSmooth[tempT]>2))
                        stormHsList = [list(group) for group in mit.consecutive_groups(stormHsInd[0])]
                        c = 0
                        c1 = 0
                        for qq in range(len(stormHsList)):
                            innerInnerCumuWE = []
                            innerInnerCumuLWP = []
                            innerInnerStormHs = []
                            if len(stormHsList[qq])>12:
                                c = c + len(stormHsList[qq])
                                innerInnerStormHs = np.append(innerInnerStormHs,hsSmooth[tempT][stormHsList[qq]])
                                innerInnerCumuWE = np.append(innerInnerCumuWE,lwpC[tempT][stormHsList[qq]])
                                innerInnerCumuLWP = np.append(innerInnerCumuLWP,weC[tempT][stormHsList[qq]])

                                c1 = c1 + 1

                        if c == 0:
                            tempStormHs = np.append(tempStormHs, 0)
                            tempStormCumuWE = np.append(tempStormCumuWE, 0)
                            tempStormCumuLWP = np.append(tempStormCumuLWP, 0)
                        else:
                            tempStormHs = np.append(tempStormHs,np.mean(innerInnerStormHs))
                            tempStormCumuWE = np.append(tempStormCumuWE,np.sum(innerInnerCumuWE))
                            tempStormCumuLWP = np.append(tempStormCumuLWP,np.sum(np.abs(innerInnerCumuLWP)))


                        calmC = len(tempT[0])-c
                        tempHs = np.append(tempHs,hsSmooth[tempT])

                        tempAvgHs = np.append(tempAvgHs,np.nanmean(hsSmooth[tempT]))
                        # tempHs = np.append(tempHs,hsCombined[tempT])
                        tempDays = np.append(tempDays,daysBetween[xx][wInd[0][tt]])
                        tempCurrentBmu = np.append(tempCurrentBmu,xx)
                        tempPrevBmu = np.append(tempPrevBmu,yy)
                        tempTp = np.append(tempTp,tpCombined[tempT])
                        tempAvgTp = np.append(tempAvgTp,np.nanmean(tpCombined[tempT]))
                        tempDm = np.append(tempDm,waveNorm[tempT])
                        tempHoverL = np.append(tempHoverL,HoverL[tempT])
                        tempLo = np.append(tempLo,Lo[tempT])
                        tempLWP = np.append(tempLWP,lwpC[tempT])
                        tempCumuLWP = np.append(tempCumuLWP, np.nansum(np.abs(lwpC[tempT])))
                        tempAvgLWP = np.append(tempAvgLWP, np.nanmean(np.abs(lwpC[tempT])))
                        tempWE = np.append(tempWE,weC[tempT])
                        tempCumuWE = np.append(tempCumuWE, np.nansum(weC[tempT]))
                        tempAvgWE = np.append(tempAvgWE, np.nanmean(weC[tempT]))
                        tempTi = np.append(tempTi,tC[tempT])
                        tempStorm = np.append(tempStorm,c)
                        tempCalm = np.append(tempCalm,calmC)
                        if len(tempT[0])>0:
                            tempRatio = np.append(tempRatio,c/calmC)
                        else:
                            tempRatio = np.append(tempRatio,0)
                        tempNumStorms = np.append(tempNumStorms,c1)
                        tempIr = np.append(tempIr,Ir[tempT])
                        tempAvgIr = np.append(tempAvgIr,np.nanmean(Ir[tempT]))

                        # tempSwellHs = np.append(tempSwellHs,hsSwellCombined[tempT])
                        # tempSwellTp = np.append(tempSwellTp,tpSwellCombined[tempT])
                        # tempSwellDm = np.append(tempSwellDm,waveNormSwell[tempT])
                        # tempSwellLWP = np.append(tempSwellLWP,lwpSwell[tempT])
                        # tempWindseaHs = np.append(tempWindseaHs,hsWindseaCombined[tempT])
                        # tempWindseaTp = np.append(tempWindseaTp,tpWindseaCombined[tempT])
                        # tempWindseaDm = np.append(tempWindseaDm,waveNormWindsea[tempT])
                        # tempWindseaLWP = np.append(tempWindseaLWP,lwpWindsea[tempT])
        else:
            tempCurrentBmu = []
            tempPrevBmu = []
            tempHs = []
            tempAvgHs = []
            tempTp = []
            tempAvgTp = []
            tempHoverL = []
            tempDm = []
            tempLo = []
            tempLWP = []
            tempCumuLWP = []
            tempAvgLWP = []
            tempWE = []
            tempCumuWE = []
            tempAvgWE = []
            tempTi = []
            tempStorm = []
            tempCalm = []
            tempRatio = []
            tempNumStorms = []
            tempIr = []
            tempDays = []
            tempStormHs = []
            tempStormCumuLWP = []
            tempStormCumuWE = []
            tempAvgIr = []

            # tempSwellHs = []
            # tempSwellTp = []
            # tempSwellDm = []
            # tempSwellLWP = []
            # tempWindseaHs = []
            # tempWindseaTp = []
            # tempWindseaDm = []
            # tempWindseaLWP = []
        innerListCurrentBmu.append(tempCurrentBmu)
        innerListPrevBmu.append(tempPrevBmu)
        innerListHs.append(tempHs)
        innerListAvgHs.append(tempAvgHs)
        innerListTp.append(tempTp)
        innerListAvgTp.append(tempAvgTp)
        innerListHoverL.append(tempHoverL)
        innerListDm.append(tempDm)
        innerListLo.append(tempLo)
        innerListLWP.append(tempLWP)
        innerListCumuLWP.append(tempCumuLWP)
        innerListAvgLWP.append(tempAvgLWP)
        innerListWE.append(tempWE)
        innerListCumuWE.append(tempCumuWE)
        innerListAvgWE.append(tempAvgWE)
        innerListT.append(tempTi)
        innerListStorm.append(tempStorm)
        innerListCalm.append(tempCalm)
        innerListRatio.append(tempRatio)
        innerListNumStorms.append(tempNumStorms)
        innerListIr.append(tempIr)
        innerListAvgIr.append(tempAvgIr)

        innerListDays.append(tempDays)
        innerListStormHs.append(tempStormHs)
        innerListStormCumuWE.append(tempStormCumuWE)
        innerListStormCumuLWP.append(tempStormCumuLWP)
        # innerListSwellHs.append(tempSwellHs)
        # innerListSwellTp.append(tempSwellTp)
        # innerListSwellDm.append(tempSwellDm)
        # innerListSwellLWP.append(tempSwellLWP)
        # innerListWindseaHs.append(tempWindseaHs)
        # innerListWindseaTp.append(tempWindseaTp)
        # innerListWindseaDm.append(tempWindseaDm)
        # innerListWindseaLWP.append(tempWindseaLWP)
    wDays.append(innerListDays)
    wCurrentBmu.append(innerListCurrentBmu)
    wPrevBmu.append(innerListPrevBmu)
    wHs.append(innerListHs)
    wAvgHs.append(innerListAvgHs)
    wTp.append(innerListTp)
    wAvgTp.append(innerListAvgTp)
    wDm.append(innerListDm)
    wHoverL.append(innerListHoverL)
    wLo.append(innerListLo)
    wLWP.append(innerListLWP)
    wCumuLWP.append(innerListCumuLWP)
    wAvgLWP.append(innerListAvgLWP)
    wWE.append(innerListWE)
    wCumuWE.append(innerListCumuWE)
    wAvgWE.append(innerListAvgWE)
    wT.append(innerListT)
    wStorm.append(innerListStorm)
    wCalm.append(innerListCalm)
    wRatio.append(innerListRatio)
    wNumStorms.append(innerListNumStorms)
    wIr.append(innerListIr)
    wAvgIr.append(innerListAvgIr)
    wStormHs.append(innerListStormHs)
    wStormCumuWE.append(innerListStormCumuWE)
    wStormCumuLWP.append(innerListStormCumuLWP)
    # swellHs.append(innerListSwellHs)
    # swellTp.append(innerListSwellTp)
    # swellDm.append(innerListSwellDm)
    # swellLWP.append(innerListSwellLWP)
    # windseaHs.append(innerListWindseaHs)
    # windseaTp.append(innerListWindseaTp)
    # windseaDm.append(innerListWindseaDm)
    # windseaLWP.append(innerListWindseaLWP)


wWaterLevel = []        # Offshore Wave Heights
wPredWaterLevel = []        # Offshore Wave Periods
wResWaterLevel = []        # Offshore Wave Directions
for xx in range(numClusters):
    innerListWaterLevel = []
    innerListPredWaterLevel = []
    innerListResWaterLevel = []
    for yy in range(numClusters):
        # if both are equal then the wave conditions didn't cause a transition and
        # everything between the last two dates should be added to a distribution
        # if yy is different...
        wInd = np.where(prevbin[xx]==yy)
        if len(wInd[0])>0:
            tempWaterLevel = []
            tempPredWaterLevel = []
            tempResWaterLevel = []
            for tt in range(len(wInd[0])):
                if daysBetween[xx][wInd[0][tt]] < 92:
                    tempT = np.where((tWaterLevelFRF < date[xx][wInd[0][tt]]) & (tWaterLevelFRF > prevdate[xx][wInd[0][tt]]))
                    if len(tempT[0]) > 0:
                        tempWaterLevel = np.append(tempWaterLevel,waterLevelFRF[tempT])
                        tempPredWaterLevel = np.append(tempPredWaterLevel,predictedWaterLevelFRF[tempT])
                        tempResWaterLevel = np.append(tempResWaterLevel,residualWaterLevelFRF[tempT])
        else:
            tempWaterLevel = []
            tempPredWaterLevel = []
            tempResWaterLevel = []
        innerListWaterLevel.append(tempWaterLevel)
        innerListPredWaterLevel.append(tempPredWaterLevel)
        innerListResWaterLevel.append(tempResWaterLevel)
    wWaterLevel.append(innerListWaterLevel)
    wPredWaterLevel.append(innerListPredWaterLevel)
    wResWaterLevel.append(innerListResWaterLevel)





wSusWindSpeed = []        # Offshore Wave Heights
wWindDirection = []        # Offshore Wave Periods
wWindSpeed = []        # Offshore Wave Directions
for xx in range(numClusters):
    innerListSusWindSpeed = []
    innerListWindDirection = []
    innerListWindSpeed = []
    for yy in range(numClusters):
        # if both are equal then the wave conditions didn't cause a transition and
        # everything between the last two dates should be added to a distribution
        # if yy is different...
        wInd = np.where(prevbin[xx]==yy)
        if len(wInd[0])>0:
            tempSusWindSpeed = []
            tempWindDirection = []
            tempWindSpeed = []
            for tt in range(len(wInd[0])):
                if daysBetween[xx][wInd[0][tt]] < 92:
                    tempT = np.where((tWinds < date[xx][wInd[0][tt]]) & (tWinds > prevdate[xx][wInd[0][tt]]))
                    if len(tempT[0]) > 0:
                        tempSusWindSpeed = np.append(tempSusWindSpeed,sustWindSpeed[tempT])
                        tempWindDirection  = np.append(tempWindDirection ,windDirection[tempT])
                        tempWindSpeed = np.append(tempWindSpeed,windSpeed[tempT])
        else:
            tempSusWindSpeed = []
            tempWindDirection  = []
            tempWindSpeed = []
        innerListSusWindSpeed.append(tempSusWindSpeed)
        innerListWindDirection .append(tempWindDirection)
        innerListWindSpeed.append(tempWindSpeed)
    wSusWindSpeed.append(innerListSusWindSpeed)
    wWindDirection .append(innerListWindDirection)
    wWindSpeed.append(innerListWindSpeed)


# can we separate into no transition, up transition, down transition


stayHs = []
stayTp = []
stayDm = []
stayLWP = []
downHs = []
downTp = []
downDm = []
downLWP = []
upHs = []
upTp = []
upDm = []
upLWP = []
#
# stayHsSwell = []
# stayTpSwell = []
# stayDmSwell = []
# stayLWPSwell = []
# downHsSwell = []
# downTpSwell = []
# downDmSwell = []
# downLWPSwell = []
# upHsSwell = []
# upTpSwell = []
# upDmSwell = []
# upLWPSwell = []
#
# stayHsWindsea = []
# stayTpWindsea = []
# stayDmWindsea = []
# stayLWPWindsea = []
# downHsWindsea = []
# downTpWindsea = []
# downDmWindsea = []
# downLWPWindsea = []
# upHsWindsea = []
# upTpWindsea = []
# upDmWindsea = []
# upLWPWindsea = []

for xx in range(numClusters):
    for yy in range(numClusters):
        if xx == yy:
            if len(wHs[xx][yy]) > 0:
                stayHs.append(wHs[xx][yy])
                stayTp.append(wTp[xx][yy])
                stayDm.append(wDm[xx][yy])
                stayLWP.append(wLWP[xx][yy])
                # stayHsSwell.append(swellHs[xx][yy])
                # stayTpSwell.append(swellTp[xx][yy])
                # stayDmSwell.append(swellDm[xx][yy])
                # stayLWPSwell.append(swellLWP[xx][yy])
                #
                # stayHsWindsea.append(windseaHs[xx][yy])
                # stayTpWindsea.append(windseaTp[xx][yy])
                # stayDmWindsea.append(windseaDm[xx][yy])
                # stayLWPWindsea.append(windseaLWP[xx][yy])
        if xx == (yy-1):
            if len(wHs[xx][yy]) > 0:
                downHs.append(wHs[xx][yy])
                downTp.append(wTp[xx][yy])
                downDm.append(wDm[xx][yy])
                downLWP.append(wLWP[xx][yy])
                # downHsSwell.append(swellHs[xx][yy])
                # downTpSwell.append(swellTp[xx][yy])
                # downDmSwell.append(swellDm[xx][yy])
                # downLWPSwell.append(swellLWP[xx][yy])
                #
                # downHsWindsea.append(windseaHs[xx][yy])
                # downTpWindsea.append(windseaTp[xx][yy])
                # downDmWindsea.append(windseaDm[xx][yy])
                # downLWPWindsea.append(windseaLWP[xx][yy])
        if xx == (yy+1):
            if len(wHs[xx][yy]) > 0:
                upHs.append(wHs[xx][yy])
                upTp.append(wTp[xx][yy])
                upDm.append(wDm[xx][yy])
                upLWP.append(wLWP[xx][yy])

        if xx == 14 and yy == 0:
            if len(wHs[xx][yy]) > 0:
                downHs.append(wHs[xx][yy])
                downTp.append(wTp[xx][yy])
                downDm.append(wDm[xx][yy])
                downLWP.append(wLWP[xx][yy])

        if xx == 0 and yy == 14:
            if len(wHs[xx][yy]) > 0:
                upHs.append(wHs[xx][yy])
                upTp.append(wTp[xx][yy])
                upDm.append(wDm[xx][yy])
                upLWP.append(wLWP[xx][yy])

                # upHsSwell.append(swellHs[xx][yy])
                # upTpSwell.append(swellTp[xx][yy])
                # upDmSwell.append(swellDm[xx][yy])
                # upLWPSwell.append(swellLWP[xx][yy])
                #
                # upHsWindsea.append(windseaHs[xx][yy])
                # upTpWindsea.append(windseaTp[xx][yy])
                # upDmWindsea.append(windseaDm[xx][yy])
                # upLWPWindsea.append(windseaLWP[xx][yy])

stayHsArray = np.concatenate(stayHs).ravel()
upHsArray = np.concatenate(upHs).ravel()
downHsArray = np.concatenate(downHs).ravel()
stayTpArray = np.concatenate(stayTp).ravel()
upTpArray = np.concatenate(upTp).ravel()
downTpArray = np.concatenate(downTp).ravel()
stayDmArray = np.concatenate(stayDm).ravel()
upDmArray = np.concatenate(upDm).ravel()
downDmArray = np.concatenate(downDm).ravel()
stayLWPArray = np.concatenate(stayLWP).ravel()
upLWPArray = np.concatenate(upLWP).ravel()
downLWPArray = np.concatenate(downLWP).ravel()

# stayHsSwellArray = np.concatenate(stayHsSwell).ravel()
# upHsSwellArray = np.concatenate(upHsSwell).ravel()
# downHsSwellArray = np.concatenate(downHsSwell).ravel()
# stayTpSwellArray = np.concatenate(stayTpSwell).ravel()
# upTpSwellArray = np.concatenate(upTpSwell).ravel()
# downTpSwellArray = np.concatenate(downTpSwell).ravel()
# stayDmSwellArray = np.concatenate(stayDmSwell).ravel()
# upDmSwellArray = np.concatenate(upDmSwell).ravel()
# downDmSwellArray = np.concatenate(downDmSwell).ravel()
# stayLWPSwellArray = np.concatenate(stayLWPSwell).ravel()
# upLWPSwellArray = np.concatenate(upLWPSwell).ravel()
# downLWPSwellArray = np.concatenate(downLWPSwell).ravel()
#
# stayHsWindseaArray = np.concatenate(stayHsWindsea).ravel()
# upHsWindseaArray = np.concatenate(upHsWindsea).ravel()
# downHsWindseaArray = np.concatenate(downHsWindsea).ravel()
# stayTpWindseaArray = np.concatenate(stayTpWindsea).ravel()
# upTpWindseaArray = np.concatenate(upTpWindsea).ravel()
# downTpWindseaArray = np.concatenate(downTpWindsea).ravel()
# stayDmWindseaArray = np.concatenate(stayDmWindsea).ravel()
# upDmWindseaArray = np.concatenate(upDmWindsea).ravel()
# downDmWindseaArray = np.concatenate(downDmWindsea).ravel()
# stayLWPWindseaArray = np.concatenate(stayLWPWindsea).ravel()
# upLWPWindseaArray = np.concatenate(upLWPWindsea).ravel()
# downLWPWindseaArray = np.concatenate(downLWPWindsea).ravel()

# Where are the extremes?

stayOver2m = np.where((stayHsArray>2.0))
percentStayOver2m = len(stayOver2m[0])/len(stayHsArray)
downOver2m = np.where((downHsArray>2.0))
percentDownOver2m = len(downOver2m[0])/len(downHsArray)
upOver2m = np.where((upHsArray>2.0))
percentUpOver2m = len(upOver2m[0])/len(upHsArray)



import matplotlib.cm as cm
import matplotlib.colors as mcolors
plt.style.use('default')

from scipy.stats.kde import gaussian_kde
dist_space = np.linspace(0, 4, 80)
fig = plt.figure(figsize=(10,10))
gs2 = gridspec.GridSpec(6, 5)

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
        normalize = mcolors.Normalize(vmin=0.6, vmax=1.8)

        ax.set_xlim([0,4])
        # ax.set_ylim([0,2])
        data = wIr[xx][yy]
        if len(data)>0:
            kde = gaussian_kde(data)
            colorparam[counter] = np.nanmean(data)
            colormap = cm.Reds
            color = colormap(normalize(colorparam[counter]))
            ax.plot(dist_space, kde(dist_space),linewidth=1,color=color)
            # ax.hist(data)
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
# cbar.set_label('Mean Hs (m)')
#cbar = ax.cax.colorbar(s_map, ticks=colorparam, format='%2.2g') # format='%2i' for integer





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
        normalize = mcolors.Normalize(vmin=0, vmax=1.5)

        ax.set_xlim([0,4])
        # ax.set_ylim([0,2])
        data = wNumStorms[xx][yy]
        if len(data)>0:
            #kde = gaussian_kde(data)
            #colorparam[counter] = np.nanmean(data)
            #colormap = cm.Reds
            #color = colormap(normalize(colorparam[counter]))
            #ax.plot(dist_space, kde(dist_space),linewidth=1,color=color)
            ax.hist(data,bins=[0,1,2,3,4])

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
# cbar.set_label('Mean Hs (m)')
#cbar = ax.cax.colorbar(s_map, ticks=colorparam, format='%2.2g') # format='%2i' for integer


plt.style.use('dark_background')
dist_space = np.linspace(-800, 800, 80)
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






plt.style.use('default')

import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.stats.kde import gaussian_kde
dist_space = np.linspace(-0.6, 0.6, 50)
fig = plt.figure(figsize=(10,10))
colorparam = np.zeros((numClusters*numClusters,))
counter = 0
for xx in range(numClusters):
    for yy in range(numClusters):
        ax = plt.subplot2grid((numClusters,numClusters), (yy,xx), rowspan=1, colspan=1)
        #normalize = mcolors.Normalize(vmin=np.min(colorparams), vmax=np.max(colorparams))
        normalize = mcolors.Normalize(vmin=-0.15, vmax=0.15)

        # ax.set_xlim([0,4])
        ax.set_ylim([0,5.5])
        data = wResWaterLevel[xx][yy]


        if len(data)>0:
            baddata = np.isnan(data)
            data = data[~baddata]+0.063
            kde = gaussian_kde(data)
            colorparam[counter] = np.nanmean(data)
            colormap = cm.coolwarm
            color = colormap(normalize(colorparam[counter]))
            ax.plot(dist_space, kde(dist_space),linewidth=1,color='k')
            y1 = kde(dist_space)
            negs = np.where((dist_space <=0))
            pos = np.where((dist_space >=0))
            ax.fill_between(dist_space[negs],0,y1[negs],color='blue')
            ax.fill_between(dist_space[pos],0,y1[pos],color='red')

            # ax.hist(data)
            ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
            ax.spines['top'].set_color([0.5, 0.5, 0.5])
            ax.spines['right'].set_color([0.5, 0.5, 0.5])
            ax.spines['left'].set_color([0.5, 0.5, 0.5])
            #ax.text(1.8, 1, np.round(colorparam*100)/100, fontweight='bold')

        else:
            ax.spines['bottom'].set_color([0.75, 0.75, 0.75])
            ax.spines['top'].set_color([0.75, 0.75, 0.75])
            ax.spines['right'].set_color([0.75, 0.75, 0.75])
            ax.spines['left'].set_color([0.75, 0.75, 0.75])
            if yy < 8:
                ax.xaxis.set_ticks([])
        if yy < 7:
            ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
        counter = counter+1
plt.show()
plt.tight_layout()
# s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
# s_map.set_array(colorparam)
# fig.subplots_adjust(right=0.92)
# cbar_ax = fig.add_axes([0.94, 0.10, 0.02, 0.7])
# cbar = fig.colorbar(s_map, cax=cbar_ax)
# cbar.set_label('Mean Hs (m)')
#cbar = ax.cax.colorbar(s_map, ticks=colorparam, format='%2.2g') # format='%2i' for integer











plt.style.use('default')
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, numClusters))
colors2 = cm.gray(np.linspace(0, 1, numClusters))

## Chord option #1

# chord diagram
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

import numpy as np

LW = 0.3

def polar2xy(r, theta):
    return np.array([r*np.cos(theta), r*np.sin(theta)])

def hex2rgb(c):
    return tuple(int(c[i:i+2], 16)/256.0 for i in (1, 3 ,5))

def IdeogramArc(start=0, end=60, radius=1.0, width=0.2, ax=None, color=(1,0,0)):
    # start, end should be in [0, 360)
    if start > end:
        start, end = end, start
    start *= np.pi/180.
    end *= np.pi/180.
    # optimal distance to the control points
    # https://stackoverflow.com/questions/1734745/how-to-create-circle-with-b%C3%A9zier-curves
    opt = 4./3. * np.tan((end-start)/ 4.) * radius
    inner = radius*(1-width)
    verts = [
        polar2xy(radius, start),
        polar2xy(radius, start) + polar2xy(opt, start+0.5*np.pi),
        polar2xy(radius, end) + polar2xy(opt, end-0.5*np.pi),
        polar2xy(radius, end),
        polar2xy(inner, end),
        polar2xy(inner, end) + polar2xy(opt*(1-width), end-0.5*np.pi),
        polar2xy(inner, start) + polar2xy(opt*(1-width), start+0.5*np.pi),
        polar2xy(inner, start),
        polar2xy(radius, start),
        ]

    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.LINETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CLOSEPOLY,
             ]

    if ax == None:
        return verts, codes
    else:
        path = Path(verts, codes)
        print(color)
        patch = patches.PathPatch(path, facecolor=color, edgecolor=color, lw=LW, alpha=0.5)
        #patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.4,), lw=LW)

        ax.add_patch(patch)


def ChordArc(start1=0, end1=60, start2=180, end2=240, radius=1.0, chordwidth=0.7, ax=None, color=(1,0,0)):
    # start, end should be in [0, 360)
    if start1 > end1:
        start1, end1 = end1, start1
    if start2 > end2:
        start2, end2 = end2, start2
    start1 *= np.pi/180.
    end1 *= np.pi/180.
    start2 *= np.pi/180.
    end2 *= np.pi/180.
    opt1 = 4./3. * np.tan((end1-start1)/ 4.) * radius
    opt2 = 4./3. * np.tan((end2-start2)/ 4.) * radius
    rchord = radius * (1-chordwidth)
    verts = [
        polar2xy(radius, start1),
        polar2xy(radius, start1) + polar2xy(opt1, start1+0.5*np.pi),
        polar2xy(radius, end1) + polar2xy(opt1, end1-0.5*np.pi),
        polar2xy(radius, end1),
        polar2xy(rchord, end1),
        polar2xy(rchord, start2),
        polar2xy(radius, start2),
        polar2xy(radius, start2) + polar2xy(opt2, start2+0.5*np.pi),
        polar2xy(radius, end2) + polar2xy(opt2, end2-0.5*np.pi),
        polar2xy(radius, end2),
        polar2xy(rchord, end2),
        polar2xy(rchord, start1),
        polar2xy(radius, start1),
        ]

    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             ]

    if ax == None:
        return verts, codes
    else:
        path = Path(verts, codes)
        #patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.4,), lw=LW)
        patch = patches.PathPatch(path, facecolor=color, edgecolor=color, lw=LW, alpha=0.5)

        ax.add_patch(patch)

def selfChordArc(start=0, end=60, radius=1.0, chordwidth=0.7, ax=None, color=(1,0,0)):
    # start, end should be in [0, 360)
    if start > end:
        start, end = end, start
    start *= np.pi/180.
    end *= np.pi/180.
    opt = 4./3. * np.tan((end-start)/ 4.) * radius
    rchord = radius * (1-chordwidth)
    verts = [
        polar2xy(radius, start),
        polar2xy(radius, start) + polar2xy(opt, start+0.5*np.pi),
        polar2xy(radius, end) + polar2xy(opt, end-0.5*np.pi),
        polar2xy(radius, end),
        polar2xy(rchord, end),
        polar2xy(rchord, start),
        polar2xy(radius, start),
        ]

    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             ]

    if ax == None:
        return verts, codes
    else:
        path = Path(verts, codes)
        #patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.4,), lw=LW)
        patch = patches.PathPatch(path, facecolor=color, edgecolor=color, lw=LW, alpha=0.5)

        ax.add_patch(patch)

def chordDiagram(X, ax, colors=None, width=0.1, pad=2, chordwidth=0.7):
    """Plot a chord diagram
    Parameters
    ----------
    X :
        flux data, X[i, j] is the flux from i to j
    ax :
        matplotlib `axes` to show the plot
    colors : optional
        user defined colors in rgb format. Use function hex2rgb() to convert hex color to rgb color. Default: d3.js category10
    width : optional
        width/thickness of the ideogram arc
    pad : optional
        gap pad between two neighboring ideogram arcs, unit: degree, default: 2 degree
    chordwidth : optional
        position of the control points for the chords, controlling the shape of the chords
    """
    # X[i, j]:  i -> j
    x = X.sum(axis = 1) # sum over rows
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)

    if colors is None:
    # use d3.js category10 https://github.com/d3/d3-3.x-api-reference/blob/master/Ordinal-Scales.md#category10
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                  '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        if len(x) > 10:
            print('x is too large! Use x smaller than 10')
        colors = [hex2rgb(colors[i]) for i in range(len(x))]

    # find position for each start and end
    y = x/np.sum(x).astype(float) * (360 - pad*len(x))

    pos = {}
    arc = []
    nodePos = []
    start = 0
    for i in range(len(x)):
        end = start + y[i]
        arc.append((start, end))
        angle = 0.5*(start+end)
        #print(start, end, angle)
        if -30 <= angle <= 210:
            angle -= 90
        else:
            angle -= 270
        nodePos.append(tuple(polar2xy(1.1, 0.5*(start+end)*np.pi/180.)) + (angle,))
        z = (X[i, :]/x[i].astype(float)) * (end - start)
        ids = np.argsort(z)
        z0 = start
        for j in ids:
            pos[(i, j)] = (z0, z0+z[j])
            z0 += z[j]
        start = end + pad

    for i in range(len(x)):
        start, end = arc[i]
        print(colors[i])
        IdeogramArc(start=start, end=end, radius=1.0, ax=ax, color=colors[i], width=width)
        start, end = pos[(i,i)]
        selfChordArc(start, end, radius=1.-width, color=colors[i], chordwidth=chordwidth*0.7, ax=ax)
        for j in range(i):
            color = colors[i]
            if X[i, j] > X[j, i]:
                color = colors[j]
            start1, end1 = pos[(i,j)]
            start2, end2 = pos[(j,i)]
            ChordArc(start1, end1, start2, end2,
                     radius=1.-width, color=colors[i], chordwidth=chordwidth, ax=ax)

    #print(nodePos)
    return nodePos



##################################
#if __name__ == "__main__":
fig = plt.figure(figsize=(6,6))
#flux = matrix
flux = flatarray2
ax = plt.axes([0,0,1,1])

    #nodePos = chordDiagram(flux, ax, colors=[hex2rgb(x) for x in ['#666666', '#66ff66', '#ff6666', '#6666ff']])
#nodePos = chordDiagram(flux, ax)
tempColorsInd = np.where((colors==1))
tempColorsInd2 = np.where((colors==0))

tempColors = colors
tempColors[tempColorsInd] = 0.9999
tempColors[tempColorsInd2] = 0.0001
tempColors = tempColors
nodePos = chordDiagram(flux,ax,colors=tempColors[:,0:3])
ax.axis('off')
prop = dict(fontsize=16*0.8, ha='center', va='center')
#plt.title('Winter')
plt.show()




########### Lets Spit out some variables for random forest generation ########
## We want bmus, time diff,
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
    lwpList.append(lwpC[tempWave])
    weList.append(weC[tempWave])
    fvList.append(fV[tempWave])
    irList.append(Ir[tempWave])
    hlList.append(HoverL[tempWave])

dataPickle = 'sandbarsSouthernTransectRandomForestInputData8.pickle'
output = {}

output['currentBmu'] = wCurrentBmu
output['prevBmu'] = wPrevBmu
output['daysBetween'] = daysBetween
output['avgHs'] = wAvgHs
output['avgTp'] = wAvgTp
output['avgWE'] = wAvgWE
output['avgLWP'] = wAvgLWP
output['avgIR'] = wAvgIr
output['stormCumuLWP'] = wStormCumuLWP
output['stormCumuWE'] = wStormCumuWE
output['stormAvgHs'] = wStormHs
output['days'] = wDays
output['sortedBmus'] = sorted_bmus
output['wWindSpeed'] = wWindSpeed
output['wWindDirection'] = wWindDirection
output['wSustWindSpeed'] = wSusWindSpeed
output['wHs'] = wHs
output['wTp'] = wTp
output['wDm'] = wDm
output['wResWaterLevel'] = wResWaterLevel
output['wLWP'] = wLWP
output['wWE'] = wWE
output['wIr'] = wIr
output['wHoverL'] = wHoverL
output['wPredWaterLevel'] = wPredWaterLevel
output['wWaterLevel'] = wWaterLevel

output['hsSplit'] = hsSplit
output['tpSplit'] = tpSplit
output['dmSplit'] = dmSplit
output['lwpSplit'] = lwpSplit
output['weSplit'] = weSplit
output['hoverlSplit'] = hoverlSplit
output['loSplit'] = loSplit
output['irSplit'] = irSplit
output['fvSplit'] = fvSplit

output['hsList'] = hsList
output['tpList'] = tpList
output['dmList'] = dmList
output['lwpList'] = lwpList
output['weList'] = weList
output['fvList'] = fvList
output['irList'] = irList
output['hlList'] = hlList

output['waterLevelSplit'] = waterLevelSplit
output['predWaterLevelSplit'] = predWaterLevelSplit
output['resWaterLevelSplit'] = resWaterLevelSplit

output['sustWindSpeedSplit'] = sustWindSpeedSplit
output['windSpeedSplit'] = windSpeedSplit
output['windDirectionSplit'] = windDirectionSplit
output['sorted_bmus'] = sorted_bmus



import pickle
with open(dataPickle,'wb') as f:
    pickle.dump(output, f)

# dataPickle = 'sandbarsSouthernTransectRandomForestInputData.pickle'
# output = {}
# output['currentBmu'] = wCurrentBmu
# output['prevBmu'] = wPrevBmu
# output['daysBetween'] = daysBetween
# output['avgHs'] = wAvgHs
# output['avgTp'] = wAvgTp
# output['avgWE'] = wAvgWE
# output['avgLWP'] = wAvgLWP
# output['avgIR'] = wAvgIr
# output['stormCumuLWP'] = wStormCumuLWP
# output['stormCumuWE'] = wStormCumuWE
# output['stormAvgHs'] = wStormHs
# output['days'] = wDays
# output['sortedBmus'] = sorted_bmus
#
# import pickle
# with open(dataPickle,'wb') as f:
#     pickle.dump(output, f)





colors = cm.hsv(np.linspace(0, 1, (numClusters+1)))

fig13 = plt.figure(figsize=(12,4))

ax1 = plt.subplot2grid((2,5),(0,0),rowspan=1,colspan=1)
ii = 6
finder = order[ii]
ax1.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)
# ax1.set_yticks([])
# ax1.set_xticks([])
ax1.set_xlim([0, 500])
ax1.set_ylim([-8, 0.5])


ax2 = plt.subplot2grid((2,5),(0,1),rowspan=1,colspan=1)
ii = 6
finder = order[ii]
ax2.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)
ax2.set_yticks([])
ax2.set_xticks([])
ax2.set_xlim([0, 500])
ax2.set_ylim([-8, 0.5])
ax3 = plt.subplot2grid((2,5),(0,2),rowspan=1,colspan=1)
ii = 7
finder = order[ii]
ax3.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)
ax3.set_yticks([])
ax3.set_xticks([])
ax3.set_xlim([0, 500])
ax3.set_ylim([-8, 0.5])
ax4 = plt.subplot2grid((2,5),(0,3),rowspan=1,colspan=1)
ii = 0
finder = order[ii]
ax4.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)
ax4.set_yticks([])
ax4.set_xticks([])
ax4.set_xlim([0, 500])
ax4.set_ylim([-8, 0.5])
ax5 = plt.subplot2grid((2,5),(0,4),rowspan=1,colspan=1)
ii = 1
finder = order[ii]
ax5.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)
ax5.set_yticks([])
ax5.set_xticks([])
ax5.set_xlim([0, 500])
ax5.set_ylim([-8, 0.5])

counter = 0
#dist_space = np.arange(0,5,0.1)
dist_space = np.linspace(0, 5, 50)

normalize = mcolors.Normalize(vmin=.5, vmax=1.5)
colorparam = np.empty((4,))
ax6 = plt.subplot2grid((2, 5), (1, 1), rowspan=1, colspan=1)
ax6.set_xlim([0, 3])
ax6.set_ylim([0, 1.5])
xx = 6
yy = 7
data = wHs[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(data)
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax6.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax6.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax6.spines['top'].set_color([0.5, 0.5, 0.5])
ax6.spines['right'].set_color([0.5, 0.5, 0.5])
ax6.spines['left'].set_color([0.5, 0.5, 0.5])
ax6.set_yticks([])
ax6.set_xlabel('Hs (m)')

counter = counter+1
ax7 = plt.subplot2grid((2, 5), (1, 2), rowspan=1, colspan=1)
ax7.set_xlim([0, 3])
ax7.set_ylim([0, 1.5])
xx = 7
yy = 7
data = wHs[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(data)
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax7.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax7.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax7.spines['top'].set_color([0.5, 0.5, 0.5])
ax7.spines['right'].set_color([0.5, 0.5, 0.5])
ax7.spines['left'].set_color([0.5, 0.5, 0.5])
ax7.set_yticks([])
ax7.set_xlabel('Hs (m)')

counter = counter+1

ax8 = plt.subplot2grid((2, 5), (1, 3), rowspan=1, colspan=1)
ax8.set_xlim([0, 3])
ax8.set_ylim([0, 1.5])
xx = 0
yy = 7
data = wHs[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(data)
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax8.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax8.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax8.spines['top'].set_color([0.5, 0.5, 0.5])
ax8.spines['right'].set_color([0.5, 0.5, 0.5])
ax8.spines['left'].set_color([0.5, 0.5, 0.5])
ax8.set_yticks([])
ax8.set_xlabel('Hs (m)')

counter = counter+1

ax9 = plt.subplot2grid((2, 5), (1, 4), rowspan=1, colspan=1)
ax9.set_xlim([0, 3])
ax9.set_ylim([0, 1.5])
xx = 1
yy = 7
data = wHs[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(data)
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax9.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax9.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax9.spines['top'].set_color([0.5, 0.5, 0.5])
ax9.spines['right'].set_color([0.5, 0.5, 0.5])
ax9.spines['left'].set_color([0.5, 0.5, 0.5])
ax9.set_yticks([])
ax9.set_xlabel('Hs (m)')

colorparam1 = colorparam
normalize1 = normalize












#########



colors = cm.hsv(np.linspace(0, 1, (numClusters+1)))

n = 6
m = 4

fig13 = plt.figure(figsize=(12,10))

ax1 = plt.subplot2grid((n,m),(0,1),rowspan=1,colspan=1)
ii = 7
finder = order[ii]
ax1.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)

ax1.set_xlim([0, 500])
ax1.set_ylim([-8, 0.5])
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['top'].set_linewidth(2)
ax1.spines['right'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)

ax2 = plt.subplot2grid((n,m),(1,0),rowspan=1,colspan=1)
ii = 6
finder = order[ii]
ax2.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)
ax2.set_yticks([])
ax2.set_xticks([])
ax2.set_xlim([0, 500])
ax2.set_ylim([-8, 0.5])
ax3 = plt.subplot2grid((n,m),(1,1),rowspan=1,colspan=1)
ii = 7
finder = order[ii]
ax3.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)
ax3.set_yticks([])
ax3.set_xticks([])
ax3.set_xlim([0, 500])
ax3.set_ylim([-8, 0.5])
ax4 = plt.subplot2grid((n,m),(1,2),rowspan=1,colspan=1)
ii = 0
finder = order[ii]
ax4.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)
ax4.set_yticks([])
ax4.set_xticks([])
ax4.set_xlim([0, 500])
ax4.set_ylim([-8, 0.5])
ax5 = plt.subplot2grid((n,m),(1,3),rowspan=1,colspan=1)
ii = 1
finder = order[ii]
ax5.plot(xinterp, profiles[finder, :], color=colors[ii, :], label=ii)
ax5.set_yticks([])
ax5.set_xticks([])
ax5.set_xlim([0, 500])
ax5.set_ylim([-8, 0.5])

counter = 0
#dist_space = np.arange(0,5,0.1)
dist_space = np.linspace(0, 5, 50)

normalize = mcolors.Normalize(vmin=.5, vmax=1.5)
colorparam = np.empty((4,))
ax6 = plt.subplot2grid((n, m), (2, 0), rowspan=1, colspan=1)
ax6.set_xlim([0, 3])
ax6.set_ylim([0, 1.5])
xx = 6
yy = 7
data = wHs[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(data)
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax6.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax6.text(1.5,1.2,r'$Hs_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax6.text(1.5,0.9,r'$Hs_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax6.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax6.spines['top'].set_color([0.5, 0.5, 0.5])
ax6.spines['right'].set_color([0.5, 0.5, 0.5])
ax6.spines['left'].set_color([0.5, 0.5, 0.5])
ax6.set_yticks([])
ax6.set_xlabel('Hs (m)')

counter = counter+1
ax7 = plt.subplot2grid((n, m), (2, 1), rowspan=1, colspan=1)
ax7.set_xlim([0, 3])
ax7.set_ylim([0, 1.5])
xx = 7
yy = 7
data = wHs[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(data)
ax7.text(1.5,1.2,r'$Hs_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax7.text(1.5,0.9,r'$Hs_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax7.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax7.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax7.spines['top'].set_color([0.5, 0.5, 0.5])
ax7.spines['right'].set_color([0.5, 0.5, 0.5])
ax7.spines['left'].set_color([0.5, 0.5, 0.5])
ax7.set_yticks([])
ax7.set_xlabel('Hs (m)')

counter = counter+1

ax8 = plt.subplot2grid((n, m), (2, 2), rowspan=1, colspan=1)
ax8.set_xlim([0, 3])
ax8.set_ylim([0, 1.5])
xx = 0
yy = 7
data = wHs[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(data)
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax8.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax8.text(1.5,1.2,r'$Hs_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax8.text(1.5,0.9,r'$Hs_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax8.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax8.spines['top'].set_color([0.5, 0.5, 0.5])
ax8.spines['right'].set_color([0.5, 0.5, 0.5])
ax8.spines['left'].set_color([0.5, 0.5, 0.5])
ax8.set_yticks([])
ax8.set_xlabel('Hs (m)')

counter = counter+1

ax9 = plt.subplot2grid((n, m), (2, 3), rowspan=1, colspan=1)
ax9.set_xlim([0, 3])
ax9.set_ylim([0, 1.5])
xx = 1
yy = 7
data = wHs[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(data)
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax9.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax9.text(1.5,1.2,r'$Hs_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax9.text(1.5,0.9,r'$Hs_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax9.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax9.spines['top'].set_color([0.5, 0.5, 0.5])
ax9.spines['right'].set_color([0.5, 0.5, 0.5])
ax9.spines['left'].set_color([0.5, 0.5, 0.5])
ax9.set_yticks([])
ax9.set_xlabel('Hs (m)')

normalize1 = normalize
colorparam1 = colorparam
# s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
# s_map.set_array(colorparam)
# fig13.subplots_adjust(right=0.89)
# cbar_ax = fig13.add_axes([0.92, 0.6, 0.02, 0.2])
# cbar = fig13.colorbar(s_map, cax=cbar_ax)
# cbar.set_label('Mean Hs (m)')
# #cbar = ax.cax.colorbar(s_map, ticks=colorparam, format='%2.2g') # format='%2i' for integer








counter = 0
#dist_space = np.arange(0,5,0.1)
dist_space = np.linspace(-800, 800, 50)

normalize = mcolors.Normalize(vmin=0, vmax=200)
colorparam = np.empty((4,))
ax6 = plt.subplot2grid((n, m), (3, 0), rowspan=1, colspan=1)
ax6.set_xlim([0, 500])
ax6.set_ylim([0, 0.008])
xx = 6
yy = 7
data = np.abs(wLWP[xx][yy])
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(np.abs(data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax6.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax6.text(200.0,0.007,r'$LWP_{50\%}$' + ' = {:.0f}'.format(np.nanmean(data)))
ax6.text(200.0,0.005,r'$LWP_{90\%}$' + ' = {:.0f}'.format(np.nanpercentile(data,90)))
ax6.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax6.spines['top'].set_color([0.5, 0.5, 0.5])
ax6.spines['right'].set_color([0.5, 0.5, 0.5])
ax6.spines['left'].set_color([0.5, 0.5, 0.5])
ax6.set_yticks([])
ax6.set_xlabel('LWP')

counter = counter+1
ax7 = plt.subplot2grid((n, m), (3, 1), rowspan=1, colspan=1)
ax7.set_xlim([0, 500])
ax7.set_ylim([0, 0.008])
xx = 7
yy = 7
data = np.abs(wLWP[xx][yy])
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(np.abs(data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax7.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax7.text(200.0,0.007,r'$LWP_{50\%}$' + ' = {:.0f}'.format(np.nanmean(data)))
ax7.text(200.0,0.005,r'$LWP_{90\%}$' + ' = {:.0f}'.format(np.nanpercentile(data,90)))
ax7.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax7.spines['top'].set_color([0.5, 0.5, 0.5])
ax7.spines['right'].set_color([0.5, 0.5, 0.5])
ax7.spines['left'].set_color([0.5, 0.5, 0.5])
ax7.set_yticks([])
ax7.set_xlabel('LWP')

counter = counter+1

ax8 = plt.subplot2grid((n, m), (3, 2), rowspan=1, colspan=1)
ax8.set_xlim([0, 500])
ax8.set_ylim([0, 0.008])
xx = 0
yy = 7
data = np.abs(wLWP[xx][yy])
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(np.abs(data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax8.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax8.text(200.0,0.007,r'$LWP_{50\%}$' + ' = {:.0f}'.format(np.nanmean(data)))
ax8.text(200.0,0.005,r'$LWP_{90\%}$' + ' = {:.0f}'.format(np.nanpercentile(data,90)))
ax8.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax8.spines['top'].set_color([0.5, 0.5, 0.5])
ax8.spines['right'].set_color([0.5, 0.5, 0.5])
ax8.spines['left'].set_color([0.5, 0.5, 0.5])
ax8.set_yticks([])
ax8.set_xlabel('LWP')

counter = counter+1

ax9 = plt.subplot2grid((n, m), (3,3), rowspan=1, colspan=1)
ax9.set_xlim([0, 500])
ax9.set_ylim([0, 0.008])
xx = 1
yy = 7
data = np.abs(wLWP[xx][yy])
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(np.abs(data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax9.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax9.text(200.0,0.007,r'$LWP_{50\%}$' + ' = {:.0f}'.format(np.nanmean(data)))
ax9.text(200.0,0.005,r'$LWP_{90\%}$' + ' = {:.0f}'.format(np.nanpercentile(data,90)))
ax9.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax9.spines['top'].set_color([0.5, 0.5, 0.5])
ax9.spines['right'].set_color([0.5, 0.5, 0.5])
ax9.spines['left'].set_color([0.5, 0.5, 0.5])
ax9.set_yticks([])
ax9.set_xlabel('LWP')

colorparam2 = colorparam
normalize2 = normalize






counter = 0
#dist_space = np.arange(0,5,0.1)
dist_space = np.linspace(0, 4, 50)

normalize = mcolors.Normalize(vmin=0, vmax=3)
colorparam = np.empty((4,))
ax6 = plt.subplot2grid((n, m), (4, 0), rowspan=1, colspan=1)
ax6.set_xlim([0, 4])
ax6.set_ylim([0, 1.2])
xx = 6
yy = 7
data = wIr[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(np.abs(data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax6.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax6.text(2.,0.95,r'$Ir_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax6.text(2.,0.75,r'$Ir_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax6.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax6.spines['top'].set_color([0.5, 0.5, 0.5])
ax6.spines['right'].set_color([0.5, 0.5, 0.5])
ax6.spines['left'].set_color([0.5, 0.5, 0.5])
# ax6.set_yticks([])
ax6.set_xlabel('Ir')

counter = counter+1
ax7 = plt.subplot2grid((n, m), (4, 1), rowspan=1, colspan=1)
ax7.set_xlim([0, 4])
ax7.set_ylim([0, 1.2])
xx = 7
yy = 7
data = wIr[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(np.abs(data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax7.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax7.text(2.,0.95,r'$Ir_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax7.text(2.,0.75,r'$Ir_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax7.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax7.spines['top'].set_color([0.5, 0.5, 0.5])
ax7.spines['right'].set_color([0.5, 0.5, 0.5])
ax7.spines['left'].set_color([0.5, 0.5, 0.5])
ax7.set_yticks([])
ax7.set_xlabel('Ir')

counter = counter+1

ax8 = plt.subplot2grid((n, m), (4, 2), rowspan=1, colspan=1)
ax8.set_xlim([0, 4])
ax8.set_ylim([0, 1.2])
xx = 0
yy = 7
data = wIr[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(np.abs(data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax8.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax8.text(2.,0.95,r'$Ir_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax8.text(2.,0.75,r'$Ir_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax8.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax8.spines['top'].set_color([0.5, 0.5, 0.5])
ax8.spines['right'].set_color([0.5, 0.5, 0.5])
ax8.spines['left'].set_color([0.5, 0.5, 0.5])
ax8.set_yticks([])
ax8.set_xlabel('Ir')

counter = counter+1

ax9 = plt.subplot2grid((n, m), (4, 3), rowspan=1, colspan=1)
ax9.set_xlim([0, 4])
ax9.set_ylim([0, 1.2])
xx = 1
yy = 7
data = wIr[xx][yy]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean(np.abs(data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax9.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax9.text(2.,0.95,r'$Ir_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax9.text(2.,0.75,r'$Ir_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax9.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax9.spines['top'].set_color([0.5, 0.5, 0.5])
ax9.spines['right'].set_color([0.5, 0.5, 0.5])
ax9.spines['left'].set_color([0.5, 0.5, 0.5])
ax9.set_yticks([])
ax9.set_xlabel('Ir')



counter = 0
#dist_space = np.arange(0,5,0.1)
dist_space = np.linspace(-0.75, 0.75, 50)

normalize = mcolors.Normalize(vmin=0, vmax=0.3)
colorparam = np.empty((4,))
ax6 = plt.subplot2grid((n, m), (5, 0), rowspan=1, colspan=1)
ax6.set_xlim([-0.75, 0.75])
ax6.set_ylim([0, 3.2])
xx = 6
yy = 7
data = wResWaterLevel[xx][yy]
baddata = np.isnan(data)
data = data[~baddata]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean((data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax6.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax6.text(-0.7,2.5,r'$NTR_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax6.text(-0.7,1.75,r'$NTR_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax6.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax6.spines['top'].set_color([0.5, 0.5, 0.5])
ax6.spines['right'].set_color([0.5, 0.5, 0.5])
ax6.spines['left'].set_color([0.5, 0.5, 0.5])
ax6.set_yticks([])
ax6.set_xlabel('NTR (m)')

counter = counter+1
ax7 = plt.subplot2grid((n, m), (5, 1), rowspan=1, colspan=1)
ax7.set_xlim([-0.75, 0.75])
ax7.set_ylim([0, 3.2])
xx = 7
yy = 7
data = wResWaterLevel[xx][yy]
baddata = np.isnan(data)
data = data[~baddata]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean((data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax7.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax7.text(-0.7,2.5,r'$NTR_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax7.text(-0.7,1.75,r'$NTR_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax7.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax7.spines['top'].set_color([0.5, 0.5, 0.5])
ax7.spines['right'].set_color([0.5, 0.5, 0.5])
ax7.spines['left'].set_color([0.5, 0.5, 0.5])
ax7.set_yticks([])
ax7.set_xlabel('NTR (m)')

counter = counter+1

ax8 = plt.subplot2grid((n, m), (5, 2), rowspan=1, colspan=1)
ax8.set_xlim([-0.75, 0.75])
ax8.set_ylim([0, 3.2])
xx = 0
yy = 7
data = wResWaterLevel[xx][yy]
baddata = np.isnan(data)
data = data[~baddata]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean((data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax8.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax8.text(-0.7,2.5,r'$NTR_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax8.text(-0.7,1.75,r'$NTR_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax8.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax8.spines['top'].set_color([0.5, 0.5, 0.5])
ax8.spines['right'].set_color([0.5, 0.5, 0.5])
ax8.spines['left'].set_color([0.5, 0.5, 0.5])
ax8.set_yticks([])
ax8.set_xlabel('NTR (m)')

counter = counter+1

ax9 = plt.subplot2grid((n, m), (5, 3), rowspan=1, colspan=1)
ax9.set_xlim([-0.75, 0.75])
ax9.set_ylim([0, 3.2])
xx = 1
yy = 7
data = wResWaterLevel[xx][yy]
baddata = np.isnan(data)
data = data[~baddata]
kde = gaussian_kde(data)
colorparam[counter] = np.nanmean((data))
colormap = cm.Reds
color = colormap(normalize(colorparam[counter]))
ax9.plot(dist_space, kde(dist_space), linewidth=1, color=color)
ax9.text(-0.7,2.5,r'$NTR_{50\%}$' + ' = {:.2f}'.format(np.nanmean(data)))
ax9.text(-0.7,1.75,r'$NTR_{90\%}$' + ' = {:.2f}'.format(np.nanpercentile(data,90)))
ax9.spines['bottom'].set_color([0.5, 0.5, 0.5])
ax9.spines['top'].set_color([0.5, 0.5, 0.5])
ax9.spines['right'].set_color([0.5, 0.5, 0.5])
ax9.spines['left'].set_color([0.5, 0.5, 0.5])
ax9.set_yticks([])
ax9.set_xlabel('NTR (m)')


plt.tight_layout()
#
# fig13.subplots_adjust(right=0.89)
#
# s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
# s_map.set_array(colorparam)
# cbar_ax = fig13.add_axes([0.92, 0.05, 0.02, 0.2])
# cbar = fig13.colorbar(s_map, cax=cbar_ax)
# cbar.set_label('Mean Ir')
#
#
# s_map = cm.ScalarMappable(norm=normalize2, cmap=colormap)
# s_map.set_array(colorparam2)
# cbar_ax = fig13.add_axes([0.92, 0.3, 0.02, 0.2])
# cbar = fig13.colorbar(s_map, cax=cbar_ax)
# cbar.set_label('Mean abs(LWP) (W)')
# #cbar = ax.cax.colorbar(s_map, ticks=colorparam, format='%2.2g') # format='%2i' for integer
#
#
#
# s_map = cm.ScalarMappable(norm=normalize1, cmap=colormap)
# s_map.set_array(colorparam1)
# cbar_ax = fig13.add_axes([0.92, 0.54, 0.02, 0.2])
# cbar = fig13.colorbar(s_map, cax=cbar_ax)
# cbar.set_label('Mean Hs (m)')
#cbar = ax.cax.colorbar(s_map, ticks=colorparam, format='%2.2g') # format='%2i' for integer


#cbar = ax.cax.colorbar(s_map, ticks=colorparam, format='%2.2g') # format='%2i' for integer




