
# Ok, first lets load a DWT struct

import scipy.io as sio
from scipy.io.matlab.mio5_params import mat_struct
import numpy as np
import more_itertools as mit
import os
import datetime
from netCDF4 import Dataset
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib import gridspec
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
import mda

from scipy.stats.kde import gaussian_kde
import datetime as DT

waterLevelDir = '/media/dylananderson/Elements1/frfWaterLevel'
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


wavedir = '/media/dylananderson/Elements/WIS_ST63218/'

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
    #timeTemp = [datenum_to_datetime(x) for x in waves['t'].flatten()]
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



wavedir26 = '/media/dylananderson/Elements1/26mArray/'
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
    Hs26m = np.append(Hs26m,waves26m['waveHs'][::2])
    Tp26m = np.append(Tp26m,waves26m['waveTp'][::2])
    Dm26m = np.append(Dm26m,waves26m['waveMeanDirection'][::2])
    timeWave26m = np.append(timeWave26m,waves26m['t'][::2])

ind = np.where((Hs26m > 0))
hs26m = Hs26m[ind]
tp26m = Tp26m[ind]
dm26m = Dm26m[ind]
t26m = timeWave26m[ind]
import datetime as DT
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
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


hsCombined = np.append(Hs,hs26m)
hsSmooth = moving_average(hsCombined,3)
tpCombined = np.append(Tp,tp26m)
dmCombined = np.append(Dm,dm26m)
tWave = [DT.datetime.fromtimestamp(x) for x in timeWave]
# tWave = [datetime.fromtimestamp(x) for x in timeWave]
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

badDirs = np.where((dmCombined > 360))
dmCombined[badDirs] = dmCombined[badDirs]*np.nan

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

# waveNormSwell = dmSwellCombined - 72
# negSwell = np.where((waveNormSwell > 180))
# waveNormSwell[negSwell[0]] = waveNormSwell[negSwell[0]]-360
# offposSwell = np.where((waveNormSwell>90))
# offnegSwell = np.where((waveNormSwell<-90))
# waveNormSwell[offposSwell[0]] = waveNormSwell[offposSwell[0]]*0
# waveNormSwell[offnegSwell[0]] = waveNormSwell[offnegSwell[0]]*0

# waveNormWindsea = dmWindseaCombined - 72
# negWindsea = np.where((waveNormWindsea > 180))
# waveNormWindsea[negWindsea[0]] = waveNormWindsea[negWindsea[0]]-360
# offposWindsea = np.where((waveNormWindsea>90))
# offnegWindsea = np.where((waveNormWindsea<-90))
# waveNormWindsea[offposWindsea[0]] = waveNormWindsea[offposWindsea[0]]*0
# waveNormWindsea[offnegWindsea[0]] = waveNormWindsea[offnegWindsea[0]]*0

lwpC = 1025*np.square(hsCombined)*tpCombined*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
weC = np.square(hsCombined)*tpCombined

# lwpSwell = 1025*np.square(hsSwellCombined)*tpSwellCombined*(9.81/(64*np.pi))*np.cos(waveNormSwell*(np.pi/180))*np.sin(waveNormSwell*(np.pi/180))
# weSwell = np.square(hsSwellCombined)*tpSwellCombined

# lwpWindsea = 1025*np.square(hsWindseaCombined)*tpWindseaCombined*(9.81/(64*np.pi))*np.cos(waveNormWindsea*(np.pi/180))*np.sin(waveNormWindsea*(np.pi/180))
# weWindsea = np.square(hsWindseaCombined)*tpWindseaCombined

# tC = np.array(tWave)


avgHs = np.nanmean(hsCombined)
hs98 = np.nanpercentile(hsCombined,98)
hs95 = np.nanpercentile(hsCombined,95)
hs90 = np.nanpercentile(hsCombined,90)
hs85 = np.nanpercentile(hsCombined,85)


years = np.arange(1979,2019)
winterHs = np.empty(len(years))
summerHs = np.empty(len(years))
summerTime = []
winterTime = []
winterHs90 = np.empty(len(years))
summerHs90 = np.empty(len(years))
winterHs95 = np.empty(len(years))
summerHs95 = np.empty(len(years))
winterHs80 = np.empty(len(years))
summerHs80 = np.empty(len(years))
for x in range(len(years)):
    t1 = datetime(years[x],10,1)
    t2 = datetime(years[x]+1,4,30)
    tempDates = np.where((tC > t1) & (tC < t2))
    winterHs[x] = np.nanmax(hsCombined[tempDates])
    winterHs90[x] = np.nanpercentile(hsCombined[tempDates],90)
    winterHs95[x] = np.nanpercentile(hsCombined[tempDates],95)
    winterHs80[x] = np.nanpercentile(hsCombined[tempDates],80)

    winterTime.append(datetime(years[x]+1,2,1))

    t3 = datetime(years[x]+1,5,1)
    t4 = datetime(years[x]+2,9,30)
    tempDates2 = np.where((tC > t3) & (tC < t4))
    summerHs[x] = np.nanmax(hsCombined[tempDates2])
    summerHs90[x] = np.nanpercentile(hsCombined[tempDates2],90)
    summerHs95[x] = np.nanpercentile(hsCombined[tempDates2],95)
    summerHs80[x] = np.nanpercentile(hsCombined[tempDates2],80)

    summerTime.append(datetime(years[x]+1,7,1))

plt.figure()
ax1 = plt.subplot2grid((2,2),(0,0),rowspan=1,colspan=2)
ax1.plot(winterTime,winterHs,'o-',label='max')
ax1.plot(winterTime,winterHs95,'o-',label='95')
ax1.plot(winterTime,winterHs90,'o-',label='90')
ax1.plot(winterTime,winterHs80,'o-',label='80')
ax1.legend()
ax1.set_title('Winter MHs')
ax2 = plt.subplot2grid((2,2),(1,0),rowspan=1,colspan=2)
ax2.plot(summerTime,summerHs,'o-',label='max')
ax2.plot(summerTime,summerHs95,'o-',label='95')
ax2.plot(summerTime,summerHs90,'o-',label='90')
ax2.plot(summerTime,summerHs80,'o-',label='80')
ax2.set_title('Summer Hs')
ax2.legend()






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


#DWT = ReadMatfile('/media/dylananderson/Elements1/NC_climate/Nags_Head_DWTs_49_w20minDates_2degree_plus5TCs3dayafterentering2.mat')
#DWT = ReadMatfile('/media/dylananderson/Elements1/NC_climate/Nags_Head_DWTS_25_w50minDates_plus5TCs_goodorder.mat')
DWT = ReadMatfile('/media/dylananderson/Elements/NC_climate/Nags_Head_DWTs_25_w50minDates_2degree_plus5TCs6daysAroundEntering.mat')
PCA = ReadMatfile('/media/dylananderson/Elements/NC_climate/Nags_Head_SLPs_2degree_memory_July2020.mat')
numDWTs = 30
mycols = ReadMatfile('/media/dylananderson/Elements/shusin6_contents/codes/mycmap_col.mat')
mycmap = mycols['mycmap_col']
# Need to get the dates for the bmus into the correct format (datetime)
def datevec2datetime(d_vec):
    '''
    Returns datetime list from a datevec matrix
    d_vec = [[y1 m1 d1 H1 M1],[y2 ,2 d2 H2 M2],..]
    '''
    return [datetime(d[0], d[1], d[2], d[3], d[4]) for d in d_vec]

dwtDates = np.array(datevec2datetime(DWT['DWT']['dates']))
dwtBmus = DWT['DWT']['bmus']
removeTooSoon = np.where(dwtDates < tC[-1]) #tC[-45])



dwtTimes = dwtDates[removeTooSoon]
dwtBMUS = dwtBmus[removeTooSoon]









stormHsInd = np.where((hsSmooth > hs90))
stormHsList = [list(group) for group in mit.consecutive_groups(stormHsInd[0])]


hsStormList = []
hsStormMaxList = []
hsStormMinList = []
tpStormList = []
tpStormMaxList = []
tpStormMinList = []
dmStormList = []
ntrStormList = []
timeStormList = []
hourStormList = []
indStormList = []
indNTRStormList = []
bmuStormList = []
for xx in range(len(stormHsList)-2):

    i1 = stormHsList[xx][0]
    i2 = stormHsList[xx][-1]
    t1 = tC[i1]
    t2 = tC[i2]
    nexti1 = stormHsList[xx+1][0]
    diff = nexti1 - i2
    # if tC[i1] > datetime.datetime(2019,1,1):
    #     numToBeat =
    if diff < 24:
        i2 = stormHsList[xx+1][-1]
        t2 = tC[i2]
        nexti1 = stormHsList[xx+2][0]
        diff2 = nexti1-i2
        if diff2 < 24:
            i2 = stormHsList[xx + 2][-1]
            t2 = tC[i2]

    tempWave = np.where((tC < t2) & (tC > t1))
    if len(tempWave[0]) > 12:
        t1 = tC[i1-12]
        t2 = tC[i2+12]
        t3 = tC[i2+36]
        t4 = tC[i1-18]

        tempWave = np.where((tC < t2) & (tC > t1))
        tempWaterLevel = np.where((tWaterLevelFRF < t2) & (tWaterLevelFRF > t1))
        tempBMU = np.where((dwtTimes < t3) & (dwtTimes > t4))
        indices2 = tempWaterLevel[0]
        indices = np.arange(i1-12,i2+12)


        bmuStormList.append(dwtBMUS[tempBMU])
        hsStormList.append(hsCombined[tempWave])
        ntrStormList.append(residualWaterLevelFRF[tempWaterLevel])
        # hsMaxList = np.append(hsMaxList,np.nanmax(hsCombined[tempWave]))
        # hsMinList = np.append(hsMinList,np.nanmin(hsCombined[tempWave]))
        tpStormList.append(tpCombined[tempWave])
        # tpMaxList = np.append(tpMaxList,np.nanmax(tpCombined[tempWave]))
        # tpMinList = np.append(tpMinList,np.nanmin(tpCombined[tempWave]))
        # dmStormList.append(dmCombined[tempWave])
        dmStormList.append(waveNorm[tempWave])

        timeStormList.append(tC[tempWave])
        indStormList.append(indices)
        indNTRStormList.append(indices2)


hsMaxStorm = []
hsMaxTimeInStorm = []
bmuMaxTimeInStorm = []
durationHoursStorm = []
dmAvgStorm = []
ntrMaxStorm = []
tpMaxStorm = []
timeStorm = []
for x in range(len(hsStormList)):
    tempHs = hsStormList[x]
    tempTp = tpStormList[x]
    timeStorm.append(timeStormList[x][12])
    tempMaxInds = np.where((np.nanmax(tempHs) == tempHs))
    hsMaxStorm.append(np.nanmax(tempHs))
    dayOf = int(np.floor(tempMaxInds[0][-1]/24))
    hsMaxTimeInStorm.append(tempMaxInds[0][-1])
    tpMaxStorm.append(tempTp[tempMaxInds[0][-1]])
    bmuMaxTimeInStorm.append(bmuStormList[x][dayOf])
    durationHoursStorm.append(len(tempHs))
    tempDm = dmStormList[x]
    dmAvgStorm.append(np.nanmean(tempDm))
    tempNTR = ntrStormList[x]
    if len(tempNTR) > 0:
        ntrMaxStorm.append(np.nanmax(tempNTR))
    else:
        ntrMaxStorm.append(np.nan)


hsStormArray = np.array(hsMaxStorm)
tpStormArray = np.array(tpMaxStorm)
dmStormArray = np.array(dmAvgStorm)
ntrStormArray = np.array(ntrMaxStorm)
durationStormArray = np.array(durationHoursStorm)
timeStormArray = np.array(timeStorm)

filteredHs = hsStormArray[~np.isnan(ntrStormArray)]
filteredTp = tpStormArray[~np.isnan(ntrStormArray)]
filteredDm = dmStormArray[~np.isnan(ntrStormArray)]
filteredDur = durationStormArray[~np.isnan(ntrStormArray)]
filteredTime = timeStormArray[~np.isnan(ntrStormArray)]
filteredNTR = ntrStormArray[~np.isnan(ntrStormArray)]

badNTRs = np.where((filteredNTR > 1.5))
filteredNTR[badNTRs] = 0.25

badTp = np.where((filteredHs > 5.5) & (filteredTp < 10))
filteredTp[badTp] = 10

plt.figure()
ax1 = plt.subplot2grid((5,5),(0,0),rowspan=1,colspan=1)
ax1.hist(filteredNTR)
ax1.set_xlabel('NTR (m)')
ax1.set_ylabel('NTR (m)')

ax1.yaxis.set_ticks([])

ax2 = plt.subplot2grid((5,5),(0,1),rowspan=1,colspan=1)
ax2.plot(filteredHs,filteredNTR,'k.')
ax3 = plt.subplot2grid((5,5),(0,2),rowspan=1,colspan=1)
ax3.plot(filteredTp,filteredNTR,'k.')
ax4 = plt.subplot2grid((5,5),(0,3),rowspan=1,colspan=1)
ax4.plot(filteredDm,filteredNTR,'k.')
ax5 = plt.subplot2grid((5,5),(0,4),rowspan=1,colspan=1)
ax5.plot(filteredDur,filteredNTR,'k.')
ax5.set_xlim([24, 300])

ax6 = plt.subplot2grid((5,5),(1,1),rowspan=1,colspan=1)
ax6.hist(filteredHs)
ax6.yaxis.set_ticks([])
ax6.set_xlabel('Hs (m)')
ax6.set_ylabel('Hs (m)')

ax7 = plt.subplot2grid((5,5),(1,2),rowspan=1,colspan=1)
ax7.plot(filteredTp,filteredHs,'k.')
ax8 = plt.subplot2grid((5,5),(1,3),rowspan=1,colspan=1)
ax8.plot(filteredDm,filteredHs,'k.')
ax9 = plt.subplot2grid((5,5),(1,4),rowspan=1,colspan=1)
ax9.plot(filteredDur,filteredHs,'k.')
ax9.set_xlim([24, 300])

ax10 = plt.subplot2grid((5,5),(2,2),rowspan=1,colspan=1)
ax10.hist(filteredTp)
ax10.set_xlabel('Period (x)')
ax10.yaxis.set_ticks([])
ax10.set_ylabel('Period (x)')
ax11 = plt.subplot2grid((5,5),(2,3),rowspan=1,colspan=1)
ax11.plot(filteredDm,filteredTp,'k.')
ax12 = plt.subplot2grid((5,5),(2,4),rowspan=1,colspan=1)
ax12.plot(filteredDur,filteredTp,'k.')
ax12.set_xlim([24, 300])

ax13 = plt.subplot2grid((5,5),(3,3),rowspan=1,colspan=1)
ax13.hist(filteredDm)
ax13.set_xlabel('Direction (avg)')
ax13.yaxis.set_ticks([])
ax13.set_ylabel('Direction (avg)')
ax14 = plt.subplot2grid((5,5),(3,4),rowspan=1,colspan=1)
ax14.plot(filteredDur,filteredDm,'k.')
ax14.set_xlim([24, 300])

ax15 = plt.subplot2grid((5,5),(4,4),rowspan=1,colspan=1)
ax15.hist(filteredDur)
ax15.set_xlim([24, 300])
ax15.yaxis.set_ticks([])
ax15.set_xlabel('Duration (hours)')



dataCop = []
for xx in range(len(filteredHs)):
    dataCop.append(list([filteredHs[xx],filteredTp[xx],filteredDm[xx],filteredNTR[xx],filteredDur[xx]]))


clusterPickle = 'filteredStorms.pickle'
output = {}
output['dataCop'] = dataCop
output['hsStormArray'] = hsStormArray
output['tpStormArray'] = tpStormArray
output['dmStormArray'] = dmStormArray
output['ntrStormArray'] = ntrStormArray
output['durationStormArray'] = durationStormArray
output['timeStormList'] = timeStormList
output['filteredHs'] = filteredHs
output['filteredTp'] = filteredTp
output['filteredDm'] = filteredDm
output['filteredNTR'] = filteredNTR
output['filteredDur'] = filteredDur
output['filteredTime'] = filteredTime

import pickle
with open(clusterPickle,'wb') as f:
    pickle.dump(output, f)


#dataCop = [list(filteredHs),list(filteredTp),list(filteredDm),list(filteredNTR),list(filteredDur)]
from copula import pyCopula
cop = pyCopula.Copula(dataCop)
samples = cop.gendata(10000)

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

data = np.empty((len(sampleHs),5))
data[:,0] = sampleHs
data[:,1] = sampleTp
data[:,2] = sampleDm
data[:,3] = sampleNTR
data[:,4] = sampleDur





ax21 = plt.subplot2grid((5,5),(1,0),rowspan=1,colspan=1)
ax21.plot(sampleNTR,sampleHs,'.',color=[0.5,0.5,0.5])
ax22 = plt.subplot2grid((5,5),(2,0),rowspan=1,colspan=1)
ax22.plot(sampleNTR,sampleTp,'.',color=[0.5,0.5,0.5])
ax23 = plt.subplot2grid((5,5),(3,0),rowspan=1,colspan=1)
ax23.plot(sampleNTR,sampleDm,'.',color=[0.5,0.5,0.5])
ax24 = plt.subplot2grid((5,5),(4,0),rowspan=1,colspan=1)
ax24.plot(sampleNTR,sampleDur,'.',color=[0.5,0.5,0.5])
ax24.set_xlabel('NTR (m)')
ax25 = plt.subplot2grid((5,5),(2,1),rowspan=1,colspan=1)
ax25.plot(sampleHs,sampleTp,'.',color=[0.5,0.5,0.5])
ax26 = plt.subplot2grid((5,5),(3,1),rowspan=1,colspan=1)
ax26.plot(sampleHs,sampleDm,'.',color=[0.5,0.5,0.5])
ax27 = plt.subplot2grid((5,5),(4,1),rowspan=1,colspan=1)
ax27.plot(sampleHs,sampleDur,'.',color=[0.5,0.5,0.5])
ax27.set_xlabel('Hs (m)')
ax28 = plt.subplot2grid((5,5),(3,2),rowspan=1,colspan=1)
ax28.plot(sampleTp,sampleDm,'.',color=[0.5,0.5,0.5])
ax29 = plt.subplot2grid((5,5),(4,2),rowspan=1,colspan=1)
ax29.plot(sampleTp,sampleDur,'.',color=[0.5,0.5,0.5])
ax29.set_xlabel('Tp (s)')
ax30 = plt.subplot2grid((5,5),(4,3),rowspan=1,colspan=1)
ax30.plot(sampleDm,sampleDur,'.',color=[0.5,0.5,0.5])
ax30.set_xlabel('Dm (s)')



mdaSub = mda.MaxDiss_Simplified_NoThreshold(data,200,[0,1,3,4],[2])

ax21.plot(mdaSub[0:100,3],mdaSub[0:100,0],'.',color='red')
ax22.plot(mdaSub[0:100,3],mdaSub[0:100,1],'.',color='red')
ax23.plot(mdaSub[0:100,3],mdaSub[0:100,2],'.',color='red')
ax24.plot(mdaSub[0:100,3],mdaSub[0:100,4],'.',color='red')
ax25.plot(mdaSub[0:100,0],mdaSub[0:100,1],'.',color='red')
ax26.plot(mdaSub[0:100,0],mdaSub[0:100,2],'.',color='red')
ax27.plot(mdaSub[0:100,0],mdaSub[0:100,4],'.',color='red')
ax28.plot(mdaSub[0:100,1],mdaSub[0:100,2],'.',color='red')
ax29.plot(mdaSub[0:100,1],mdaSub[0:100,4],'.',color='red')
ax30.plot(mdaSub[0:100,2],mdaSub[0:100,4],'.',color='red')

ax21.plot(mdaSub[100:200,3],mdaSub[100:200,0],'.',color='orange')
ax22.plot(mdaSub[100:200,3],mdaSub[100:200,1],'.',color='orange')
ax23.plot(mdaSub[100:200,3],mdaSub[100:200,2],'.',color='orange')
ax24.plot(mdaSub[100:200,3],mdaSub[100:200,4],'.',color='orange')
ax25.plot(mdaSub[100:200,0],mdaSub[100:200,1],'.',color='orange')
ax26.plot(mdaSub[100:200,0],mdaSub[100:200,2],'.',color='orange')
ax27.plot(mdaSub[100:200,0],mdaSub[100:200,4],'.',color='orange')
ax28.plot(mdaSub[100:200,1],mdaSub[100:200,2],'.',color='orange')
ax29.plot(mdaSub[100:2000,1],mdaSub[100:200,4],'.',color='orange')
ax30.plot(mdaSub[100:200,2],mdaSub[100:200,4],'.',color='orange')


clusterPickle = 'hypotheticalStorms.pickle'
output = {}
output['hypotheticalStorms'] = mdaSub

import pickle
with open(clusterPickle,'wb') as f:
    pickle.dump(output, f)
# ax15.hist(sampleDur)



order = DWT['DWT']['order']


import matplotlib.cm as cm

etcolors = cm.rainbow(np.linspace(0, 1, numDWTs-5))
tccolors = np.flipud(cm.gray(np.linspace(0,1,6)))

dwtcolors = np.vstack((etcolors,tccolors[1:,:]))

plt.figure()
for xx in range(len(order)):
    num = order[xx]
    tempInd = np.where((bmuMaxTimeInStorm==num))
    for ff in range(len(tempInd[0])):
        plt.scatter(durationHoursStorm[tempInd[0][ff]],hsMaxStorm[tempInd[0][ff]],c=dwtcolors[num-1])
plt.xlabel('Duration (hours)')
plt.ylabel('Max Wave Height')


plt.figure()
for xx in range(len(order)):
    num = order[xx]
    tempInd = np.where((bmuMaxTimeInStorm==num))
    for ff in range(len(tempInd[0])):
        plt.scatter(hsMaxStorm[tempInd[0][ff]],ntrMaxStorm[tempInd[0][ff]],c=dwtcolors[num-1])
plt.xlabel('Max Hs (m)')
plt.ylabel('Max NTR (m)')


# t1 = 275000
# t2 = 281000
t1 = 320000
t2 = 340000
plt.figure()
ax1 = plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=4)
ax1.plot(tC[t1:t2],hsCombined[t1:t2],'k')
# ax1.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])
ax1.set_xlim([datetime(2017,6,1), datetime(2017,12,1)])
indices = np.arange(715,750)
for xx in indices:
    tempIndex = indStormList[xx]
    ax1.plot(tC[tempIndex], hsCombined[tempIndex], 'r.')

ax2 = plt.subplot2grid((4,4),(1,0),rowspan=1,colspan=4)
ax2.plot(tC[t1:t2],tpCombined[t1:t2],'k')
# ax2.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])
ax2.set_xlim([datetime(2017,6,1), datetime(2017,12,1)])
for xx in indices:
    tempIndex = indStormList[xx]
    ax2.plot(tC[tempIndex], tpCombined[tempIndex], 'r.')

ax3 = plt.subplot2grid((4,4),(2,0),rowspan=1,colspan=4)
ax3.plot(tC[t1:t2],dmCombined[t1:t2],'k')
# ax3.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])
ax3.set_xlim([datetime(2017,6,1), datetime(2017,12,1)])
for xx in indices:
    tempIndex = indStormList[xx]
    ax3.plot(tC[tempIndex], dmCombined[tempIndex], 'r.')

t1 = 3130000
t2 = 3400000
ax4 = plt.subplot2grid((4,4),(3,0),rowspan=1,colspan=4)
ax4.plot(tWaterLevelFRF[t1:t2],residualWaterLevelFRF[t1:t2],'k')
# ax3.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])
ax4.set_xlim([datetime(2017,6,1), datetime(2017,12,1)])
for xx in indices:
    tempIndex = indNTRStormList[xx]
    ax4.plot(tWaterLevelFRF[tempIndex], residualWaterLevelFRF[tempIndex], 'r.')




t1 = 275000
t2 = 281000
t1 = 340000
t2 = 355000
plt.figure()
ax1 = plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=4)
ax1.plot(tC[t1:t2],hsCombined[t1:t2],'k')
# ax1.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])
ax4.set_xlim([datetime(2019,6,1), datetime(2019,12,1)])
# ax1.set_xlim([datetime(2017,6,1), datetime(2017,12,1)])
indices = np.arange(786,800)
for xx in indices:
    tempIndex = indStormList[xx]
    ax1.plot(tC[tempIndex], hsCombined[tempIndex], 'r.')

ax2 = plt.subplot2grid((4,4),(1,0),rowspan=1,colspan=4)
ax2.plot(tC[t1:t2],tpCombined[t1:t2],'k')
# ax2.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])
ax4.set_xlim([datetime(2019,6,1), datetime(2019,12,1)])
# ax2.set_xlim([datetime(2017,6,1), datetime(2017,12,1)])
for xx in indices:
    tempIndex = indStormList[xx]
    ax2.plot(tC[tempIndex], tpCombined[tempIndex], 'r.')

ax3 = plt.subplot2grid((4,4),(2,0),rowspan=1,colspan=4)
ax3.plot(tC[t1:t2],dmCombined[t1:t2],'k')
# ax3.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])
ax4.set_xlim([datetime(2019,6,1), datetime(2019,12,1)])
# ax3.set_xlim([datetime(2017,6,1), datetime(2017,12,1)])
for xx in indices:
    tempIndex = indStormList[xx]
    ax3.plot(tC[tempIndex], dmCombined[tempIndex], 'r.')

t1 = 2600000
t2 = 2750000
ax4 = plt.subplot2grid((4,4),(3,0),rowspan=1,colspan=4)
ax4.plot(tWaterLevelFRF[t1:t2],residualWaterLevelFRF[t1:t2],'k')
ax4.set_xlim([datetime(2019,6,1), datetime(2019,12,1)])
# ax4.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])
# ax4.set_xlim([datetime(2017,6,1), datetime(2017,12,1)])
for xx in indices:
    tempIndex = indNTRStormList[xx]
    ax4.plot(tWaterLevelFRF[tempIndex], residualWaterLevelFRF[tempIndex], 'r.')



asdfg

storm1 = dict()
storm1['waveTime'] = tC[indStormList[615]]
storm1['waveHs'] = hsCombined[indStormList[615]]
storm1['waveTp'] = tpCombined[indStormList[615]]
storm1['waveDm'] = dmCombined[indStormList[615]]
storm1['waterLevelTime'] = tWaterLevelFRF[indNTRStormList[615]]
storm1['waterLevelResidual'] = residualWaterLevelFRF[indNTRStormList[615]]
storm1['waterLevelTide'] = predictedWaterLevelFRF[indNTRStormList[615]]

storm2 = dict()
storm2['waveTime'] = tC[indStormList[616]]
storm2['waveHs'] = hsCombined[indStormList[616]]
storm2['waveTp'] = tpCombined[indStormList[616]]
storm2['waveDm'] = dmCombined[indStormList[616]]
storm2['waterLevelTime'] = tWaterLevelFRF[indNTRStormList[616]]
storm2['waterLevelResidual'] = residualWaterLevelFRF[indNTRStormList[616]]
storm2['waterLevelTide'] = predictedWaterLevelFRF[indNTRStormList[616]]

storm3 = dict()
storm3['waveTime'] = tC[indStormList[617]]
storm3['waveHs'] = hsCombined[indStormList[617]]
storm3['waveTp'] = tpCombined[indStormList[617]]
storm3['waveDm'] = dmCombined[indStormList[617]]
storm3['waterLevelTime'] = tWaterLevelFRF[indNTRStormList[617]]
storm3['waterLevelResidual'] = residualWaterLevelFRF[indNTRStormList[617]]
storm3['waterLevelTide'] = predictedWaterLevelFRF[indNTRStormList[617]]

storm4 = dict()
storm4['waveTime'] = tC[indStormList[618]]
storm4['waveHs'] = hsCombined[indStormList[618]]
storm4['waveTp'] = tpCombined[indStormList[618]]
storm4['waveDm'] = dmCombined[indStormList[618]]
storm4['waterLevelTime'] = tWaterLevelFRF[indNTRStormList[618]]
storm4['waterLevelResidual'] = residualWaterLevelFRF[indNTRStormList[618]]
storm4['waterLevelTide'] = predictedWaterLevelFRF[indNTRStormList[618]]

storm5 = dict()
storm5['waveTime'] = tC[indStormList[619]]
storm5['waveHs'] = hsCombined[indStormList[619]]
storm5['waveTp'] = tpCombined[indStormList[619]]
storm5['waveDm'] = dmCombined[indStormList[619]]
storm5['waterLevelTime'] = tWaterLevelFRF[indNTRStormList[619]]
storm5['waterLevelResidual'] = residualWaterLevelFRF[indNTRStormList[619]]
storm5['waterLevelTide'] = predictedWaterLevelFRF[indNTRStormList[619]]



clusterPickle = 'juneToNovemberStorms.pickle'
output = {}
output['storm1'] = storm1
output['storm2'] = storm2
output['storm3'] = storm3
output['storm4'] = storm4
output['storm5'] = storm5

import pickle
with open(clusterPickle,'wb') as f:
    pickle.dump(output, f)


# maxHsStorm = np.empty(len(hsStormList),)
# durationStorm = np.empty(len(hsStormList),)
#
# for xx in range(len(hsStormList)):
#     maxHsStorm[xx] = np.nanmax(hsStormList[xx])
#     durationStorm[xx] = len(hsStormList[xx])

















# Next, we need to split the waves up into each DWT
numDWTs = 30
dwtHs = []
dwtMaxHs = []
# dwtHshydro = []
# dwtHsSwell = []
# dwtHsWindsea = []
dwtTp = []
# dwtTpSwell = []
# dwtTpWindsea = []
dwtDm = []
# dwtDmSwell = []
# dwtDmWindsea = []
dwtLWP = []
dwtWE = []
dwtT = []
for xx in range(numDWTs):
    tempHs = []
    tempMaxHs = []
    tempHydro = []
    # tempHsSwell = []
    # tempHsWindsea = []
    # tempTpSwell = []
    # tempTpWindsea = []
    # tempDmSwell = []
    # tempDmWindsea = []
    tempTp = []
    tempDm = []
    tempLWP = []
    tempWE = []
    tempTi = []

    wInd = np.where((dwtBMUS[0:-1] == (xx+1)))
    for tt in range(len(wInd[0])):
        tempT = np.where((tC < dwtTimes[wInd[0][tt]+1]) & (tC > dwtTimes[wInd[0][tt]]))
        if len(tempT[0]) > 0:
            tempMaxHs = np.append(tempMaxHs, np.nanmax(hsCombined[tempT]))

        tempHs = np.append(tempHs, hsCombined[tempT])
        # tempHydro.append(hsCombined[tempT])
        # tempHsSwell = np.append(tempHsSwell, hsSwellCombined[tempT])
        # tempHsWindsea = np.append(tempHsWindsea, hsWindseaCombined[tempT])
        # tempTpSwell = np.append(tempTpSwell, tpSwellCombined[tempT])
        # tempTpWindsea = np.append(tempTpWindsea, tpWindseaCombined[tempT])
        # tempDmSwell = np.append(tempDmSwell, waveNormSwell[tempT])
        # tempDmWindsea = np.append(tempDmWindsea, waveNormWindsea[tempT])
        tempTp = np.append(tempTp, tpCombined[tempT])
        tempDm = np.append(tempDm, waveNorm[tempT])
        tempLWP = np.append(tempLWP, lwpC[tempT])
        tempWE = np.append(tempWE, weC[tempT])
        tempTi = np.append(tempTi, tC[tempT])

    # dwtHshydro.append(tempHydro)
    dwtHs.append(tempHs)
    dwtMaxHs.append(tempMaxHs)
    # dwtHsWindsea.append(tempHsWindsea)
    # dwtHsSwell.append(tempHsSwell)
    # dwtTpWindsea.append(tempTpWindsea)
    # dwtTpSwell.append(tempTpSwell)
    # dwtDmWindsea.append(tempDmWindsea)
    # dwtDmSwell.append(tempDmSwell)
    dwtTp.append(tempTp)
    dwtDm.append(tempDm)
    dwtLWP.append(tempLWP)
    dwtWE.append(tempWE)
    dwtT.append(tempTi)

# Next, we need to split the waves up into each DWT
# numDWTs = 30
dwtMaxHsHydro = []
dwtMinHsHydro = []
dwtHsHydro = []
dwtTpHydro = []
dwtMaxTpHydro = []
dwtMinTpHydro = []
dwtDmHydro = []
dwtMeanDmHydro = []
dwtTimeHydro = []

for xx in range(numDWTs):
    tempMaxHsHydro = []
    tempMinHsHydro = []
    tempHsHydro = []
    tempTpHydro = []
    tempMaxTpHydro = []
    tempMinTpHydro = []
    tempDmHydro = []
    tempMeanDmHydro = []
    tempTimeHydro = []

    wInd = np.where((dwtBMUS[0:-1] == (xx+1)))
    wHydroList = [list(group) for group in mit.consecutive_groups(wInd[0])]

    for tt in range(len(wHydroList)):
        tempT = np.where((tC < dwtTimes[wHydroList[tt][-1]+1]) & (tC > dwtTimes[wHydroList[tt][0]]))
        if len(tempT[0]) > 0:
            tempMaxHsHydro.append(np.nanmax(hsCombined[tempT]))
            tempMinHsHydro.append(np.nanmin(hsCombined[tempT]))
            tempMaxTpHydro.append(np.nanmax(tpCombined[tempT]))
            tempMinTpHydro.append(np.nanmin(tpCombined[tempT]))
            tempMeanDmHydro.append(np.nanmean(dmCombined[tempT]))

        tempHsHydro.append(hsCombined[tempT])
        tempTpHydro.append(tpCombined[tempT])
        tempDmHydro.append(dmCombined[tempT])
        tempTimeHydro.append(tC[tempT])

    dwtHsHydro.append(tempHsHydro)
    dwtMaxHsHydro.append(tempMaxHsHydro)
    dwtMinHsHydro.append(tempMinHsHydro)
    dwtTpHydro.append(tempTpHydro)
    dwtMaxTpHydro.append(tempMaxTpHydro)
    dwtMinTpHydro.append(tempMinTpHydro)
    dwtDmHydro.append(tempDmHydro)
    dwtTimeHydro.append(tempTimeHydro)




hydroList = [list(group) for group in mit.consecutive_groups(dwtBMUS)]

from itertools import groupby
grouped_L = [(k, sum(1 for i in g)) for k,g in groupby(dwtBMUS)]

cursor = 0
hydrographs = []
for k, l in grouped_L:
    hydrographs.append([np.arange(cursor, cursor + l)])
    cursor += l
    # if l < 5:
    #     hydrographs.append([np.arange(cursor, cursor + l)])
    #     cursor += l
    # else:
    #     j = l/5
    #     rem = np.remainder(l,5)
    #     fl = np.floor(j)-1
    #     for f in range(int(fl)):
    #         hydrographs.append([np.arange(cursor, cursor + 5)])
    #         cursor += 5
    #     if rem > 0:
    #         hydrographs.append([np.arange(cursor, cursor + rem)])
    #         cursor += rem

hsList = []
hsMaxList = []
hsMinList = []
tpList = []
tpMaxList = []
tpMinList = []
dmList = []
timeList = []
dwtList = []
dayList = []
bmuList = []
for xx in range(len(hydrographs)-1):
    i1 = hydrographs[xx][0][0]
    i2 = hydrographs[xx][0][-1]
    t1 = dwtTimes[i1]
    t2 = dwtTimes[i2+1]
    tempWave = np.where((tC < t2) & (tC > t1))
    # hsCombinedTemp = hsCombined[tempWave].astype(int)
    # array_of_tuples = map(tuple, hsCombinedTemp)
    # tuple_of_tuples = tuple(array_of_tuples)
    #
    # hsList.append(tuple_of_tuples)
    hsList.append(hsCombined[tempWave])
    # hsMaxList = np.append(hsMaxList,np.nanmax(hsCombined[tempWave]))
    # hsMinList = np.append(hsMinList,np.nanmin(hsCombined[tempWave]))
    tpList.append(tpCombined[tempWave])
    # tpMaxList = np.append(tpMaxList,np.nanmax(tpCombined[tempWave]))
    # tpMinList = np.append(tpMinList,np.nanmin(tpCombined[tempWave]))
    dmList.append(dmCombined[tempWave])
    timeList.append(tC[tempWave])
    dwtList.append(dwtBMUS[i1])
    bmuList = np.append(bmuList, dwtBMUS[i1])
    dayList = np.append(dayList, len(hydrographs[xx][0]))


from scipy.interpolate import interp1d

x = np.linspace(0, 1, num=50)

plt.figure()
gs2 = gridspec.GridSpec(6, 5)

for xx in range(numDWTs):
    ax = plt.subplot(gs2[xx])
    getInd = np.where((bmuList == (xx+1)))

    c = 0
    for k in getInd[0]:
        og = hsList[k]

        if len(og)>0:

            ox = np.linspace(0,1,len(og))
            nx = np.linspace(0,1,40)
            f = interp1d(ox,og,'cubic')
            nf = f(nx)
            ax.plot(nx,nf)
            if c == 0:
                newHydro = nf
            else:
                newHydro = np.vstack((newHydro,nf))
            c += 1
    ax.plot(nx,np.nanmean(newHydro,axis=0),'k',linewidth=3)




# # Next, we need to split the water levels up into each DWT
# dwtWL = []
# dwtPredWL = []
# dwtResWL = []
#
# for xx in range(numDWTs):
#     tempWL = []
#     tempPredWL = []
#     tempResWL = []
#
#     wInd = np.where((dwtBMUS[0:-1] == (xx+1)))
#     for tt in range(len(wInd[0])):
#         tempT = np.where((tWaterLevelFRF < dwtTimes[wInd[0][tt]+1]) & (tWaterLevelFRF > dwtTimes[wInd[0][tt]]))
#         tempWL = np.append(tempWL, waterLevelFRF[tempT])
#         tempPredWL = np.append(tempPredWL, predictedWaterLevelFRF[tempT])
#         tempResWL= np.append(tempResWL, residualWaterLevelFRF[tempT])
#
#     dwtWL.append(tempWL)
#     dwtPredWL.append(tempPredWL)
#     dwtResWL.append(tempResWL)
#

meanDWTHs = np.zeros((np.shape(order)))
for xx in range(numDWTs):
    data = dwtHs[xx]
    meanDWTHs[xx] = np.nanmean(data)

newOrder = np.argsort(meanDWTHs)

orderedDWTs = np.zeros((np.shape(dwtBMUS)))
for xx in range(numDWTs):
    dateInd = np.where((dwtBMUS == xx))
    orderedDWTs[dateInd] = newOrder[xx]



#plt.style.use('dark_background')



dist_space = np.linspace(0, 4, 80)
fig = plt.figure(figsize=(10,10))

gs2 = gridspec.GridSpec(6, 5)

colorparam = np.zeros((numDWTs,))
counter = 0
plotIndx = 0
plotIndy = 0
for xx in range(numDWTs):
    #dwtInd = xx
    dwtInd = order[xx]-1
    #dwtInd = newOrder[xx]

    #ax = plt.subplot2grid((6, 5), (plotIndx, plotIndy), rowspan=1, colspan=1)
    ax = plt.subplot(gs2[xx])

    # normalize = mcolors.Normalize(vmin=np.min(colorparams), vmax=np.max(colorparams))
    normalize = mcolors.Normalize(vmin=.5, vmax=1.8)

    ax.set_xlim([0, 3])
    ax.set_ylim([0, 2])
    data = dwtHs[dwtInd]
    if len(data) > 0:
        kde = gaussian_kde(data)
        colorparam[counter] = np.nanmean(data)
        colormap = cm.Reds
        color = colormap(normalize(colorparam[counter]))
        ax.plot(dist_space, kde(dist_space), linewidth=1, color=color)
        ax.spines['bottom'].set_color([0.5, 0.5, 0.5])
        ax.spines['top'].set_color([0.5, 0.5, 0.5])
        ax.spines['right'].set_color([0.5, 0.5, 0.5])
        ax.spines['left'].set_color([0.5, 0.5, 0.5])
        # ax.text(1.8, 1, np.round(colorparam*100)/100, fontweight='bold')

    else:
        ax.spines['bottom'].set_color([0.3, 0.3, 0.3])
        ax.spines['top'].set_color([0.3, 0.3, 0.3])
        ax.spines['right'].set_color([0.3, 0.3, 0.3])
        ax.spines['left'].set_color([0.3, 0.3, 0.3])

    if plotIndx < 6:
        ax.xaxis.set_ticks([])
        ax.xaxis.set_ticklabels([])
    if plotIndy > 0:
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])

    counter = counter + 1
    if plotIndy < 7:
        plotIndy = plotIndy + 1
    else:
        plotIndy = 0
        plotIndx = plotIndx + 1

    ax.yaxis.set_ticklabels([])
    ax.yaxis.set_ticks([])

    print(plotIndy, plotIndx)

plt.show()
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map.set_array(colorparam)
fig.subplots_adjust(right=0.86)
cbar_ax = fig.add_axes([0.89, 0.15, 0.02, 0.7])
cbar = fig.colorbar(s_map, cax=cbar_ax)
cbar.set_label('Mean Hs (m)')






repmatDesviacion = np.tile(DWT['PCA']['Desviacion'], (25,1))
repmatMedia = np.tile(DWT['PCA']['Media'], (25,1))
Km_ = np.multiply(DWT['KMA']['centroids'],repmatDesviacion) + repmatMedia

[mK, nK] = np.shape(Km_)

Km_slp = Km_[:,0:int(nK/2)]
Km_grd = Km_[:,int(nK/2):]
X_B = DWT['X_B']
Y_B = DWT['Y_B']
SLP_C = DWT['SLP_C']
Km_slp = Km_slp[:,0:len(X_B)]
Km_grd = Km_grd[:,1:len(X_B)]

Xs = np.arange(np.min(X_B),np.max(X_B),2)
#Xs = np.arange(-(360-np.min(X_B)),(np.max(X_B)-360))
Ys = np.arange(np.min(Y_B),np.max(Y_B),2)
lenXB = len(X_B)
[XR,YR] = np.meshgrid(Xs,Ys)
#sea_nodes = np.arange(0,len(X_B))
#sea_nodes = np.zeros((np.shape(X_B)))
sea_nodes = []
for qq in range(lenXB-1):
    sea_nodes.append(np.where((XR == X_B[qq]) & (YR == Y_B[qq])))
    #print(indexTest[0])

flat_list = [item for sublist in sea_nodes for item in sublist]


plt.style.use('default')
from mpl_toolkits.basemap import Basemap
fig2 = plt.figure(figsize=(8,10))
gs1 = gridspec.GridSpec(6, 5)
gs1.update(wspace=0.00, hspace=0.00) # set the spacing between axes.
counter = 0
plotIndx = 0
plotIndy = 0
for qq in range(numDWTs):
#qq = 0

    #ax = plt.subplot2grid((7, 7), (plotIndx, plotIndy), rowspan=1, colspan=1)
    ax = plt.subplot(gs1[qq])
    #num = qq
    num = order[qq]
    #num = newOrder[qq]+1
    wt = np.ones((np.shape(XR))) * np.nan

    if num < 25:
        num_index = np.where((dwtBMUS == num))
        temp = Km_slp[(num-1),:]/100 - np.nanmean(SLP_C,axis=0)/100
    else:
        num_index = np.where((dwtBmus == num))
        temp = np.nanmean(SLP_C[num_index,:],axis=1)/100 - np.nanmean(SLP_C,axis=0)/100
        temp= temp.flatten()
    for tt in range(len(sea_nodes)):
        wt[sea_nodes[tt]] = temp[tt]
#reshaped = np.tile(wt,(np.shape(XR)))


#m = Basemap(projection='cyl',lon_0=320,lat_0=0)
    m = Basemap(projection='merc',llcrnrlat=-40,urcrnrlat=50,llcrnrlon=275,urcrnrlon=370,lat_ts=10,resolution='c')
#m = Basemap(projection='stere',lon_0=-120,lat_0=20,llcrnrlat=-10,llcrnrlon=50,urcrnrlat=50, urcrnrlon=190)
#m = Basemap(width=12000000,height=8000000,resolution='l',projection='stere',lat_ts=30,lat_0=20,lon_0=-70.)
    m.fillcontinents(color=dwtcolors[qq])
    m.drawcoastlines()   #畫海岸線
#parallels = np.arange(-90.,90,30.)
#meridians = np.arange(0.,360.,20.)
#m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)  # 緯度度線、在左右兩邊加緯度標籤
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    clevels = np.arange(-17,17,1)
    cx,cy =m(XR,YR)  # convert to map projection coordinate
#CS = m.pcolormesh(cx,cy,wt,clevels,cmap=cm.jet,shading='gouraud')
    CS = m.contourf(cx,cy,wt,clevels,vmin=-7.5,vmax=7.5,cmap=cm.RdBu_r,shading='gouraud')

    #plt.colorbar(CS,orientation='horizontal')
    if plotIndx < 7:
        ax.xaxis.set_ticks([])
        ax.xaxis.set_ticklabels([])
    if plotIndy > 0:
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
    counter = counter + 1
    if plotIndy < 7:
        plotIndy = plotIndy + 1
    else:
        plotIndy = 0
        plotIndx = plotIndx + 1
#date=time_d[0].strftime('%Y/%m/%d')
#plt.title('Cylindrical, T at 1000 hPa, GFS'+date)
#plt.tight_layout()
colormap = cm.RdBu_r
normalize = mcolors.Normalize(vmin=-7.5, vmax=7.5)

s_map2 = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map2.set_array(colorparam)
fig.subplots_adjust(right=0.85)
cbar_ax2 = fig2.add_axes([0.91, 0.15, 0.02, 0.7])
cbar2 = fig2.colorbar(s_map2, cax=cbar_ax2)
cbar2.set_label('slp anom (mbar)')

plt.show()




t1 = 275000
t2 = 281000
plt.figure()
ax1 = plt.subplot2grid((4,4),(2,0),rowspan=1,colspan=4)
ax1.plot(tC[t1:t2],hsCombined[t1:t2],'k')
ax1.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])

ax2 = plt.subplot2grid((4,4),(1,0),rowspan=1,colspan=4)
ax2.plot(tC[t1:t2],tpCombined[t1:t2],'k')
ax2.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])

ax3 = plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=4)
ax3.plot(tC[t1:t2],dmCombined[t1:t2],'k')
ax3.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])

ax4 = plt.subplot2grid((4,4),(3,0),rowspan=1,colspan=4)
# ax3.set_xlim([datetime(2011,6,1), datetime(2012,1,1)])

#indices = np.arange(2730,2800)
indices = np.arange(2540,2800)

for xx in indices:
    t1 = timeList[xx][0]
    t2 = timeList[xx][-1]
    # # convert to matplotlib date representation
    start = mdates.date2num(t1)
    end = mdates.date2num(t2)
    width = end - start
    # Plot rectangle
    qq = dwtList[xx]
    if qq > 30:
        rect = Rectangle((start, 0), width, 1, color='black')
    else:
        rect = Rectangle((start, 0), width, 1, color=dwtcolors[qq-1])
    ax4.add_patch(rect)

ax4.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])


# locator = mdates.AutoDateLocator(minticks=3)
# formatter = mdates.AutoDateFormatter(locator)
# ax4.xaxis.set_major_locator(locator)
# ax4.xaxis.set_major_formatter(formatter)












