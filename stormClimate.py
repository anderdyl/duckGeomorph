
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
#bmuMaxTimeInStorm = []
durationHoursStorm = []
dmAvgStorm = []
ntrMaxStorm = []
tpMaxStorm = []
timeStorm = []
timeStormEnd = []
for x in range(len(hsStormList)):
    tempHs = hsStormList[x]
    tempTp = tpStormList[x]
    timeStorm.append(timeStormList[x][12])
    timeStormEnd.append(timeStormList[x][-12])
    tempMaxInds = np.where((np.nanmax(tempHs) == tempHs))
    hsMaxStorm.append(np.nanmax(tempHs))
    dayOf = int(np.floor(tempMaxInds[0][-1]/24))
    hsMaxTimeInStorm.append(tempMaxInds[0][-1])
    tpMaxStorm.append(tempTp[tempMaxInds[0][-1]])
    #bmuMaxTimeInStorm.append(bmuStormList[x][dayOf])
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
timeStormEndArray = np.array(timeStormEnd)

filteredHs = hsStormArray[~np.isnan(ntrStormArray)]
filteredTp = tpStormArray[~np.isnan(ntrStormArray)]
filteredDm = dmStormArray[~np.isnan(ntrStormArray)]
filteredDur = durationStormArray[~np.isnan(ntrStormArray)]
filteredTime = timeStormArray[~np.isnan(ntrStormArray)]
filteredTimeEnd = timeStormEndArray[~np.isnan(ntrStormArray)]
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








asdfghjkl

dataCop = []
for xx in range(len(filteredHs)):
    dataCop.append(list([filteredHs[xx],filteredTp[xx],filteredDm[xx],filteredNTR[xx],filteredDur[xx]]))


clusterPickle = 'filteredStormsLatest.pickle'
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
output['filteredTimeEnd'] = filteredTimeEnd
import pickle
with open(clusterPickle,'wb') as f:
    pickle.dump(output, f)


def datetime2datevec(dtime):
    'Return matlab date vector from datetimes'
    return [dtime.year, dtime.month, dtime.day]

mdateVecTime = [datetime2datevec(x) for x in filteredTime]
mdateVecTimeEnd = [datetime2datevec(x) for x in filteredTimeEnd]
mdateVecTC= [datetime2datevec(x) for x in tC]


filteredOutput = dict()
filteredOutput['filteredHs'] = filteredHs
filteredOutput['filteredTp'] = filteredTp
filteredOutput['filteredDm'] = filteredDm
filteredOutput['filteredNTR'] = filteredNTR
filteredOutput['filteredTimeEnd'] = mdateVecTimeEnd
filteredOutput['filteredTime'] = mdateVecTime
filteredOutput['filteredDur'] = filteredDur
#filteredOutput['timeStormList'] = timeStormList
filteredOutput['hsCombined'] = hsCombined
filteredOutput['tpCombined'] = tpCombined
filteredOutput['dmCombined'] = dmCombined
filteredOutput['tCombined'] = mdateVecTC

import scipy.io
scipy.io.savemat('stormClimateDetails.mat',filteredOutput)


