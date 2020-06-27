
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


dbfile = open('sandbarsSouthernTransect.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()
alllines = data['alllines']
xinterp = data['xinterp']
time = data['time']


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



mat = scipy.io.loadmat('matlabsSouthernTransectClusters.mat')

numCenters = int(mat['clusters']['NumberCenters'][0][0].flatten())
group = mat['clusters']['Group'][0][0]
bmus = mat['clusters']['PatternsGroup'][0][0].flatten()
order = mat['order'].flatten()
#sorted_bmus = np.zeros((len(bmus),))
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




# def transitionDictionary(bmus,dates):

bins = dict()
date = dict()
nextbin = dict()
nextdate = dict()
prevbin = dict()
prevdate = dict()
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



m,m2 = transition_matrix(sorted_bmus)
for row in m: print(' '.join('{0:.2f}'.format(x) for x in row))

flat_list = [item for sublist in m for item in sublist]
flatarray = np.asarray(flat_list)
flatarray.resize(numClusters, numClusters)

flat_list2 = [item for sublist in m2 for item in sublist]
flatarray2 = np.asarray(flat_list2)
flatarray2.resize(numClusters, numClusters)


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


hsCombined = Hs
tpCombined = Tp
dmCombined = Dm
hsSwellCombined = hsSwell
tpSwellCombined = tpSwell
dmSwellCombined = dmSwell
hsWindseaCombined = hsWindsea
tpWindseaCombined = tpWindsea
dmWindseaCombined = dmWindsea

badDirs = np.where((dmCombined > 360))
dmCombined[badDirs] = dmCombined[badDirs]*np.nan

badDirsSwell = np.where((dmSwellCombined > 360))
dmSwellCombined[badDirsSwell] = dmSwellCombined[badDirsSwell]*np.nan
badDirsWindsea = np.where((dmWindseaCombined > 360))
dmWindseaCombined[badDirsWindsea] = dmWindseaCombined[badDirsWindsea]*np.nan

waveNorm = dmCombined - 72
neg = np.where((waveNorm > 180))
waveNorm[neg[0]] = waveNorm[neg[0]]-360
offpos = np.where((waveNorm>90))
offneg = np.where((waveNorm<-90))
waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
waveNorm[offneg[0]] = waveNorm[offneg[0]]*0

waveNormSwell = dmSwellCombined - 72
negSwell = np.where((waveNormSwell > 180))
waveNormSwell[negSwell[0]] = waveNormSwell[negSwell[0]]-360
offposSwell = np.where((waveNormSwell>90))
offnegSwell = np.where((waveNormSwell<-90))
waveNormSwell[offposSwell[0]] = waveNormSwell[offposSwell[0]]*0
waveNormSwell[offnegSwell[0]] = waveNormSwell[offnegSwell[0]]*0

waveNormWindsea = dmWindseaCombined - 72
negWindsea = np.where((waveNormWindsea > 180))
waveNormWindsea[negWindsea[0]] = waveNormWindsea[negWindsea[0]]-360
offposWindsea = np.where((waveNormWindsea>90))
offnegWindsea = np.where((waveNormWindsea<-90))
waveNormWindsea[offposWindsea[0]] = waveNormWindsea[offposWindsea[0]]*0
waveNormWindsea[offnegWindsea[0]] = waveNormWindsea[offnegWindsea[0]]*0

lwpC = 1025*np.square(hsCombined)*tpCombined*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
weC = np.square(hsCombined)*tpCombined

lwpSwell = 1025*np.square(hsSwellCombined)*tpSwellCombined*(9.81/(64*np.pi))*np.cos(waveNormSwell*(np.pi/180))*np.sin(waveNormSwell*(np.pi/180))
weSwell = np.square(hsSwellCombined)*tpSwellCombined

lwpWindsea = 1025*np.square(hsWindseaCombined)*tpWindseaCombined*(9.81/(64*np.pi))*np.cos(waveNormWindsea*(np.pi/180))*np.sin(waveNormWindsea*(np.pi/180))
weWindsea = np.square(hsWindseaCombined)*tpWindseaCombined
tWave = [DT.datetime.fromtimestamp(x) for x in timeWave]

tC = np.array(tWave)

from itertools import groupby
from operator import itemgetter
import more_itertools as mit

wHs = []
wTp = []
wDm = []
wLWP = []
wWE = []
wT = []
swellHs = []
swellTp = []
swellDm = []
swellLWP = []
windseaHs = []
windseaTp = []
windseaDm = []
windseaLWP = []
for xx in range(numClusters):
    innerListHs = []
    innerListTp = []
    innerListDm = []
    innerListLWP = []
    innerListWE = []
    innerListT = []
    innerListSwellHs = []
    innerListSwellTp = []
    innerListSwellDm = []
    innerListSwellLWP = []
    innerListWindseaHs = []
    innerListWindseaTp = []
    innerListWindseaDm = []
    innerListWindseaLWP = []
    for yy in range(numClusters):
        # if both are equal then the wave conditions didn't cause a transition and
        # everything between the last two dates should be added to a distribution

        # if yy is different...
        wInd = np.where(prevbin[xx]==yy)
        if len(wInd[0])>0:
            tempHs = []
            tempTp = []
            tempDm = []
            tempLWP = []
            tempWE = []
            tempTi = []
            tempSwellHs = []
            tempSwellTp = []
            tempSwellDm = []
            tempSwellLWP = []
            tempWindseaHs = []
            tempWindseaTp = []
            tempWindseaDm = []
            tempWindseaLWP = []
            for tt in range(len(wInd[0])):
                tempT = np.where((tC < date[xx][wInd[0][tt]]) & (tC > prevdate[xx][wInd[0][tt]]))



                # lets look for storms in the 26 m record
                stormHsInd = np.where((hsCombined[tempT]>2))
                stormHsList = [list(group) for group in mit.consecutive_groups(stormHsInd[0])]



                tempHs = np.append(tempHs,hsCombined[tempT])
                tempTp = np.append(tempTp,tpCombined[tempT])
                tempDm = np.append(tempDm,waveNorm[tempT])
                tempLWP = np.append(tempLWP,lwpC[tempT])
                tempWE = np.append(tempWE,weC[tempT])
                tempTi = np.append(tempTi,tC[tempT])
                tempSwellHs = np.append(tempSwellHs,hsSwellCombined[tempT])
                tempSwellTp = np.append(tempSwellTp,tpSwellCombined[tempT])
                tempSwellDm = np.append(tempSwellDm,waveNormSwell[tempT])
                tempSwellLWP = np.append(tempSwellLWP,lwpSwell[tempT])
                tempWindseaHs = np.append(tempWindseaHs,hsWindseaCombined[tempT])
                tempWindseaTp = np.append(tempWindseaTp,tpWindseaCombined[tempT])
                tempWindseaDm = np.append(tempWindseaDm,waveNormWindsea[tempT])
                tempWindseaLWP = np.append(tempWindseaLWP,lwpWindsea[tempT])
        else:
            tempHs = []
            tempTp = []
            tempDm = []
            tempLWP =[]
            tempWE = []
            tempTi = []
            tempSwellHs = []
            tempSwellTp = []
            tempSwellDm = []
            tempSwellLWP = []
            tempWindseaHs = []
            tempWindseaTp = []
            tempWindseaDm = []
            tempWindseaLWP = []
        innerListHs.append(tempHs)
        innerListTp.append(tempTp)
        innerListDm.append(tempDm)
        innerListLWP.append(tempLWP)
        innerListWE.append(tempWE)
        innerListT.append(tempTi)
        innerListSwellHs.append(tempSwellHs)
        innerListSwellTp.append(tempSwellTp)
        innerListSwellDm.append(tempSwellDm)
        innerListSwellLWP.append(tempSwellLWP)
        innerListWindseaHs.append(tempWindseaHs)
        innerListWindseaTp.append(tempWindseaTp)
        innerListWindseaDm.append(tempWindseaDm)
        innerListWindseaLWP.append(tempWindseaLWP)
    wHs.append(innerListHs)
    wTp.append(innerListTp)
    wDm.append(innerListDm)
    wLWP.append(innerListLWP)
    wWE.append(innerListWE)
    wT.append(innerListT)
    swellHs.append(innerListSwellHs)
    swellTp.append(innerListSwellTp)
    swellDm.append(innerListSwellDm)
    swellLWP.append(innerListSwellLWP)
    windseaHs.append(innerListWindseaHs)
    windseaTp.append(innerListWindseaTp)
    windseaDm.append(innerListWindseaDm)
    windseaLWP.append(innerListWindseaLWP)


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

stayHsSwell = []
stayTpSwell = []
stayDmSwell = []
stayLWPSwell = []
downHsSwell = []
downTpSwell = []
downDmSwell = []
downLWPSwell = []
upHsSwell = []
upTpSwell = []
upDmSwell = []
upLWPSwell = []

stayHsWindsea = []
stayTpWindsea = []
stayDmWindsea = []
stayLWPWindsea = []
downHsWindsea = []
downTpWindsea = []
downDmWindsea = []
downLWPWindsea = []
upHsWindsea = []
upTpWindsea = []
upDmWindsea = []
upLWPWindsea = []

for xx in range(numClusters):
    for yy in range(numClusters):
        if xx == yy:
            if len(wHs[xx][yy]) > 0:
                stayHs.append(wHs[xx][yy])
                stayTp.append(wTp[xx][yy])
                stayDm.append(wDm[xx][yy])
                stayLWP.append(wLWP[xx][yy])
                stayHsSwell.append(swellHs[xx][yy])
                stayTpSwell.append(swellTp[xx][yy])
                stayDmSwell.append(swellDm[xx][yy])
                stayLWPSwell.append(swellLWP[xx][yy])

                stayHsWindsea.append(windseaHs[xx][yy])
                stayTpWindsea.append(windseaTp[xx][yy])
                stayDmWindsea.append(windseaDm[xx][yy])
                stayLWPWindsea.append(windseaLWP[xx][yy])
        if xx < yy:
            if len(wHs[xx][yy]) > 0:
                downHs.append(wHs[xx][yy])
                downTp.append(wTp[xx][yy])
                downDm.append(wDm[xx][yy])
                downLWP.append(wLWP[xx][yy])
                downHsSwell.append(swellHs[xx][yy])
                downTpSwell.append(swellTp[xx][yy])
                downDmSwell.append(swellDm[xx][yy])
                downLWPSwell.append(swellLWP[xx][yy])

                downHsWindsea.append(windseaHs[xx][yy])
                downTpWindsea.append(windseaTp[xx][yy])
                downDmWindsea.append(windseaDm[xx][yy])
                downLWPWindsea.append(windseaLWP[xx][yy])
        if yy < xx:
            if len(wHs[xx][yy]) > 0:
                upHs.append(wHs[xx][yy])
                upTp.append(wTp[xx][yy])
                upDm.append(wDm[xx][yy])
                upLWP.append(wLWP[xx][yy])
                upHsSwell.append(swellHs[xx][yy])
                upTpSwell.append(swellTp[xx][yy])
                upDmSwell.append(swellDm[xx][yy])
                upLWPSwell.append(swellLWP[xx][yy])

                upHsWindsea.append(windseaHs[xx][yy])
                upTpWindsea.append(windseaTp[xx][yy])
                upDmWindsea.append(windseaDm[xx][yy])
                upLWPWindsea.append(windseaLWP[xx][yy])

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

stayHsSwellArray = np.concatenate(stayHsSwell).ravel()
upHsSwellArray = np.concatenate(upHsSwell).ravel()
downHsSwellArray = np.concatenate(downHsSwell).ravel()
stayTpSwellArray = np.concatenate(stayTpSwell).ravel()
upTpSwellArray = np.concatenate(upTpSwell).ravel()
downTpSwellArray = np.concatenate(downTpSwell).ravel()
stayDmSwellArray = np.concatenate(stayDmSwell).ravel()
upDmSwellArray = np.concatenate(upDmSwell).ravel()
downDmSwellArray = np.concatenate(downDmSwell).ravel()
stayLWPSwellArray = np.concatenate(stayLWPSwell).ravel()
upLWPSwellArray = np.concatenate(upLWPSwell).ravel()
downLWPSwellArray = np.concatenate(downLWPSwell).ravel()

stayHsWindseaArray = np.concatenate(stayHsWindsea).ravel()
upHsWindseaArray = np.concatenate(upHsWindsea).ravel()
downHsWindseaArray = np.concatenate(downHsWindsea).ravel()
stayTpWindseaArray = np.concatenate(stayTpWindsea).ravel()
upTpWindseaArray = np.concatenate(upTpWindsea).ravel()
downTpWindseaArray = np.concatenate(downTpWindsea).ravel()
stayDmWindseaArray = np.concatenate(stayDmWindsea).ravel()
upDmWindseaArray = np.concatenate(upDmWindsea).ravel()
downDmWindseaArray = np.concatenate(downDmWindsea).ravel()
stayLWPWindseaArray = np.concatenate(stayLWPWindsea).ravel()
upLWPWindseaArray = np.concatenate(upLWPWindsea).ravel()
downLWPWindseaArray = np.concatenate(downLWPWindsea).ravel()

# Where are the extremes?

stayOver2m = np.where((stayHsArray>2.0))
percentStayOver2m = len(stayOver2m[0])/len(stayHsArray)
downOver2m = np.where((downHsArray>2.0))
percentDownOver2m = len(downOver2m[0])/len(downHsArray)
upOver2m = np.where((upHsArray>2.0))
percentUpOver2m = len(upOver2m[0])/len(upHsArray)



import matplotlib.cm as cm
import matplotlib.colors as mcolors
plt.style.use('dark_background')

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




