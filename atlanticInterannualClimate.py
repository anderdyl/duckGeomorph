import scipy.io as sio
from scipy.io.matlab.mio5_params import mat_struct
import numpy as np
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
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import datetime
from dateutil.relativedelta import relativedelta
import random
import xarray as xr
import geopy.distance
from dipy.segment.metric import Metric
from dipy.segment.metric import ResampleFeature
from dipy.segment.clustering import QuickBundles
from mpl_toolkits.basemap import Basemap


def randomcolor():
    colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0,14)]
    return "#"+color


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


def dt2cal(dt):
    """
    Convert array of datetime64 to a calendar array of year, month, day, hour,
    minute, seconds, microsecond with these quantites indexed on the last axis.

    Parameters
    ----------
    dt : datetime64 array (...)
        numpy.ndarray of datetimes of arbitrary shape

    Returns
    -------
    cal : uint32 array (..., 7)
        calendar array with last axis representing year, month, day, hour,
        minute, second, microsecond
    """

    # allocate output
    out = np.empty(dt.shape + (7,), dtype="u4")
    # decompose calendar floors
    Y, M, D, h, m, s = [dt.astype(f"M8[{x}]") for x in "YMDhms"]
    out[..., 0] = Y + 1970 # Gregorian Year
    out[..., 1] = (M - Y) + 1 # month
    out[..., 2] = (D - M) + 1 # dat
    out[..., 3] = (dt - D).astype("m8[h]") # hour
    out[..., 4] = (dt - h).astype("m8[m]") # minute
    out[..., 5] = (dt - m).astype("m8[s]") # second
    out[..., 6] = (dt - s).astype("m8[us]") # microsecond
    return out

class GPSDistance(Metric):
    """computer the average GPS distance between two streamlines"""
    def __init__(self):
        super(GPSDistance, self).__init__(feature=ResampleFeature(nb_points=20))

    def are_compatible(self, shape1, shape2):
        return len(shape1) == len(shape2)

    # def dist(self,v1,v2):
    #     x = [geopy.distance.vincenty([p[0][0],p[0][1]], [p[1][0],p[1][1]]).km for p in list(zip(v1,v2))]
    #     currD = np.mean(x)
    #     return currD
    def dist(self, v1, v2):
        x = [geopy.distance.distance([p[0][0], p[0][1]], [p[1][0], p[1][1]]).kilometers for p in list(zip(v1, v2))]
        currD = np.mean(x)
        return currD

#DWT = ReadMatfile('/media/dylananderson/Elements1/NC_climate/Nags_Head_DWTs_49_w20minDates_2degree_plus5TCs3dayafterentering2.mat')
#DWT = ReadMatfile('/media/dylananderson/Elements1/NC_climate/Nags_Head_DWTS_25_w50minDates_plus5TCs_goodorder.mat')
DWT = ReadMatfile('/media/dylananderson/Elements1/NC_climate/Nags_Head_DWTs_25_w50minDates_2degree_plus5TCs6daysAroundEntering.mat')
PCAmat = ReadMatfile('/media/dylananderson/Elements1/NC_climate/Nags_Head_SLPs_2degree_memory_July2020.mat')
SLPs = ReadMatfile('/media/dylananderson/Elements1/NC_climate/NorthAtlanticSLPs_June2021.mat')
numDWTs = 30
mycols = ReadMatfile('/media/dylananderson/Elements1/shusin6_contents/codes/mycmap_col.mat')
mycmap = mycols['mycmap_col']
# Need to get the dates for the bmus into the correct format (datetime)
def datevec2datetime(d_vec):
    '''
    Returns datetime list from a datevec matrix
    d_vec = [[y1 m1 d1 H1 M1],[y2 ,2 d2 H2 M2],..]
    '''
    return [datetime.datetime(d[0], d[1], d[2], d[3], d[4]) for d in d_vec]

dwtDates = np.array(datevec2datetime(DWT['DWT']['dates']))
dwtBmus = DWT['DWT']['bmus']
dwtBMUS = dwtBmus
order = DWT['DWT']['order']

import matplotlib.cm as cm

etcolors = cm.rainbow(np.linspace(0, 1, numDWTs-5))
tccolors = np.flipud(cm.gray(np.linspace(0,1,6)))

dwtcolors = np.vstack((etcolors,tccolors[1:,:]))





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
sea_nodes = []
for qq in range(lenXB-1):
    sea_nodes.append(np.where((XR == X_B[qq]) & (YR == Y_B[qq])))
    #print(indexTest[0])


flat_list = [item for sublist in sea_nodes for item in sublist]


plt.style.use('default')
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
    m.drawcoastlines()
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
#s_map2.set_array(colorparam)
fig2.subplots_adjust(right=0.85)
cbar_ax2 = fig2.add_axes([0.91, 0.15, 0.02, 0.7])
cbar2 = fig2.colorbar(s_map2, cax=cbar_ax2)
cbar2.set_label('slp anom (mbar)')

plt.show()



DailyPCs = DWT['PCA']['PC']
DailySortedBmus = DWT['DWT']['bmus']-1
DailyTime = DWT['DWT']['dates']
DailyDates = np.array(datevec2datetime(DWT['DWT']['dates']))
DailyDatesMatrix = np.reshape(np.ravel(DailyTime),(len(DailyPCs)+1,6))
#DailyPCs = xr.open_dataset(os.path.join(data_dir, "kma.nc"))

correctedBmus = np.nan*np.ones((np.shape(DailySortedBmus)))
for i in range(30):
    indDWTbmus = order[i]-1
    whereInd = np.where((DailySortedBmus == i))[0]
    correctedBmus[whereInd] = indDWTbmus


s1 = np.full([len(np.unique(DailySortedBmus)), len(np.unique(DailyDatesMatrix[:,0])) - 1], np.nan)

# June/June #!!!!
for i in range(len(np.unique(DailyDatesMatrix[:,0])) - 1):
    s = np.where((DailyDatesMatrix[:,0] == np.unique(DailyDatesMatrix[:,0])[i]) & (DailyDatesMatrix[:,1] == 6))
    ss = np.where((DailyDatesMatrix[:,0] == np.unique(DailyDatesMatrix[:,0])[i] + 1) & (DailyDatesMatrix[:,1] == 5))
    for j in range(len(np.unique(correctedBmus))):
        s1[j, i] = len(np.where(correctedBmus[s[0][0]:ss[0][-1]] == j)[0])


PC1 = np.full([len(np.unique(DailySortedBmus)), len(np.unique(DailyDatesMatrix[:,0])) - 1], np.nan)
PC2 = np.full([len(np.unique(DailySortedBmus)), len(np.unique(DailyDatesMatrix[:,0])) - 1], np.nan)
PC3 = np.full([len(np.unique(DailySortedBmus)), len(np.unique(DailyDatesMatrix[:,0])) - 1], np.nan)
PC4 = np.full([len(np.unique(DailySortedBmus)), len(np.unique(DailyDatesMatrix[:,0])) - 1], np.nan)
PC5 = np.full([len(np.unique(DailySortedBmus)), len(np.unique(DailyDatesMatrix[:,0])) - 1], np.nan)

# June/June #!!!!
for i in range(len(np.unique(DailyDatesMatrix[:,0])) - 1):
    s = np.where((DailyDatesMatrix[:,0] == np.unique(DailyDatesMatrix[:,0])[i]) & (DailyDatesMatrix[:,1] == 6))
    ss = np.where((DailyDatesMatrix[:,0] == np.unique(DailyDatesMatrix[:,0])[i] + 1) & (DailyDatesMatrix[:,1] == 5))
    for j in range(len(np.unique(DailySortedBmus))):
        indDWT = order[j]-1                                 # need the order of this DWT sample
        yearlyBmus = DailySortedBmus[s[0][0]:ss[0][-1]]     # lets get the bmus for this year
        yearlyPCs = DailyPCs[s[0][0]:ss[0][-1],:]           # lets get the PCs for this year
        indBmus = np.where(yearlyBmus == indDWT)[0]         # which of those DWTs do we care about in this iteration
        # PC1[indDWT,i] = np.nanmean(DailyPCs[0,indBmus])
        # PC2[indDWT,i] = np.nanmean(DailyPCs[1,indBmus])
        # PC3[indDWT,i] = np.nanmean(DailyPCs[2,indBmus])
        PC1[indDWT,i] = np.nansum(yearlyPCs[indBmus,0])
        PC2[indDWT,i] = np.nansum(yearlyPCs[indBmus,1])
        PC3[indDWT,i] = np.nansum(yearlyPCs[indBmus,2])
        PC4[indDWT,i] = np.nansum(yearlyPCs[indBmus,3])
        PC5[indDWT,i] = np.nansum(yearlyPCs[indBmus,4])


npercent = DWT['Npercent']
tempPC1 = np.nansum(PC1,axis=0)
tempPC2 = np.nansum(PC2,axis=0)
tempPC3 = np.nansum(PC3,axis=0)
tempPC4 = np.nansum(PC4,axis=0)
tempPC5 = np.nansum(PC5,axis=0)

normPC1 = np.divide(tempPC1,np.nanmax(tempPC1))*npercent[0]
normPC2 = np.divide(tempPC2,np.nanmax(tempPC2))*npercent[1]
normPC3 = np.divide(tempPC3,np.nanmax(tempPC3))*npercent[2]
normPC4 = np.divide(tempPC4,np.nanmax(tempPC4))*npercent[3]
normPC5 = np.divide(tempPC5,np.nanmax(tempPC5))*npercent[4]

pcAggregates = np.full((len(normPC1),5),np.nan)
pcAggregates[:,0] = normPC1
pcAggregates[:,1] = normPC2
pcAggregates[:,2] = normPC3
pcAggregates[:,3] = normPC4
pcAggregates[:,4] = normPC5

n_clusters = 6

kmeans = KMeans(n_clusters, init='k-means++', random_state=100)  # 80

n_components = 5  # !!!!
data = pcAggregates#[:, 0:n_components]

#    data1=data/np.std(data,axis=0)

awt_bmus = kmeans.fit_predict(data)


fig = plt.figure(figsize=[14, 9])
gs2 = gridspec.GridSpec(n_components + 1, 1)
for nn in range(n_components):
    ax2 = fig.add_subplot(gs2[nn])
    ax2.plot(np.unique(DailyDatesMatrix[:,0])[:-1], pcAggregates[:, nn], 'k.-', linewidth=1.8, markersize=8)

    ax2.set_ylabel('PC-' + str(nn + 1), fontsize=13)
    ax2.grid('minor')
    ax2.set_xticklabels([])

ax2 = fig.add_subplot(gs2[nn + 1])
ax2.plot(np.unique(DailyDatesMatrix[:,0])[:-1], awt_bmus + 1, 'k.:', linewidth=1.8, markersize=10, color='grey')


import csv
with open('atlanticMultiDecadalOscillation.txt', 'r') as fd:
    c = 0
    dataAMO = list()
    for line in fd:
        splitLine = line.split()
        print(splitLine[1:])
        for t in splitLine[1:]:
            dataAMO.append(float(t))

amo = np.asarray(dataAMO)
amo = amo[0:-7]
dt = datetime.date(1856, 1, 1)
end = datetime.date(2021, 6, 1)
#step = datetime.timedelta(months=1)
step = relativedelta(months=1)
amoTime = []
while dt < end:
    amoTime.append(dt)#.strftime('%Y-%m-%d'))
    dt += step



with open('NAO2021.txt', 'r') as fd:
    c = 0
    dataNAO = list()
    for line in fd:
        splitLine = line.split(',')
        secondSplit = splitLine[1].split('/')
        dataNAO.append(float(secondSplit[0]))
nao = np.asarray(dataNAO)

dt = datetime.date(1950, 1, 1)
end = datetime.date(2021, 6, 1)
#step = datetime.timedelta(months=1)
step = relativedelta(months=1)
naoTime = []
while dt < end:
    naoTime.append(dt)#.strftime('%Y-%m-%d'))
    dt += step


with open('aceNorthAtlanticHurricanes.csv', 'r') as fd:
    c = 0
    dataACE = list()
    for line in fd:
        splitLine = line.split(',')
        secondSplit = splitLine[3].split('/')
        print(secondSplit[0])
        dataACE.append(float(secondSplit[0]))
ace = np.asarray(dataACE)
dt = datetime.date(1851, 1, 1)
end = datetime.date(2020, 5, 1)
#step = datetime.timedelta(months=1)
step = relativedelta(years=1)
aceTime = []
while dt < end:
    aceTime.append(dt)#.strftime('%Y-%m-%d'))
    dt += step


dataMJO = ReadMatfile('/media/dylananderson/Elements1/SERDP/Data/MJO/mjo_australia_2021.mat')
mjoPhase = dataMJO['phase']
dt = datetime.date(1974, 6, 1)
end = datetime.date(2021, 6, 18)
step = relativedelta(days=1)
mjoTime = []
while dt < end:
    mjoTime.append(dt)#.strftime('%Y-%m-%d'))
    dt += step


figClimate = plt.figure()
ax1Cl = plt.subplot2grid((4,1),(0,0),rowspan=1,colspan=1)
ax1Cl.plot(mjoTime,mjoPhase,'.')
ax1Cl.set_xlim([datetime.date(1979,1,1),datetime.date(2021,5,1)])
ax1Cl.set_ylim([0,9])
ax1Cl.set_ylabel('MJO')
ax2Cl = plt.subplot2grid((4,1),(1,0),rowspan=1,colspan=1)
ax2Cl.plot(aceTime,ace)
ax2Cl.set_xlim([datetime.date(1979,1,1),datetime.date(2021,5,1)])
ax2Cl.set_ylabel('ACE')
ax3Cl = plt.subplot2grid((4,1),(2,0),rowspan=1,colspan=1)
ax3Cl.plot(naoTime,nao)
ax3Cl.set_xlim([datetime.date(1979,1,1),datetime.date(2021,5,1)])
ax3Cl.set_ylabel('NAO')
ax4Cl = plt.subplot2grid((4,1),(3,0),rowspan=1,colspan=1)
ax4Cl.plot(amoTime,amo)
ax4Cl.set_xlim([datetime.date(1979,1,1),datetime.date(2021,5,1)])
ax4Cl.set_ylabel('AMO')



file = '/media/dylananderson/Elements1/SERDP/Data/TC_tracks/IBTrACS.NA.v04r00.nc'
data = xr.open_dataset(file)
#merged = data['merged'].values
TCtime = data['time'].values
TClon = data['lon'].values
TClat = data['lat'].values
TCpres = data['wmo_pres']
TCwind = data['wmo_wind']

indexTC = np.arange(1,len(TCtime))

streams = []
streamTime = []
streamWind = []
streamPres = []
for hh in range(len(TCtime)):
    indexReal = np.where(TClat[hh,:] > 0)
    if len(indexReal[0])>1:
        lat_lng_data = np.c_[TClat[hh,indexReal[0]], TClon[hh,indexReal[0]]]
        streams.append(lat_lng_data)
        streamTime.append(TCtime[hh,indexReal[0]])
        streamPres.append(TCpres[hh,indexReal[0]])
        streamWind.append(TCwind[hh,indexReal[0]])

metric = GPSDistance()
qb = QuickBundles(threshold=2750,metric=metric)
clusters = qb.cluster(streams)
print("Nb. clusters:",len(clusters))

fig = plt.figure()
#m = Basemap(projection='merc', llcrnrlat=0, urcrnrlat=60, llcrnrlon=250, urcrnrlon=-10, lat_ts=10,resolution='c')
#m = Basemap(projection='ortho',lat_0=0, lon_0=0)
# #m.fillcontinents(color=dwtcolors[qq])
plotIndy=0
plotIndx=0
for clustersIndex in range(6):
    p1 = plt.subplot2grid((2, 3), (plotIndx, plotIndy), rowspan=1, colspan=1)
    m = Basemap(llcrnrlon=-120.7, llcrnrlat=0., urcrnrlon=-10.1, urcrnrlat=60, projection='merc', lat_1=30., lat_2=60.,
                lat_0=34.83158, lon_0=-98.)

    color = randomcolor()
    for i in clusters[clustersIndex].indices:
        cx, cy = m(streams[i][:,1], streams[i][:,0])  # convert to map projection coordinate
        # cx = streams[i][:,1]
        # cy = streams[i][:,0]
        m.plot(cx, cy, marker=None, color=color)
    m.drawcoastlines()
    if plotIndy < 2:
        plotIndy = plotIndy + 1
    else:
        plotIndy = 0
        plotIndx = plotIndx + 1
for i in clusters[6].indices:
    cx, cy = m(streams[i][:,1], streams[i][:,0])  # convert to map projection coordinate
        # cx = streams[i][:,1]
        # cy = streams[i][:,0]
    m.plot(cx, cy, marker=None, color=color)






fig = plt.figure()
#m = Basemap(projection='merc', llcrnrlat=0, urcrnrlat=60, llcrnrlon=250, urcrnrlon=-10, lat_ts=10,resolution='c')
#m = Basemap(projection='ortho',lat_0=0, lon_0=0)
# #m.fillcontinents(color=dwtcolors[qq])
plotIndy=0
plotIndx=0
for clustersIndex in range(6):
    p1 = plt.subplot2grid((2, 3), (plotIndx, plotIndy), rowspan=1, colspan=1)
    m = Basemap(llcrnrlon=-120.7, llcrnrlat=0., urcrnrlon=-10.1, urcrnrlat=60, projection='merc', lat_1=30., lat_2=60.,
                lat_0=34.83158, lon_0=-98.)

    color = randomcolor()
    for i in clusters[clustersIndex].indices:
        cx, cy = m(streams[i][:,1], streams[i][:,0])  # convert to map projection coordinate
        m.plot(cx[0], cy[0], marker='.', color=color)
    m.drawcoastlines()
    if plotIndy < 2:
        plotIndy = plotIndy + 1
    else:
        plotIndy = 0
        plotIndx = plotIndx + 1
for i in clusters[6].indices:
    cx, cy = m(streams[i][:, 1], streams[i][:, 0])  # convert to map projection coordinate
    m.plot(cx[0], cy[0], marker='.', color=color)









fig = plt.figure()

p1 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
m = Basemap(llcrnrlon=-120.7, llcrnrlat=0., urcrnrlon=-10.1, urcrnrlat=60, projection='merc', lat_1=30., lat_2=60.,lat_0=34.83158, lon_0=-98.)
for clustersIndex in range(6):
    color = randomcolor()
    cx, cy = m(clusters.centroids[clustersIndex][:,1], clusters.centroids[clustersIndex][:,0])  # convert to map projection coordinate
    m.plot(cx, cy, marker=None, color=color)
# cx, cy = m(clusters.centroids[6][:,1], clusters.centroids[6][:,0])  # convert to map projection coordinate
# m.plot(cx, cy, marker=None, color=color)
m.drawcoastlines()

cluster1minTime = []
cluster1Time = []
for i in clusters[0].indices:
    arrayYear = dt2cal(streamTime[i][0])[0]
    if arrayYear > 1978:
        arrayTime = [dt2cal(dt) for dt in streamTime[i]]
        presTemp = streamPres[i].values
        minPresind = np.where((np.nanmin(presTemp)==presTemp))
        if len(minPresind[0]) > 0:
            datesToAppendMin = arrayTime[minPresind[0][0]]
            print('chose {}, started {}, ended {}'.format(datesToAppendMin,arrayTime[0],arrayTime[-1]))
            cluster1minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
        else:
            windTemp = streamWind[i].values
            minWindind = np.where((np.nanmax(windTemp)==windTemp))
            if len(minWindind[0]) > 0:
                datesToAppendMin = arrayTime[minWindind[0][0]]
                cluster1minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
            else:
                midStorm = np.round((len(arrayTime)/2))
                datesToAppendMin = arrayTime[int(midStorm)]
                cluster1minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))

        arrayDay = np.asarray([list((hh[0],hh[1],hh[2])) for hh in arrayTime])
        uniqueDay = np.unique(arrayDay[:,2],return_index=True)
        datesToAppend = arrayDay[uniqueDay[1],:]
        cluster1Time.append(datesToAppend)


cluster2Time = []
cluster2minTime = []
for i in clusters[1].indices:
    arrayYear = dt2cal(streamTime[i][0])[0]
    if arrayYear > 1978:
        arrayTime = [dt2cal(dt) for dt in streamTime[i]]
        presTemp = streamPres[i].values
        minPresind = np.where((np.nanmin(presTemp)==presTemp))
        if len(minPresind[0]) > 0:
            datesToAppendMin = arrayTime[minPresind[0][0]]
            cluster2minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
        else:
            windTemp = streamWind[i].values
            minWindind = np.where((np.nanmax(windTemp)==windTemp))
            if len(minWindind[0]) > 0:
                datesToAppendMin = arrayTime[minWindind[0][0]]
                cluster2minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
            else:
                midStorm = np.round((len(arrayTime)/2))
                datesToAppendMin = arrayTime[int(midStorm)]
                cluster2minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))

        arrayDay = np.asarray([list((hh[0],hh[1],hh[2])) for hh in arrayTime])
        uniqueDay = np.unique(arrayDay[:,2],return_index=True)
        datesToAppend = arrayDay[uniqueDay[1],:]
        cluster2Time.append(datesToAppend)


cluster3Time = []
cluster3minTime = []
for i in clusters[2].indices:
    arrayYear = dt2cal(streamTime[i][0])[0]
    if arrayYear > 1978:
        arrayTime = [dt2cal(dt) for dt in streamTime[i]]
        presTemp = streamPres[i].values
        minPresind = np.where((np.nanmin(presTemp)==presTemp))
        if len(minPresind[0]) > 0:
            datesToAppendMin = arrayTime[minPresind[0][0]]
            cluster3minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
        else:
            windTemp = streamWind[i].values
            minWindind = np.where((np.nanmax(windTemp)==windTemp))
            if len(minWindind[0]) > 0:
                datesToAppendMin = arrayTime[minWindind[0][0]]
                cluster3minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
            else:
                midStorm = np.round((len(arrayTime)/2))
                datesToAppendMin = arrayTime[int(midStorm)]
                cluster3minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))

        arrayDay = np.asarray([list((hh[0],hh[1],hh[2])) for hh in arrayTime])
        uniqueDay = np.unique(arrayDay[:,2],return_index=True)
        datesToAppend = arrayDay[uniqueDay[1],:]
        cluster3Time.append(datesToAppend)


cluster4Time = []
cluster4minTime = []
for i in clusters[3].indices:
    arrayYear = dt2cal(streamTime[i][0])[0]
    if arrayYear > 1978:
        arrayTime = [dt2cal(dt) for dt in streamTime[i]]
        presTemp = streamPres[i].values
        minPresind = np.where((np.nanmin(presTemp)==presTemp))
        if len(minPresind[0]) > 0:
            datesToAppendMin = arrayTime[minPresind[0][0]]
            cluster4minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
        else:
            windTemp = streamWind[i].values
            minWindind = np.where((np.nanmax(windTemp)==windTemp))
            if len(minWindind[0]) > 0:
                datesToAppendMin = arrayTime[minWindind[0][0]]
                cluster4minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
            else:
                midStorm = np.round((len(arrayTime)/2))
                datesToAppendMin = arrayTime[int(midStorm)]
                cluster4minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))

        arrayDay = np.asarray([list((hh[0],hh[1],hh[2])) for hh in arrayTime])
        uniqueDay = np.unique(arrayDay[:,2],return_index=True)
        datesToAppend = arrayDay[uniqueDay[1],:]
        cluster4Time.append(datesToAppend)

cluster5Time = []
cluster5minTime = []
for i in clusters[4].indices:
    arrayYear = dt2cal(streamTime[i][0])[0]
    if arrayYear > 1978:
        arrayTime = [dt2cal(dt) for dt in streamTime[i]]
        presTemp = streamPres[i].values
        minPresind = np.where((np.nanmin(presTemp)==presTemp))
        if len(minPresind[0]) > 0:
            datesToAppendMin = arrayTime[minPresind[0][0]]
            cluster5minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
        else:
            windTemp = streamWind[i].values
            minWindind = np.where((np.nanmax(windTemp)==windTemp))
            if len(minWindind[0]) > 0:
                datesToAppendMin = arrayTime[minWindind[0][0]]
                cluster5minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
            else:
                midStorm = np.round((len(arrayTime)/2))
                datesToAppendMin = arrayTime[int(midStorm)]
                cluster5minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))

        arrayDay = np.asarray([list((hh[0],hh[1],hh[2])) for hh in arrayTime])
        uniqueDay = np.unique(arrayDay[:,2],return_index=True)
        datesToAppend = arrayDay[uniqueDay[1],:]
        cluster5Time.append(datesToAppend)

cluster6Time = []
cluster6minTime = []
for i in clusters[5].indices:
    arrayYear = dt2cal(streamTime[i][0])[0]
    if arrayYear > 1978:
        arrayTime = [dt2cal(dt) for dt in streamTime[i]]
        presTemp = streamPres[i].values
        minPresind = np.where((np.nanmin(presTemp)==presTemp))
        if len(minPresind[0]) > 0:
            datesToAppendMin = arrayTime[minPresind[0][0]]
            cluster6minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
        else:
            windTemp = streamWind[i].values
            minWindind = np.where((np.nanmax(windTemp)==windTemp))
            if len(minWindind[0]) > 0:
                datesToAppendMin = arrayTime[minWindind[0][0]]
                cluster6minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))
            else:
                midStorm = np.round((len(arrayTime)/2))
                datesToAppendMin = arrayTime[int(midStorm)]
                cluster6minTime.append(np.array([datesToAppendMin[0],datesToAppendMin[1],datesToAppendMin[2]]))

        arrayDay = np.asarray([list((hh[0],hh[1],hh[2])) for hh in arrayTime])
        uniqueDay = np.unique(arrayDay[:,2],return_index=True)
        datesToAppend = arrayDay[uniqueDay[1],:]
        cluster6Time.append(datesToAppend)

X_in = SLPs['X_in']
Y_in = SLPs['Y_in']
SLP = SLPs['slp_mem']
SLPtime = SLPs['time']
#
# Dindices1 = []
# for hh in range(len(cluster1Time)):
#     times = cluster1Time[hh]
#     for qq in range(len(times)):
#         dIndex = np.where((times[qq][0]==SLPtime[:,0]) & (times[qq][1]==SLPtime[:,1]) & (times[qq][2]==SLPtime[:,2]))
#         Dindices1.append(dIndex[0][0])
# tc1Dindices = np.asarray(Dindices1)
# cluster1SLPIndex = np.unique(tc1Dindices)
# cluster1SLPs = np.nanmean(SLP[cluster1SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
# Dindices2 = []
# for hh in range(len(cluster2Time)):
#     times = cluster2Time[hh]
#     for qq in range(len(times)):
#         dIndex = np.where((times[qq][0]==SLPtime[:,0]) & (times[qq][1]==SLPtime[:,1]) & (times[qq][2]==SLPtime[:,2]))
#         Dindices2.append(dIndex[0][0])
# tc2Dindices = np.asarray(Dindices2)
# cluster2SLPIndex = np.unique(tc2Dindices)
# cluster2SLPs = np.nanmean(SLP[cluster2SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
# Dindices3 = []
# for hh in range(len(cluster3Time)):
#     times = cluster3Time[hh]
#     for qq in range(len(times)):
#         dIndex = np.where((times[qq][0]==SLPtime[:,0]) & (times[qq][1]==SLPtime[:,1]) & (times[qq][2]==SLPtime[:,2]))
#         Dindices3.append(dIndex[0][0])
# tc3Dindices = np.asarray(Dindices3)
# cluster3SLPIndex = np.unique(tc3Dindices)
# cluster3SLPs = np.nanmean(SLP[cluster3SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
# Dindices4 = []
# for hh in range(len(cluster4Time)):
#     times = cluster4Time[hh]
#     for qq in range(len(times)):
#         dIndex = np.where((times[qq][0]==SLPtime[:,0]) & (times[qq][1]==SLPtime[:,1]) & (times[qq][2]==SLPtime[:,2]))
#         Dindices4.append(dIndex[0][0])
# tc4Dindices = np.asarray(Dindices4)
# cluster4SLPIndex = np.unique(tc4Dindices)
# cluster4SLPs = np.nanmean(SLP[cluster4SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
# Dindices5 = []
# for hh in range(len(cluster5Time)):
#     times = cluster5Time[hh]
#     for qq in range(len(times)):
#         dIndex = np.where((times[qq][0]==SLPtime[:,0]) & (times[qq][1]==SLPtime[:,1]) & (times[qq][2]==SLPtime[:,2]))
#         Dindices5.append(dIndex[0][0])
# tc5Dindices = np.asarray(Dindices5)
# cluster5SLPIndex = np.unique(tc5Dindices)
# cluster5SLPs = np.nanmean(SLP[cluster5SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
# Dindices6 = []
# for hh in range(len(cluster6Time)):
#     times = cluster6Time[hh]
#     for qq in range(len(times)):
#         dIndex = np.where((times[qq][0]==SLPtime[:,0]) & (times[qq][1]==SLPtime[:,1]) & (times[qq][2]==SLPtime[:,2]))
#         Dindices6.append(dIndex[0][0])
# tc6Dindices = np.asarray(Dindices6)
# cluster6SLPIndex = np.unique(tc6Dindices)
# cluster6SLPs = np.nanmean(SLP[cluster6SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100

Dindices1 = []
for hh in range(len(cluster1minTime)):
    times = cluster1minTime[hh]
    dIndex = np.where((times[0]==SLPtime[:,0]) & ((times[1])==SLPtime[:,1]) & ((int(times[2]))==SLPtime[:,2]))
    Dindices1.append(np.arange((dIndex[0][0]-1),(dIndex[0][0]+2)))
tc1Dindices = np.asarray(Dindices1)
cluster1SLPIndex = np.unique(tc1Dindices)
cluster1SLPs = np.nanmean(SLP[cluster1SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
Dindices2 = []
for hh in range(len(cluster2minTime)):
    times = cluster2minTime[hh]
    dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
    Dindices2.append(np.arange((dIndex[0][0]-1),(dIndex[0][0]+2)))
tc2Dindices = np.asarray(Dindices2)
cluster2SLPIndex = np.unique(tc2Dindices)
cluster2SLPs = np.nanmean(SLP[cluster2SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100

Dindices3 = []
for hh in range(len(cluster3minTime)):
    times = cluster3minTime[hh]
    dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
    Dindices3.append(np.arange((dIndex[0][0]-1),(dIndex[0][0]+2)))
tc3Dindices = np.asarray(Dindices3)
cluster3SLPIndex = np.unique(tc3Dindices)
cluster3SLPs = np.nanmean(SLP[cluster3SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100

Dindices4 = []
for hh in range(len(cluster4minTime)):
    times = cluster4minTime[hh]
    dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
    Dindices4.append(np.arange((dIndex[0][0]-1),(dIndex[0][0]+2)))
tc4Dindices = np.asarray(Dindices4)
cluster4SLPIndex = np.unique(tc4Dindices)
cluster4SLPs = np.nanmean(SLP[cluster4SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100

Dindices5 = []
for hh in range(len(cluster5minTime)):
    times = cluster5minTime[hh]
    dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
    Dindices5.append(np.arange((dIndex[0][0]-1),(dIndex[0][0]+2)))
tc5Dindices = np.asarray(Dindices5)
cluster5SLPIndex = np.unique(tc5Dindices)
cluster5SLPs = np.nanmean(SLP[cluster5SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100

Dindices6 = []
for hh in range(len(cluster6minTime)):
    times = cluster6minTime[hh]
    dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
    Dindices6.append(np.arange((dIndex[0][0]-1),(dIndex[0][0]+2)))
tc6Dindices = np.asarray(Dindices6)
cluster6SLPIndex = np.unique(tc6Dindices)
cluster6SLPs = np.nanmean(SLP[cluster6SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100


#
#
# Dindices1 = []
# for hh in range(len(cluster1minTime)):
#     times = cluster1minTime[hh]
#     dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
#     Dindices1.append(dIndex[0][0])
# tc1Dindices = np.asarray(Dindices1)
# cluster1SLPIndex = np.unique(tc1Dindices)
# cluster1SLPs = np.nanmean(SLP[cluster1SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
# Dindices2 = []
# for hh in range(len(cluster2minTime)):
#     times = cluster2minTime[hh]
#     dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
#     Dindices2.append(dIndex[0][0])
# tc2Dindices = np.asarray(Dindices2)
# cluster2SLPIndex = np.unique(tc2Dindices)
# cluster2SLPs = np.nanmean(SLP[cluster2SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
#
# Dindices3 = []
# for hh in range(len(cluster3minTime)):
#     times = cluster3minTime[hh]
#     dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
#     Dindices3.append(dIndex[0][0])
# tc3Dindices = np.asarray(Dindices3)
# cluster3SLPIndex = np.unique(tc3Dindices)
# cluster3SLPs = np.nanmean(SLP[cluster3SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
#
# Dindices4 = []
# for hh in range(len(cluster4minTime)):
#     times = cluster4minTime[hh]
#     dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
#     Dindices4.append(dIndex[0][0])
# tc4Dindices = np.asarray(Dindices4)
# cluster4SLPIndex = np.unique(tc4Dindices)
# cluster4SLPs = np.nanmean(SLP[cluster4SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
#
# Dindices5 = []
# for hh in range(len(cluster5minTime)):
#     times = cluster5minTime[hh]
#     dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
#     Dindices5.append(dIndex[0][0])
# tc5Dindices = np.asarray(Dindices5)
# cluster5SLPIndex = np.unique(tc5Dindices)
# cluster5SLPs = np.nanmean(SLP[cluster5SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100
#
# Dindices6 = []
# for hh in range(len(cluster6minTime)):
#     times = cluster6minTime[hh]
#     dIndex = np.where((times[0]==SLPtime[:,0]) & (times[1]==SLPtime[:,1]) & (times[2]==SLPtime[:,2]))
#     Dindices6.append(dIndex[0][0])
# tc6Dindices = np.asarray(Dindices6)
# cluster6SLPIndex = np.unique(tc6Dindices)
# cluster6SLPs = np.nanmean(SLP[cluster6SLPIndex,:],axis=0).reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100


#temp = np.nanmean(SLP_C[num_index, :], axis=1) / 100 - np.nanmean(SLP_C, axis=0) / 100

wt = SLP[2000,:].reshape(73,43)/100 - np.nanmean(SLP, axis=0).reshape(73,43) / 100

fig3 = plt.figure()
clevels = np.arange(-20,20,1)
cx,cy =m(X_in,Y_in)  # convert to map projection coordinate

ax1 = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
m = Basemap(projection='merc',llcrnrlat=-10,urcrnrlat=70,llcrnrlon=245,urcrnrlon=370,lat_ts=10,resolution='c')
#m.fillcontinents(color=dwtcolors[qq])
m.drawcoastlines()
CS = m.contourf(cx,cy,cluster1SLPs.T,clevels,vmin=-5,vmax=7,cmap=cm.RdBu_r,shading='gouraud')

ax2 = plt.subplot2grid((2,3),(0,1),rowspan=1,colspan=1)
m = Basemap(projection='merc',llcrnrlat=-10,urcrnrlat=70,llcrnrlon=245,urcrnrlon=370,lat_ts=10,resolution='c')
m.drawcoastlines()
CS = m.contourf(cx,cy,cluster2SLPs.T,clevels,vmin=-5,vmax=5,cmap=cm.RdBu_r,shading='gouraud')

ax3 = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
m = Basemap(projection='merc',llcrnrlat=-10,urcrnrlat=70,llcrnrlon=245,urcrnrlon=370,lat_ts=10,resolution='c')
m.drawcoastlines()
CS = m.contourf(cx,cy,cluster3SLPs.T,clevels,vmin=-5,vmax=5,cmap=cm.RdBu_r,shading='gouraud')

ax4 = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
m = Basemap(projection='merc',llcrnrlat=-10,urcrnrlat=70,llcrnrlon=245,urcrnrlon=370,lat_ts=10,resolution='c')
m.drawcoastlines()
CS = m.contourf(cx,cy,cluster4SLPs.T,clevels,vmin=-5,vmax=5,cmap=cm.RdBu_r,shading='gouraud')

ax5 = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
m = Basemap(projection='merc',llcrnrlat=-10,urcrnrlat=70,llcrnrlon=245,urcrnrlon=370,lat_ts=10,resolution='c')
m.drawcoastlines()
CS = m.contourf(cx,cy,cluster5SLPs.T,clevels,vmin=-5,vmax=5,cmap=cm.RdBu_r,shading='gouraud')

ax6 = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
m = Basemap(projection='merc',llcrnrlat=-10,urcrnrlat=70,llcrnrlon=245,urcrnrlon=370,lat_ts=10,resolution='c')
m.drawcoastlines()
CS = m.contourf(cx,cy,cluster6SLPs.T,clevels,vmin=-5,vmax=5,cmap=cm.RdBu_r,shading='gouraud')

# colormap = cm.RdBu_r
# normalize = mcolors.Normalize(vmin=-7.5, vmax=7.5)
#
# s_map2 = cm.ScalarMappable(norm=normalize, cmap=colormap)
# #s_map2.set_array(colorparam)
# fig3.subplots_adjust(right=0.85)
# cbar_ax2 = fig3.add_axes([0.91, 0.15, 0.02, 0.7])
# cbar2 = fig3.colorbar(s_map2, cax=cbar_ax2)
# cbar2.set_label('slp anom (mbar)')
#
# plt.show()
mjoArrayTime = np.nan*np.ones((len(mjoTime),3))
c = 0
for qq in mjoTime:
    mjoArrayTime[c,0] = qq.year
    mjoArrayTime[c,1] = qq.month
    mjoArrayTime[c,2] = qq.day
    c += 1


mjoDindices1 = []
for hh in range(len(cluster1Time)):
    times = cluster1Time[hh]
    for qq in range(len(times)):
        dIndex = np.where((times[qq][0]==mjoArrayTime[:,0]) & (times[qq][1]==mjoArrayTime[:,1]) & (times[qq][2]==mjoArrayTime[:,2]))
        mjoDindices1.append(dIndex[0][0])
mjo1Dindices = np.asarray(mjoDindices1)
cluster1MJOIndex = np.unique(mjo1Dindices)
cluster1MJO = mjoPhase[cluster1MJOIndex]

mjoDindices2 = []
for hh in range(len(cluster2Time)):
    times = cluster2Time[hh]
    for qq in range(len(times)):
        dIndex = np.where((times[qq][0]==mjoArrayTime[:,0]) & (times[qq][1]==mjoArrayTime[:,1]) & (times[qq][2]==mjoArrayTime[:,2]))
        mjoDindices2.append(dIndex[0][0])
mjo2Dindices = np.asarray(mjoDindices2)
cluster2MJOIndex = np.unique(mjo2Dindices)
cluster2MJO = mjoPhase[cluster2MJOIndex]

mjoDindices3 = []
for hh in range(len(cluster3Time)):
    times = cluster3Time[hh]
    for qq in range(len(times)):
        dIndex = np.where((times[qq][0]==mjoArrayTime[:,0]) & (times[qq][1]==mjoArrayTime[:,1]) & (times[qq][2]==mjoArrayTime[:,2]))
        mjoDindices3.append(dIndex[0][0])
mjo3Dindices = np.asarray(mjoDindices3)
cluster3MJOIndex = np.unique(mjo3Dindices)
cluster3MJO = mjoPhase[cluster3MJOIndex]

mjoDindices4 = []
for hh in range(len(cluster4Time)):
    times = cluster4Time[hh]
    for qq in range(len(times)):
        dIndex = np.where((times[qq][0]==mjoArrayTime[:,0]) & (times[qq][1]==mjoArrayTime[:,1]) & (times[qq][2]==mjoArrayTime[:,2]))
        mjoDindices4.append(dIndex[0][0])
mjo4Dindices = np.asarray(mjoDindices4)
cluster4MJOIndex = np.unique(mjo4Dindices)
cluster4MJO = mjoPhase[cluster4MJOIndex]

mjoDindices5 = []
for hh in range(len(cluster5Time)):
    times = cluster5Time[hh]
    for qq in range(len(times)):
        dIndex = np.where((times[qq][0]==mjoArrayTime[:,0]) & (times[qq][1]==mjoArrayTime[:,1]) & (times[qq][2]==mjoArrayTime[:,2]))
        mjoDindices5.append(dIndex[0][0])
mjo5Dindices = np.asarray(mjoDindices5)
cluster5MJOIndex = np.unique(mjo5Dindices)
cluster5MJO = mjoPhase[cluster5MJOIndex]

mjoDindices6 = []
for hh in range(len(cluster6Time)):
    times = cluster6Time[hh]
    for qq in range(len(times)):
        dIndex = np.where((times[qq][0]==mjoArrayTime[:,0]) & (times[qq][1]==mjoArrayTime[:,1]) & (times[qq][2]==mjoArrayTime[:,2]))
        mjoDindices6.append(dIndex[0][0])
mjo6Dindices = np.asarray(mjoDindices6)
cluster6MJOIndex = np.unique(mjo6Dindices)
cluster6MJO = mjoPhase[cluster6MJOIndex]



fig4 = plt.figure()
ax1 = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
ax1.hist(cluster1MJO,8)
ax1.set_xlabel('MJO Phase')
ax1.set_title('Cluster 1')
ax2 = plt.subplot2grid((2,3),(0,1),rowspan=1,colspan=1)
ax2.hist(cluster2MJO,8)
ax2.set_xlabel('MJO Phase')
ax2.set_title('Cluster 2')
ax3 = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
ax3.hist(cluster3MJO,8)
ax3.set_xlabel('MJO Phase')
ax3.set_title('Cluster 3')
ax4 = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
ax4.hist(cluster4MJO,8)
ax4.set_xlabel('MJO Phase')
ax4.set_title('Cluster 4')
ax5 = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
ax5.hist(cluster5MJO,8)
ax5.set_xlabel('MJO Phase')
ax5.set_title('Cluster 5')
ax6 = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
ax6.hist(cluster6MJO,8)
ax6.set_xlabel('MJO Phase')
ax6.set_title('Cluster 6')






#import geopandas as gpd





# import similaritymeasures
# def compute_distance_matrix(trajectories, method="Frechet"):
#     """
#     :param method: "Frechet" or "Area"
#     """
#     n = len(trajectories)
#     dist_m = np.zeros((n, n))
#     for i in range(n - 1):
#         p = trajectories[i]
#         for j in range(i + 1, n):
#             q = trajectories[j]
#             if method == "Frechet":
#                 dist_m[i, j] = similaritymeasures.frechet_dist(p, q)
#             else:
#                 dist_m[i, j] = similaritymeasures.area_between_two_curves(p, q)
#             dist_m[j, i] = dist_m[i, j]
#     return dist_m



# Need to pre-process these tracks...

# from sklearn import metrics
# from sklearn.mixture import GMM
# gmm = GMM(n_components=no_of_cluster,  n_iter=1000)
# labels = gmm.fit_predict(TC)




a=np.tile(np.sum(s1,axis=1),(np.shape(s1)[1],1))
a=a.T
s2=s1/a

pca = PCA()
principalComponents = pca.fit_transform(s2.T)
eofs=pca.components_ #(n_components, n_features)
ev=pca.explained_variance_ratio_




x=np.flipud(range(0, 6))


plt.figure()
plt.pcolormesh(np.rot90((eofs[:,0].reshape(5,6))),cmap='RdBu_r',vmin=-0.5, vmax=0.5)

fig = plt.figure(figsize=[14,9])
gs2=gridspec.GridSpec(1,1)
ax2=fig.add_subplot(gs2[0])
ax2.pcolormesh(np.rot90((eofs[:,0].reshape(5,6))),cmap='RdBu_r',vmin=-0.5, vmax=0.5)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_aspect('equal')

n_clusters = 4

kmeans = KMeans(n_clusters, init='k-means++', random_state=100)  # 80

n_components = 3  # !!!!
data = principalComponents[:, 0:n_components]

#    data1=data/np.std(data,axis=0)

awt_bmus = kmeans.fit_predict(data)

fig = plt.figure(figsize=[14, 9])
gs2 = gridspec.GridSpec(n_components + 1, 1)
for nn in range(n_components):
    ax2 = fig.add_subplot(gs2[nn])
    ax2.plot(np.unique(DailyDatesMatrix[:,0])[:-1], principalComponents[:, nn], 'k.-', linewidth=1.8, markersize=8)

    # ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_vso], principalComponents[:, nn][s_vso], '.', markersize=12,
    #          label='V.Str.-EN', color='darkred')
    # ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_so], principalComponents[:, nn][s_so], '.', markersize=12,
    #          label='Str.-EN', color='crimson')
    # ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_mo], principalComponents[:, nn][s_mo], '.', markersize=12,
    #          label='Mod.-EN', color='salmon')
    # ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_ma], principalComponents[:, nn][s_ma], '.', markersize=12,
    #          label='Mod.-LN', color='royalblue')
    # ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_sa], principalComponents[:, nn][s_sa], '.', markersize=12,
    #          label='Str.-LN', color='darkblue')

    ax2.set_ylabel('PC-' + str(nn + 1), fontsize=13)
    ax2.grid('minor')
    ax2.set_xticklabels([])

ax2 = fig.add_subplot(gs2[nn + 1])
ax2.plot(np.unique(DailyDatesMatrix[:,0])[:-1], awt_bmus + 1, 'k.:', linewidth=1.8, markersize=10, color='grey')

# ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_vso], awt_bmus[s_vso] + 1, '.', markersize=17, label='V.Str.-EN',
#          color='darkred')
# ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_so], awt_bmus[s_so] + 1, '.', markersize=17, label='Str.-EN',
#          color='crimson')
# ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_mo], awt_bmus[s_mo] + 1, '.', markersize=17, label='Mod.-EN',
#          color='salmon')
# ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_ma], awt_bmus[s_ma] + 1, '.', markersize=17, label='Mod.-LN',
#          color='royalblue')
# ax2.plot(np.unique(DailyPCs.time.dt.year)[:-1][s_sa], awt_bmus[s_sa] + 1, '.', markersize=17, label='Str.-LN',
#          color='darkblue')

ax2.grid('minor')
ax2.set_xlabel('Time', fontsize=13)
ax2.set_ylabel('AWT', fontsize=13)
ax2.set_ylim([0, 6.5])
#    plt.legend(fontsize=11,loc=2)

#
# pcs_awt = xr.Dataset(
#     {'PCs': (['time', 'n_components'], principalComponents[:, 0:n_components]), 'AWT': (['time'], awt_bmus)},
#     coords={'time': np.unique(DailyPCs.time.dt.year)[:-1], 'n_components': range(n_components)})
# #    pcs_awt.to_netcdf(path=os.path.join(data_dir, 'AWT_dailyProb.nc'))

gs2.tight_layout(fig, rect=[0.22, [], 0.99, []])

gs3 = gridspec.GridSpec(n_components + 1, 1)
for mm in range(n_components):
    ax3 = fig.add_subplot(gs3[mm])
    ax3.pcolormesh(np.rot90((eofs[:, mm].reshape(5, 6))), cmap='RdBu_r', vmin=-0.6, vmax=0.6)
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.set_aspect('equal')

gs3.tight_layout(fig, rect=[0.09, 0.05, 0.2, []])
