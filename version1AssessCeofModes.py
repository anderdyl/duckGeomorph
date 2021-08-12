

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
import sandBarTool.morphLib as mL



import datetime as DT
import numpy as np
import peakutils, os
from matplotlib import pyplot as plt
from scipy import signal


def findSandBarAndTrough1D(xFRF, profile, plotFname=None, **kwargs):
    """ Finds multiple bars on a single profile, will also plot QA/QC plot.  The process begins by finding a trend (see
    keyword args for details) removing it from the data.  Then smooths the profile and trend and finds peaks on that
    detrended line.  It removes portions of the line that are above mean water level (Default=0), that are towards the
    offshore portion of the domain (default = 10%), and sandbars that are within x number of cells (default = 150).  It
    will return the indices of the bars and troughs of interest.

    Assumptions:
        profiles are shoreward to seaward, left to right

    Args:
        xFRF: cross-shore coordinates to describe profile values
        profile: bathymetry values (positive up)
        plotFname: if not None, will generate a QA/QC plot that (default=None)

    Keyword Args
        'waterLevel': finds only bars below this value (Default = 0)
        'trendOrder': removes trend of order (default = 3)
        'deepWaterPercentile': remove sandbars found in the last x percent of the domain (default = 0.9)
        'minSandBarSeparation':  separation in peaks to be found in cells(Default=150 Cells)
        'meanProfile': input a trend to remove data from
        'deepWaterCoordCutoff': a cutoff to say sandbars can't exist beyond this cross-shore coordinate (value must in
                xFRF coordinates)
        'verbose': if on, will print more information to screen (default = False)

    Returns:
        Peak (list): cross-shore location of sandbar crest positions
        Trough (list): cross-shore location of trough positions

    """
    #print('got this far')

    verbose = kwargs.get('verbose', False)
    waterLevel = kwargs.get('waterLevel', 0)
    polyFit = kwargs.get('trendOrder', 3)
    deepwaterPercentile = kwargs.get('deepWaterPercentile', None)
    deepwaterCoordinate = kwargs.get('deepWaterCoordCutoff', 450)
    minSandBarSeparation = kwargs.get('minSandBarSeparation', 5)
    smoothLine = kwargs.get('smoothLengthScale', 3)
    smoothBackground = kwargs.get('smoothBackLengthScale', 5)
    profileTrend = kwargs.get('profileTrend',
                              peakutils.baseline(profile, deg=polyFit))  # find the profiles general trend
    ################################################################
    # start working on data
    assert np.size(profile) == np.size(profileTrend) == np.size(xFRF), 'ProfileTrend must be same size as input profile data, and xFRF'
    #profile_smooth = sb.running_mean(profile, smoothLine)                                                 # first smooth the data to clean out noise
    #xFRF_smooth = sb.running_mean(xFRF, smoothLine)                                                       # smooth the  cross-shore coordinates
    #profileTrend_smooth = sb.running_mean(profileTrend, smoothLine)
    filterDegree = 2   # int(np.ceil(smoothLine/4)*2-1)  # always round to odd number (half of smoothline)
    smoothLine = int(np.ceil(smoothLine/2)*2+1)    # always round to odd number (up from smoothline)
    smoothBackground = int(np.ceil(smoothBackground/2)*2+1)
    profile_smooth = signal.savgol_filter(profile, smoothLine, filterDegree)
    xFRF_smooth = signal.savgol_filter(xFRF, smoothLine, filterDegree)
    profileTrend_smooth = signal.savgol_filter(profileTrend, smoothBackground, filterDegree)
    ### check profile Trend to make sure
    findPeaksOnThisLine = profile_smooth - profileTrend_smooth                                            # remove trend
    findPeaksOnThisLine[profile_smooth >= waterLevel] = np.nan                                            # find only the one's that are below the water level
    if not np.isnan(findPeaksOnThisLine).all():                                                           # if the whole thing is nans' don't do anything
        ###################################
        # 1. find sandbar first cut peaks #
        ###################################
        peakIdx = peakutils.indexes(findPeaksOnThisLine[~np.isnan(findPeaksOnThisLine)], min_dist=minSandBarSeparation)
        # peakIdx = peakutils.indexes(findPeaksOnThisLine[~np.isnan(findPeaksOnThisLine)], min_dist=minSandBarSeparation,
        #                    thres_abs=True, thres=0.25)
        peakIdx = np.argmin(np.isnan(findPeaksOnThisLine)) + peakIdx               # SHIFT back to the original baseline (with nans)
        #if deepwaterPercentile is not None:
        #    peakIdx = peakIdx[peakIdx < len(profile_smooth)*deepwaterPercentile]       # remove any peaks found oceanward of 90% of line
        #else:
        #    peakIdx = peakIdx[xFRF_smooth[peakIdx] < deepwaterCoordinate]

        peakIdx = peakIdx[::-1]  # flip peaks to move offshore to onshore
        ############################################
        # 1a. refine peaks to point of inflection  #
        ############################################
        # make sure now that each peak is actually a local maximum by looking at the slope at each peakIDX and finding
        # the nearest zero slope towards shore
        peakIdxNew, troughIdx = [], []
        shorelineIDX = np.nanargmin(np.abs(profile_smooth))  # identify shoreline on smoothed profile
        for pp, peak in enumerate(peakIdx):
            # find point most shoreward of point with slope greater than zero (add additional point to find point before)
            dElevation = np.diff(profile_smooth[shorelineIDX:peak])                   # take the derivative of the smoothed profile
            if (dElevation > -0.0).any():                                    # are any of the slopes positive
                idxDiff = np.argwhere(dElevation > -0.0).squeeze()   # find all values that have positive slope
                idxMax = np.max(idxDiff) + shorelineIDX              # find max cross-shore location add shoreline location
                peakIdxNew.append(idxMax + 1)                     # add one to find before point of inflection
        peakIdxNew = np.unique(peakIdxNew)
        #######################
        # 2. now find troughs #
        #######################
        for peak in peakIdxNew:
            #if profile_smooth[np.argmin(profile_smooth[:peak]).squeeze()] <  profile_smooth[peak]:      # check that its shallower than the sand bar
            troughIdx.append(np.argmin(profile_smooth[:peak.squeeze()]))
        ########################################################################################
        # 3. check to see if peaks are really peaks, find local maximum between bar and trough #
        ########################################################################################
        # if np.size(peakIdx) == np.size(troughIdx):                     # we have same number of peaks and troughs
        #     for pp, peak in enumerate(peakIdx):
        #         peakIdx[pp] = troughIdx[pp] + np.argmax(profile_smooth[troughIdx[pp]:peak])
        # else:                                                          # we found a peak, but no troughs to the sand bar
        #     # now reiterate the find method, but start from found peak move towards shore
        #     for peak in peakIdx:
        #         findPeaksOnThisLine[range(len(findPeaksOnThisLine)) > peak] = np.nan
        #####################################
        # Last: plot for QA/QC              #
        #####################################
        plt.ioff()                           # turn off plot visible
        if plotFname is not None:            # now plot if interested
            if verbose: print("plotting {}".format(plotFname))
            plt.figure(figsize=(8, 5))
            try:
                plt.suptitle(DT.datetime.strptime(plotFname.split('_')[-1].split('.')[0], '%Y%m%d'))
            except ValueError: # happens when there's not a date in the filename
                plt.suptitle(os.path.basename(plotFname))
            plt.plot(xFRF, profile, 'C1.', label='Raw')
            plt.plot(xFRF_smooth, profile_smooth, 'c.', ms=2, label='smoothed')
            plt.plot(xFRF_smooth, findPeaksOnThisLine, label='Find peaks on this line')
            plt.plot(xFRF_smooth, profileTrend_smooth, label='Trend')
            plt.plot([0, len(profile)], [0,0], 'k--')
            if np.size(peakIdx) > 0:
                plt.plot(xFRF_smooth[peakIdx], profile_smooth[peakIdx], 'r.', ms=5, label='Inital bar Location')
                plt.plot(xFRF_smooth[peakIdx], findPeaksOnThisLine[peakIdx], 'r.', ms=5)
            if np.size(peakIdxNew) >0:
                plt.plot(xFRF_smooth[peakIdxNew], profile_smooth[peakIdxNew], 'rd', ms=7, label='Refined Bar location')
            if np.size(troughIdx) >0:
                plt.plot(xFRF_smooth[troughIdx], profile_smooth[troughIdx], 'bo', ms=6, label='Trough Location')
            plt.legend(loc='upper right')
            plt.ylabel('Elevation NAVD88 [m]')
            plt.xlabel('cross-shore location [m]')
            plt.savefig(plotFname)
            plt.close()
        if np.size(peakIdxNew) > 0 and np.size(troughIdx) > 0:
            return xFRF_smooth[peakIdxNew], xFRF_smooth[troughIdx]
        elif np.size(peakIdxNew) > 0 and np.size(troughIdx) == 0:
            return xFRF_smooth[peakIdxNew], None
        elif np.size(peakIdxNew) == 0 and np.size(troughIdx) > 0:
            return None, xFRF_smooth[troughIdx]
        else:
            return None, None
    else:
        return None, None



dbfile = open('sandbarsSouthernTransect_referencedMHHW_5lineAvg.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()
# alllines = data['alllines']
# xinterp = data['xinterp']
# time = data['time']


dbfile2 = open('ceofsSouthernTransectLatest.pickle', 'rb')
data2 = pickle.load(dbfile2)
dbfile2.close()

timeS = data2['time']
alllinesS = data2['alllines']
xinterpS = data2['xinterp']
SS = data2['S']
RS = data2['Rt']
thetaS = data2['thetaRadians']
theta2S = data2['thetaDegrees']
phitS = data2['phiRadian']
phit2S = data2['phiDegrees']
totalVS = data2['totalV']
lambaS = data2['lamda']
percentVS = data2['percentV']

dbfile3 = open('ceofsNorthernTransect.pickle', 'rb')
data3 = pickle.load(dbfile3)
dbfile3.close()

timeN = data3['time']
alllinesN = data3['alllines']
xinterpN = data3['xinterp']
SN = data3['S']
RN = data3['Rt']
thetaN = data3['thetaRadians']
theta2N = data3['thetaDegrees']
phitN = data3['phiRadian']
phit2N = data3['phiDegrees']
totalVN = data3['totalV']
lambaN = data3['lamda']
percentVN = data3['percentV']


dbfile = open('sandbarsSouthernTransect_referencedMHHW_5lineAvg.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()
alllines = data['alllines']
xinterp = data['xinterp']
time = data['time']



dbfile2 = open('moveInfoOver.pickle', 'rb')
data2 = pickle.load(dbfile2)
dbfile2.close()

colorsC = data2['colors']
flatarrayC = data2['flatarray']
orderAngleC = data2['orderAngle']
orderAngle2C = data2['orderAngle2']
orderC = data2['order']

# for ii in range(len(alllinesS)):
#
#    if ii > 2:
#        fig = plt.figure(figsize=(12, 8))
#        a1 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
#        a1.plot(xinterp, alllines[ii - 2, :], color=[0.8, 0.8, 0.8], linewidth=2,
#                label='{}/{}/{}'.format(time[ii - 2].month, time[ii - 2].day, time[ii - 2].year))
#        a1.plot(xinterp, alllines[ii - 1, :], color=[0.5, 0.5, 0.5], linewidth=3,
#                label='{}/{}/{}'.format(time[ii - 1].month, time[ii - 1].day, time[ii - 1].year))
#        a1.plot(xinterp, alllines[ii, :], 'k', linewidth=3,
#                label='{}/{}/{}'.format(time[ii].month, time[ii].day, time[ii].year))
#        a1.set_title('{}/{}/{}'.format(time[ii].month, time[ii].day, time[ii].year))
#        a1.set_xlim([-20,530])
#        a1.set_ylim([-7.5, 0.5])
#        plt.legend()
#
#        a1.set_xlabel('cross-shore')
#        a1.set_ylabel('depth')
#        plt.savefig('/home/dylananderson/projects/duckGeomorph/alignedProfiles/interpolatedLines{}_date{}_{}'.format(ii,time[ii].month,time[ii].year))
#        plt.close()

dt = time[1:]-time[0:-1]

days = [t.days for t in dt]#[item for sublist in m for item in sublist]
days2 = np.array(days)
plt.figure()
plt.plot(days)

dbfile = open('sandbarsSouthernTransectRandomForestInputData8.pickle','rb')
data = pickle.load(dbfile)
dbfile.close()

currentBmu = data['currentBmu']
prevBmu = data['prevBmu']
daysBetween = data['daysBetween']
avgHs = data['avgHs']
avgTp = data['avgTp']
avgWE = data['avgWE']
avgLWP = data['avgLWP']
avgIr = data['avgIR']
stormCumuLWP = data['stormCumuLWP']
stormCumuWE = data['stormCumuWE']
days = data['days']
sortedBmus = data['sortedBmus']
windSpeed = data['wWindSpeed']
windDirection = data['wWindDirection']
sustWindSpeed = data['wSustWindSpeed']
hs = data['wHs']
tp = data['wTp']
dm = data['wDm']
resWaterLevel = data['wResWaterLevel']
lwp = data['wLWP']
we = data['wWE']
ir = data['wIr']
hoverl = data['wHoverL']
predWaterLevel = data['wPredWaterLevel']
waterLevel = data['wWaterLevel']
hsSplit = data['hsSplit']
tpSplit = data['tpSplit']
dmSplit = data['dmSplit']
lwpSplit = data['lwpSplit']
weSplit = data['weSplit']
hoverlSplit = data['hoverlSplit']
loSplit = data['loSplit']
irSplit = data['irSplit']
waterLevelSplit = data['waterLevelSplit']
predWaterLevelSplit = data['predWaterLevelSplit']
resWaterLevelSplit = data['resWaterLevelSplit']
sustWindSpeedSplit = data['sustWindSpeedSplit']
windSpeedSplit = data['windSpeedSplit']
windDirectionSplit = data['windDirectionSplit']
fvSplit = data['fvSplit']

hsList =data['hsList']
tpList =data['tpList']
dmList =data['dmList']
lwpList =data['lwpList']
weList =data['weList']
fvList =data['fvList']
irList =data['irList']
hlList =data['hlList']


hs90 = np.empty((len(hsList),))
hs50 = np.empty((len(hsList),))

we90 = np.empty((len(hsList),))
lwp90 = np.empty((len(hsList),))
cumWe = np.empty((len(hsList),))
normWe = np.empty((len(hsList),))
fv50 = np.empty((len(hsList),))
hl50 = np.empty((len(hsList),))
for xx in range(len(hsList)):
   hs90[xx] = np.nanpercentile(hsList[xx],90)
   hs50[xx] = np.nanpercentile(hsList[xx],50)

   we90[xx] = np.nanpercentile(weList[xx],90)
   lwp90[xx] = np.nanpercentile(np.abs(lwpList[xx]),90)
   cumWe[xx] = np.nansum(weList[xx])
   normWe[xx] = np.nansum(weList[xx])/days2[xx]
   fv50[xx] = np.nanmean(fvList[xx])
   hl50[xx] = np.nanmean(hlList[xx])


#sortedbmus = data['sortedbmus']

#
# def listToArray(data):
#
#     param05 = [item for sublist in data['param05'] for item in sublist]
#     param05Removed = [ele for ele in param05 if ele != []]
#     param05Array = np.concatenate(param05Removed).ravel()
#
#     param25 = [item for sublist in data['param25'] for item in sublist]
#     param25Removed = [ele for ele in param25 if ele != []]
#     param25Array = np.concatenate(param25Removed).ravel()
#
#     param50 = [item for sublist in data['param50'] for item in sublist]
#     param50Removed = [ele for ele in param50 if ele != []]
#     param50Array = np.concatenate(param50Removed).ravel()
#
#     param75 = [item for sublist in data['param75'] for item in sublist]
#     param75Removed = [ele for ele in param75 if ele != []]
#     param75Array = np.concatenate(param75Removed).ravel()
#
#     param90 = [item for sublist in data['param90'] for item in sublist]
#     param90Removed = [ele for ele in param90 if ele != []]
#     param90Array = np.concatenate(param90Removed).ravel()
#
#     param95 = [item for sublist in data['param95'] for item in sublist]
#     param95Removed = [ele for ele in param95 if ele != []]
#     param95Array = np.concatenate(param95Removed).ravel()
#
#     output = dict()
#     output['param05Array'] = param05Array
#     output['param25Array'] = param25Array
#     output['param50Array'] = param50Array
#     output['param75Array'] = param75Array
#     output['param90Array'] = param90Array
#     output['param95Array'] = param95Array
#     return output
#
#
#
# hsArrays = listToArray(hsSplit)
# hsTran05 = hsArrays['param05Array']
# hsTran25 = hsArrays['param25Array']
# hsTran50 = hsArrays['param50Array']
# hsTran75 = hsArrays['param75Array']
# hsTran90 = hsArrays['param90Array']
# hsTran95 = hsArrays['param95Array']
#
# tpArrays = listToArray(tpSplit)
# tpTran05 = tpArrays['param05Array']
# tpTran25 = tpArrays['param25Array']
# tpTran50 = tpArrays['param50Array']
# tpTran75 = tpArrays['param75Array']
# tpTran90 = tpArrays['param90Array']
# tpTran95 = tpArrays['param95Array']
#
# dmArrays = listToArray(dmSplit)
# dmTran05 = dmArrays['param05Array']
# dmTran25 = dmArrays['param25Array']
# dmTran50 = dmArrays['param50Array']
# dmTran75 = dmArrays['param75Array']
# dmTran90 = dmArrays['param90Array']
# dmTran95 = dmArrays['param95Array']
#
# lwpArrays = listToArray(lwpSplit)
# lwpTran05 = lwpArrays['param05Array']
# lwpTran25 = lwpArrays['param25Array']
# lwpTran50 = lwpArrays['param50Array']
# lwpTran75 = lwpArrays['param75Array']
# lwpTran90 = lwpArrays['param90Array']
# lwpTran95 = lwpArrays['param95Array']
#
# weArrays = listToArray(weSplit)
# weTran05 = weArrays['param05Array']
# weTran25 = weArrays['param25Array']
# weTran50 = weArrays['param50Array']
# weTran75 = weArrays['param75Array']
# weTran90 = weArrays['param90Array']
# weTran95 = weArrays['param95Array']
#
# hoverlArrays = listToArray(hoverlSplit)
# hoverlTran05 = hoverlArrays['param05Array']
# hoverlTran25 = hoverlArrays['param25Array']
# hoverlTran50 = hoverlArrays['param50Array']
# hoverlTran75 = hoverlArrays['param75Array']
# hoverlTran90 = hoverlArrays['param90Array']
# hoverlTran95 = hoverlArrays['param95Array']
#
# irArrays = listToArray(irSplit)
# irTran05 = irArrays['param05Array']
# irTran25 = irArrays['param25Array']
# irTran50 = irArrays['param50Array']
# irTran75 = irArrays['param75Array']
# irTran90 = irArrays['param90Array']
# irTran95 = irArrays['param95Array']
#
# fvArrays = listToArray(fvSplit)
# fvTran05 = fvArrays['param05Array']
# fvTran25 = fvArrays['param25Array']
# fvTran50 = fvArrays['param50Array']
# fvTran75 = fvArrays['param75Array']
# fvTran90 = fvArrays['param90Array']
# fvTran95 = fvArrays['param95Array']
#
# windSpeedArrays = listToArray(windSpeedSplit)
# windSpeedTran05 = windSpeedArrays['param05Array']
# windSpeedTran25 = windSpeedArrays['param25Array']
# windSpeedTran50 = windSpeedArrays['param50Array']
# windSpeedTran75 = windSpeedArrays['param75Array']
# windSpeedTran90 = windSpeedArrays['param90Array']
# windSpeedTran95 = windSpeedArrays['param95Array']
#
# windDirectionArrays = listToArray(windDirectionSplit)
# windDirectionTran05 = windDirectionArrays['param05Array']
# windDirectionTran25 = windDirectionArrays['param25Array']
# windDirectionTran50 = windDirectionArrays['param50Array']
# windDirectionTran75 = windDirectionArrays['param75Array']
# windDirectionTran90 = windDirectionArrays['param90Array']
# windDirectionTran95 = windDirectionArrays['param95Array']
#
# resWaterLevelArrays = listToArray(resWaterLevelSplit)
# resWaterLevelTran05 = resWaterLevelArrays['param05Array']
# resWaterLevelTran25 = resWaterLevelArrays['param25Array']
# resWaterLevelTran50 = resWaterLevelArrays['param50Array']
# resWaterLevelTran75 = resWaterLevelArrays['param75Array']
# resWaterLevelTran90 = resWaterLevelArrays['param90Array']
# resWaterLevelTran95 = resWaterLevelArrays['param95Array']



import matplotlib.colors as col
colors = np.array([[0.91510904, 0.55114749, 0.67037311],
   [0.91696411, 0.55081563, 0.66264366],
   [0.91870995, 0.55055664, 0.65485881],
   [0.92034498, 0.55037149, 0.64702356],
   [0.92186763, 0.55026107, 0.63914306],
   [0.92327636, 0.55022625, 0.63122259],
   [0.9245696 , 0.55026781, 0.62326754],
   [0.92574582, 0.5503865 , 0.6152834 ],
   [0.92680349, 0.55058299, 0.6072758 ],
   [0.92774112, 0.55085789, 0.59925045],
   [0.9285572 , 0.55121174, 0.59121319],
   [0.92925027, 0.551645  , 0.58316992],
   [0.92981889, 0.55215808, 0.57512667],
   [0.93026165, 0.55275127, 0.56708953],
   [0.93057716, 0.5534248 , 0.55906469],
   [0.93076407, 0.55417883, 0.55105838],
   [0.93082107, 0.55501339, 0.54307696],
   [0.93074689, 0.55592845, 0.53512681],
   [0.9305403 , 0.55692387, 0.52721438],
   [0.93020012, 0.55799943, 0.51934621],
   [0.92972523, 0.55915477, 0.51152885],
   [0.92911454, 0.56038948, 0.50376893],
   [0.92836703, 0.56170301, 0.49607312],
   [0.92748175, 0.56309471, 0.48844813],
   [0.9264578 , 0.56456383, 0.48090073],
   [0.92529434, 0.56610951, 0.47343769],
   [0.92399062, 0.56773078, 0.46606586],
   [0.92254595, 0.56942656, 0.45879209],
   [0.92095971, 0.57119566, 0.4516233 ],
   [0.91923137, 0.5730368 , 0.44456642],
   [0.91736048, 0.57494856, 0.4376284 ],
   [0.91534665, 0.57692945, 0.43081625],
   [0.91318962, 0.57897785, 0.42413698],
   [0.91088917, 0.58109205, 0.41759765],
   [0.90844521, 0.58327024, 0.41120533],
   [0.90585771, 0.58551053, 0.40496711],
   [0.90312676, 0.5878109 , 0.3988901 ],
   [0.90025252, 0.59016928, 0.39298143],
   [0.89723527, 0.5925835 , 0.38724821],
   [0.89407538, 0.59505131, 0.38169756],
   [0.89077331, 0.59757038, 0.37633658],
   [0.88732963, 0.60013832, 0.37117234],
   [0.88374501, 0.60275266, 0.36621186],
   [0.88002022, 0.6054109 , 0.36146209],
   [0.87615612, 0.60811044, 0.35692989],
   [0.87215369, 0.61084868, 0.352622  ],
   [0.86801401, 0.61362295, 0.34854502],
   [0.86373824, 0.61643054, 0.34470535],
   [0.85932766, 0.61926872, 0.3411092 ],
   [0.85478365, 0.62213474, 0.3377625 ],
   [0.85010767, 0.6250258 , 0.33467091],
   [0.84530131, 0.62793914, 0.3318397 ],
   [0.84036623, 0.63087193, 0.32927381],
   [0.8353042 , 0.63382139, 0.32697771],
   [0.83011708, 0.63678472, 0.32495541],
   [0.82480682, 0.63975913, 0.32321038],
   [0.81937548, 0.64274185, 0.32174556],
   [0.81382519, 0.64573011, 0.32056327],
   [0.80815818, 0.6487212 , 0.31966522],
   [0.80237677, 0.65171241, 0.31905244],
   [0.79648336, 0.65470106, 0.31872531],
   [0.79048044, 0.65768455, 0.31868352],
   [0.78437059, 0.66066026, 0.31892606],
   [0.77815645, 0.66362567, 0.31945124],
   [0.77184076, 0.66657827, 0.32025669],
   [0.76542634, 0.66951562, 0.3213394 ],
   [0.75891609, 0.67243534, 0.32269572],
   [0.75231298, 0.67533509, 0.32432138],
   [0.74562004, 0.6782126 , 0.32621159],
   [0.73884042, 0.68106567, 0.32836102],
   [0.73197731, 0.68389214, 0.33076388],
   [0.72503398, 0.68668995, 0.33341395],
   [0.7180138 , 0.68945708, 0.33630465],
   [0.71092018, 0.69219158, 0.33942908],
   [0.70375663, 0.69489159, 0.34278007],
   [0.69652673, 0.69755529, 0.34635023],
   [0.68923414, 0.70018097, 0.35013201],
   [0.6818826 , 0.70276695, 0.35411772],
   [0.67447591, 0.70531165, 0.3582996 ],
   [0.667018  , 0.70781354, 0.36266984],
   [0.65951284, 0.71027119, 0.36722061],
   [0.65196451, 0.71268322, 0.37194411],
   [0.64437719, 0.71504832, 0.37683259],
   [0.63675512, 0.71736525, 0.38187838],
   [0.62910269, 0.71963286, 0.38707389],
   [0.62142435, 0.72185004, 0.39241165],
   [0.61372469, 0.72401576, 0.39788432],
   [0.60600841, 0.72612907, 0.40348469],
   [0.59828032, 0.72818906, 0.40920573],
   [0.59054536, 0.73019489, 0.41504052],
   [0.58280863, 0.73214581, 0.42098233],
   [0.57507535, 0.7340411 , 0.42702461],
   [0.5673509 , 0.7358801 , 0.43316094],
   [0.55964082, 0.73766224, 0.43938511],
   [0.55195081, 0.73938697, 0.44569104],
   [0.54428677, 0.74105381, 0.45207286],
   [0.53665478, 0.74266235, 0.45852483],
   [0.52906111, 0.74421221, 0.4650414 ],
   [0.52151225, 0.74570306, 0.47161718],
   [0.5140149 , 0.74713464, 0.47824691],
   [0.506576  , 0.74850672, 0.48492552],
   [0.49920271, 0.74981912, 0.49164808],
   [0.49190247, 0.75107171, 0.4984098 ],
   [0.48468293, 0.75226438, 0.50520604],
   [0.47755205, 0.7533971 , 0.51203229],
   [0.47051802, 0.75446984, 0.5188842 ],
   [0.46358932, 0.75548263, 0.52575752],
   [0.45677469, 0.75643553, 0.53264815],
   [0.45008317, 0.75732863, 0.5395521 ],
   [0.44352403, 0.75816207, 0.54646551],
   [0.43710682, 0.758936  , 0.55338462],
   [0.43084133, 0.7596506 , 0.56030581],
   [0.42473758, 0.76030611, 0.56722555],
   [0.41880579, 0.76090275, 0.5741404 ],
   [0.41305637, 0.76144081, 0.58104704],
   [0.40749984, 0.76192057, 0.58794226],
   [0.40214685, 0.76234235, 0.59482292],
   [0.39700806, 0.7627065 , 0.60168598],
   [0.39209414, 0.76301337, 0.6085285 ],
   [0.38741566, 0.76326334, 0.6153476 ],
   [0.38298304, 0.76345681, 0.62214052],
   [0.37880647, 0.7635942 , 0.62890454],
   [0.37489579, 0.76367593, 0.63563704],
   [0.37126045, 0.76370246, 0.64233547],
   [0.36790936, 0.76367425, 0.64899736],
   [0.36485083, 0.76359176, 0.6556203 ],
   [0.36209245, 0.76345549, 0.66220193],
   [0.359641  , 0.76326594, 0.66873999],
   [0.35750235, 0.76302361, 0.67523226],
   [0.35568141, 0.76272903, 0.68167659],
   [0.35418202, 0.76238272, 0.68807086],
   [0.3530069 , 0.76198523, 0.69441305],
   [0.35215761, 0.7615371 , 0.70070115],
   [0.35163454, 0.76103888, 0.70693324],
   [0.35143685, 0.76049114, 0.71310742],
   [0.35156253, 0.75989444, 0.71922184],
   [0.35200839, 0.75924936, 0.72527472],
   [0.3527701 , 0.75855647, 0.73126429],
   [0.3538423 , 0.75781637, 0.73718884],
   [0.3552186 , 0.75702964, 0.7430467 ],
   [0.35689171, 0.75619688, 0.74883624],
   [0.35885353, 0.75531868, 0.75455584],
   [0.36109522, 0.75439565, 0.76020396],
   [0.36360734, 0.75342839, 0.76577905],
   [0.36637995, 0.75241752, 0.77127961],
   [0.3694027 , 0.75136364, 0.77670417],
   [0.37266493, 0.75026738, 0.7820513 ],
   [0.37615579, 0.74912934, 0.78731957],
   [0.37986429, 0.74795017, 0.79250759],
   [0.38377944, 0.74673047, 0.797614  ],
   [0.38789026, 0.74547088, 0.80263746],
   [0.3921859 , 0.74417203, 0.80757663],
   [0.39665568, 0.74283455, 0.81243022],
   [0.40128912, 0.74145908, 0.81719695],
   [0.406076  , 0.74004626, 0.82187554],
   [0.41100641, 0.73859673, 0.82646476],
   [0.41607073, 0.73711114, 0.83096336],
   [0.4212597 , 0.73559013, 0.83537014],
   [0.42656439, 0.73403435, 0.83968388],
   [0.43197625, 0.73244447, 0.8439034 ],
   [0.43748708, 0.73082114, 0.84802751],
   [0.44308905, 0.72916502, 0.85205505],
   [0.44877471, 0.72747678, 0.85598486],
   [0.45453694, 0.72575709, 0.85981579],
   [0.46036897, 0.72400662, 0.8635467 ],
   [0.4662644 , 0.72222606, 0.86717646],
   [0.47221713, 0.72041608, 0.87070395],
   [0.47822138, 0.71857738, 0.87412804],
   [0.4842717 , 0.71671065, 0.87744763],
   [0.4903629 , 0.71481659, 0.88066162],
   [0.49649009, 0.71289591, 0.8837689 ],
   [0.50264864, 0.71094931, 0.88676838],
   [0.50883417, 0.70897752, 0.88965898],
   [0.51504253, 0.70698127, 0.89243961],
   [0.52126981, 0.70496128, 0.8951092 ],
   [0.52751231, 0.70291829, 0.89766666],
   [0.53376652, 0.70085306, 0.90011093],
   [0.54002912, 0.69876633, 0.90244095],
   [0.54629699, 0.69665888, 0.90465565],
   [0.55256715, 0.69453147, 0.90675397],
   [0.55883679, 0.69238489, 0.90873487],
   [0.56510323, 0.69021993, 0.9105973 ],
   [0.57136396, 0.68803739, 0.91234022],
   [0.57761655, 0.68583808, 0.91396258],
   [0.58385872, 0.68362282, 0.91546336],
   [0.59008831, 0.68139246, 0.91684154],
   [0.59630323, 0.67914782, 0.9180961 ],
   [0.60250152, 0.67688977, 0.91922603],
   [0.60868128, 0.67461918, 0.92023033],
   [0.61484071, 0.67233692, 0.921108  ],
   [0.62097809, 0.67004388, 0.92185807],
   [0.62709176, 0.66774097, 0.92247957],
   [0.63318012, 0.66542911, 0.92297153],
   [0.63924166, 0.66310923, 0.92333301],
   [0.64527488, 0.66078227, 0.92356308],
   [0.65127837, 0.65844919, 0.92366082],
   [0.65725076, 0.65611096, 0.92362532],
   [0.66319071, 0.65376857, 0.92345572],
   [0.66909691, 0.65142302, 0.92315115],
   [0.67496813, 0.64907533, 0.92271076],
   [0.68080311, 0.64672651, 0.92213374],
   [0.68660068, 0.64437763, 0.92141929],
   [0.69235965, 0.64202973, 0.92056665],
   [0.69807888, 0.6396839 , 0.91957507],
   [0.70375724, 0.63734122, 0.91844386],
   [0.70939361, 0.63500279, 0.91717232],
   [0.7149869 , 0.63266974, 0.91575983],
   [0.72053602, 0.63034321, 0.91420578],
   [0.72603991, 0.62802433, 0.9125096 ],
   [0.7314975 , 0.62571429, 0.91067077],
   [0.73690773, 0.62341425, 0.9086888 ],
   [0.74226956, 0.62112542, 0.90656328],
   [0.74758193, 0.61884899, 0.90429382],
   [0.75284381, 0.6165862 , 0.90188009],
   [0.75805413, 0.61433829, 0.89932181],
   [0.76321187, 0.6121065 , 0.89661877],
   [0.76831596, 0.6098921 , 0.89377082],
   [0.77336536, 0.60769637, 0.89077786],
   [0.77835901, 0.6055206 , 0.88763988],
   [0.78329583, 0.6033661 , 0.88435693],
   [0.78817477, 0.60123418, 0.88092913],
   [0.79299473, 0.59912616, 0.87735668],
   [0.79775462, 0.59704339, 0.87363986],
   [0.80245335, 0.59498722, 0.86977904],
   [0.8070898 , 0.592959  , 0.86577468],
   [0.81166284, 0.5909601 , 0.86162732],
   [0.81617134, 0.5889919 , 0.8573376 ],
   [0.82061414, 0.58705579, 0.85290625],
   [0.82499007, 0.58515315, 0.84833413],
   [0.82929796, 0.58328538, 0.84362217],
   [0.83353661, 0.58145389, 0.83877142],
   [0.8377048 , 0.57966009, 0.83378306],
   [0.8418013 , 0.57790538, 0.82865836],
   [0.84582486, 0.57619119, 0.82339871],
   [0.84977422, 0.57451892, 0.81800565],
   [0.85364809, 0.57289   , 0.8124808 ],
   [0.85744519, 0.57130585, 0.80682595],
   [0.86116418, 0.56976788, 0.80104298],
   [0.86480373, 0.56827749, 0.79513394],
   [0.86836249, 0.56683612, 0.789101  ],
   [0.87183909, 0.56544515, 0.78294645],
   [0.87523214, 0.56410599, 0.77667274],
   [0.87854024, 0.56282002, 0.77028247],
   [0.88176195, 0.56158863, 0.76377835],
   [0.88489584, 0.56041319, 0.75716326],
   [0.88794045, 0.55929505, 0.75044023],
   [0.89089432, 0.55823556, 0.74361241],
   [0.89375596, 0.55723605, 0.73668312],
   [0.89652387, 0.55629781, 0.72965583],
   [0.89919653, 0.55542215, 0.72253414],
   [0.90177242, 0.55461033, 0.71532181],
   [0.90425   , 0.55386358, 0.70802274],
   [0.90662774, 0.55318313, 0.70064098],
   [0.90890408, 0.55257016, 0.69318073],
   [0.91107745, 0.55202582, 0.68564633],
   [0.91314629, 0.55155124, 0.67804225]])

cmap = col.ListedColormap(colors)
n = 50
# plt.figure()
# plt.plot(timeS[0:n],(phitS[0:n,1]))
# plt.plot(timeN[0:n],(phitN[0:n,1]))
# plt.plot(timeS[0:n],np.unwrap(phitS[0:n,1],discont=3*np.pi/2),'-.',color='r')
# plt.plot(timeN[0:n],np.unwrap(phitN[0:n,1],discont=3*np.pi/2),'--',color='b')
#
#
# plt.figure()
# plt.plot(timeS,np.unwrap(phitS[:,1],discont=2.*np.pi/2),'-.',color='r')
# plt.plot(timeN,np.unwrap(phitN[:,1],discont=2.*np.pi/2),'--',color='b')
#
# plt.plot(timeS,np.unwrap(phitS[:,2],discont=2.*np.pi/2),'-.',color='r')
# plt.plot(timeN,np.unwrap(phitN[:,2],discont=2.*np.pi/2),'--',color='b')
#





def P2R(radii, angles):
    return radii * np.exp(1j*angles)

def R2P(x):
    return np.abs(x), np.angle(x)



timeind = np.arange(0,len(time))
#timeind = np.arange(3, 50)
# RtSubset = Rt[timeind, :]
# phitSubset = phit[timeind, :]
# phit2Subset = phit2[timeind, :]
# timeSubset = time[timeind]
alllinesSubset = alllines[timeind, :]

eofPredog = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred2og = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred3og = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred4og = np.nan * np.ones((np.shape(alllinesSubset)))

for timestep in range(len(timeind)):
    mode = 0
    eofPredog[timestep, :] = RS[timestep, mode]* SS[:, mode] * np.cos(phitS[timestep, mode] - thetaS[:, mode])
    mode = 1
    eofPred2og[timestep, :] = RS[timestep, mode]* SS[:, mode] * np.cos(phitS[timestep, mode] - thetaS[:, mode])

    # mode = 2
    # eofPred3[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
    #
    # mode = 3
    # eofPred4[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])




t1 = 0
t2 = -1
fig = plt.figure(figsize=(14,10))
#plt.set_cmap('RdBu')#bwr')
plt.set_cmap('RdBu_r')

tg, xg = np.meshgrid(time, xinterp)
# ax1 = plt.subplot2grid((4,4),(0,0),rowspan=4,colspan=1)
# plt0 = ax1.pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
# fig.colorbar(plt0, ax=ax1, orientation='horizontal')
# ax1.set_ylim([time[t1], time[t2]])
# ax1.set_title('Surveys (dev.)')

ax2 = plt.subplot2grid((7,4),(0,0),rowspan=7,colspan=1)
plt1 = ax2.pcolor(xg,tg,eofPredog.T, vmin=-.75, vmax=0.75)
ax2.set_ylim([time[t1], time[t2]])
#fig.colorbar(plt1, ax=ax2, orientation='horizontal')
ax2.set_title('CEOF1 {:.2f}'.format(percentVS[0]))
ax2.get_yaxis().set_ticks([])

ax3 = plt.subplot2grid((7,4),(0,1),rowspan=7,colspan=1)
plt2 = ax3.pcolor(xg,tg,eofPred2og.T, vmin=-.85, vmax=.85)
ax3.set_ylim([time[t1], time[t2]])
#fig.colorbar(plt2, ax=ax3, orientation='horizontal')
ax3.set_title('CEOF2 {:.2f}'.format(percentVS[1]))
ax3.get_yaxis().set_ticks([])

plt.show()



fig, ax = plt.subplots(2,2)
mode = 0
ax[0,0].scatter(xinterpS, SS[:,mode],16,theta2S[:,mode],cmap=cmap)
ax[0,0].set_title('Mode 1 (40% of variance)',fontsize=24)
ax[0,0].set_ylabel('Spatial Magnitude (m)',fontsize=18)
ax[0,0].set_xlabel('Cross-shore (m)',fontsize=14)
ax[1,0].scatter(xinterpS, theta2S[:,mode],12,theta2S[:,mode],cmap=cmap)
ax[1,0].set_ylabel('Spatial Phase (deg)',fontsize=18)
ax[1,0].set_xlabel('Cross-shore (m)',fontsize=14)

mode = 1
ax[0,1].scatter(xinterpS, SS[:,mode],16,theta2S[:,mode],cmap=cmap)
ax[0,1].set_title('Mode 2 (25% of variance)',fontsize=24)
ax[0,1].set_ylabel('Spatial Magnitude (m)',fontsize=18)
ax[0,1].set_xlabel('Cross-shore (m)',fontsize=14)
ax[1,1].scatter(xinterpS, theta2S[:,mode],12,theta2S[:,mode],cmap=cmap)
ax[1,1].set_ylabel('Spatial Phase (deg)',fontsize=18)
ax[1,1].set_xlabel('Cross-shore (m)',fontsize=14)

plt.show()

# ax[0,1].scatter(timeS,RS[:,mode],12,phit2S[:,mode])
# ax[0,1].set_ylabel('Temporal Magnitude (m)')
# ax[0,1].set_xlabel('Time')
# ax[1,1].scatter(timeS,phit2S[:,mode],12,phit2S[:,mode])
# ax[1,1].set_ylabel('Temporal Phase (deg)')
# ax[1,1].set_xlabel('Time')

#
# from pandas import rolling_median
#
# threshold = 3
# df['pandas'] = rolling_median(df['u'], window=3, center=True).fillna(method='bfill').fillna(method='ffill')
#
# difference = np.abs(df['u'] - df['pandas'])
# outlier_idx = difference > threshold
#
# fig, ax = plt.subplots(figsize=figsize)
# df['u'].plot()
# df['u'][outlier_idx].plot(**kw)
# _ = ax.set_ylim(-50, 50)




fig10 = plt.figure(figsize=(10,10))
thetaMode1 = -90*np.pi/180
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
thetaMode2 = np.arange(0,180,10)
ax1 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1.plot(xinterp,combined,color=colors[ii,:])

thetaMode1 = 90*np.pi/180
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
thetaMode2 = np.arange(0,180,10)
ax2 = plt.subplot2grid((1,2),(0,1),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax2.plot(xinterp,combined,color=colors[ii,:])

# ax1.plot(xinterp,np.mean(alllinesS,axis=0)+eofPred)



fig10 = plt.figure(figsize=(10,10))
thetaMode2 = -90*np.pi/180
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2 - thetaS[:,1])

thetaMode1 = np.arange(0,180,10)
ax1 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])

   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1.plot(xinterp,combined,color=colors[ii,:])








fig10 = plt.figure(figsize=(15,10))

thetaMode2 = 150
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])

thetaMode1 = np.arange(90,280,20)
ax1c = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])

   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))

ax1c.legend()
ax1c.set_title(r'Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))

thetaMode2 = 0
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])

thetaMode1 = np.arange(-90,100,20)
ax1 = plt.subplot2grid((2,3),(0,1),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])

   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))

ax1.legend()
ax1.set_title(r'Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))


thetaMode2 = 180
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])

thetaMode1 = np.arange(-90,100,20)
ax1d = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])

   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1d.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))

ax1d.legend()
ax1d.set_title(r'Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))


thetaMode1 = -90
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])
thetaMode2 = np.arange(0,190,20)
ax1b = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1b.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_2$ = {}'.format(thetaMode2[ii]))
ax1b.legend()
ax1b.set_title(r'Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))

thetaMode1 = 75
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])
thetaMode2 = np.arange(-100,50,20)
ax2 = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax2.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_2$ = {}'.format(thetaMode2[ii]))
ax2.legend()
ax2.set_title(r'Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))

thetaMode1 = -10
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])
thetaMode2 = np.arange(-70,50,10)
ax2c = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax2c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_2$ = {}'.format(thetaMode2[ii]))
ax2c.legend()
ax2c.set_title(r'Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))

plt.show()



fig10 = plt.figure(figsize=(15,10))


thetaMode2 =5
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])

thetaMode1 = np.arange(-110,30,10)
ax1c = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])

   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))

#ax1c.legend()
ax1c.set_title(r'2a. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))



thetaMode1 = -90
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred1 = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])

thetaMode2 = np.arange(0,180,10)
ax2a = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:,1])

   combined = np.mean(alllinesS,axis=0)+eofPred1+eofPred2
   ax2a.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode2[ii]))

#ax2a.legend()
ax2a.set_title(r'2b. Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))



thetaMode2 =15
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])

thetaMode1 = np.arange(10,120,5)
ax1c = plt.subplot2grid((2,3),(0,1),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])

   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))

#ax1c.legend()
ax1c.set_title(r'3a. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))


# thetaMode2 =200
# # eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# # eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
# eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])
#
# thetaMode1 = np.arange(-40,75,5)
# ax1c = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
# import matplotlib.cm as cm
# colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
# for ii in range(len(thetaMode1)):
#    # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
#    eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
#
#    combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
#    ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))
#
# #ax1c.legend()
# ax1c.set_title(r'3b. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))

thetaMode1 = 60
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred1 = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])

thetaMode2 = np.arange(40,140,10)
ax2a = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:,1])

   combined = np.mean(alllinesS,axis=0)+eofPred1+eofPred2
   ax2a.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode2[ii]))

#ax2a.legend()
ax2a.set_title(r'3c. Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))


thetaMode1 = -90
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred1 = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])

thetaMode2 = np.arange(145,290,10)
ax2a = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:,1])

   combined = np.mean(alllinesS,axis=0)+eofPred1+eofPred2
   ax2a.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode2[ii]))

#ax2a.legend()
ax2a.set_title(r'3c. Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))



thetaMode2 =25
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
# eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1 - thetaS[:,0])
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])

thetaMode1 = np.arange(-0,135,5)
ax1c = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   # = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])

   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))

#ax1c.legend()
ax1c.set_title(r'3b. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))

plt.show()












fig10 = plt.figure(figsize=(15,10))


thetaMode2 =5
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])

thetaMode1 = np.arange(-110,30,10)
ax1c = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))
#ax1c.legend()
ax1c.set_title(r'2a. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))

thetaMode2 =25
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])
thetaMode1 = np.arange(10,120,5)
ax1c = plt.subplot2grid((2,3),(0,1),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))
#ax1c.legend()
ax1c.set_title(r'3a. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))


thetaMode2 = 150
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])
thetaMode1 = np.arange(80,230,10)
ax1d = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1d.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))
# ax1d.legend()
ax1d.set_title(r'Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))


thetaMode1 = 60
eofPred1 = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])

thetaMode2 = np.arange(40,140,10)
ax2a = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:,1])
   combined = np.mean(alllinesS,axis=0)+eofPred1+eofPred2
   ax2a.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode2[ii]))
#ax2a.legend()
ax2a.set_title(r'3c. Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))


thetaMode1 = 75
# eofPred = RtSubset[timestep, mode] * RS[:, 0] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])
thetaMode2 = np.arange(-100,50,20)
ax2 = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax2.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_2$ = {}'.format(thetaMode2[ii]))
# ax2.legend()
ax2.set_title(r'Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))

thetaMode1 = -90
eofPred = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])
thetaMode2 = np.arange(20,190,10)
ax1b = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:, 1])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1b.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_2$ = {}'.format(thetaMode2[ii]))
# ax1b.legend()
ax1b.set_title(r'Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))
plt.tight_layout()
plt.show()





thetaMode1 = -90
eofPred1 = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])
thetaMode2 = np.arange(0,180,10)
ax2a = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:,1])
   combined = np.mean(alllinesS,axis=0)+eofPred1+eofPred2
   ax2a.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode2[ii]))
#ax2a.legend()
ax2a.set_title(r'2b. Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))



thetaMode1 = 60
eofPred1 = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])

thetaMode2 = np.arange(40,140,10)
ax2a = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:,1])
   combined = np.mean(alllinesS,axis=0)+eofPred1+eofPred2
   ax2a.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode2[ii]))
#ax2a.legend()
ax2a.set_title(r'3c. Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))

thetaMode1 = -90
eofPred1 = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode1*np.pi/180 - thetaS[:,0])
thetaMode2 = np.arange(145,290,10)
ax2a = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode2)))
for ii in range(len(thetaMode2)):
   eofPred2 = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode2[ii]*np.pi/180 - thetaS[:,1])
   combined = np.mean(alllinesS,axis=0)+eofPred1+eofPred2
   ax2a.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode2[ii]))
#ax2a.legend()
ax2a.set_title(r'3c. Stationary $\phi_1$ = {} + Propagating Onshore Signal'.format(thetaMode1))


thetaMode2 =25
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])
thetaMode1 = np.arange(-0,135,5)
ax1c = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))
#ax1c.legend()
ax1c.set_title(r'3b. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))

plt.show()




mean1 = np.nanmean(RS[:, 0])
mean2 = np.nanmean(RS[:, 1])
# angleMode1a = [-40,-20,  0, 20, 40, 60, 80,100,120]
# angleMode2a = [-40,-40,-40,-40,-40,-40,-40,-40,-40]
angleMode1a = [-40,-20,  0, 20]
angleMode2a = [-40,-40,-40,-40]
amp1a = [mean1,mean1,mean1,mean1]
amp2a = np.multiply(0.75,[mean2,mean2,mean2,mean2])
# angleMode1b = [120,120,120,120,120,120,120]
# angleMode2b = [ -20, 0, 20, 40, 60, 80,100]
angleMode1b = [20,20,20,20]#,20,20,20]
angleMode2b = [ -20, 0, 20, 40]#, 60, 80,100]
amp1b = np.multiply(0.75,[mean1,mean1,mean1,mean1,mean1,mean1,mean1])
amp2b = [mean2,mean2,mean2,mean2,mean2,mean2,mean2]
angleMode1c = np.hstack((angleMode1a,angleMode1b))
angleMode2c = np.hstack((angleMode2a,angleMode2b))
amp1c = np.hstack((amp1a,amp1b))
amp2c = np.hstack((amp2a,amp2b))
angleMode1d = [20,40,60,80,100]#,120,140,160,180,200]
angleMode2d = [40,40,40,40, 40]#,40,40,40,40,40]
amp1d = [mean1,mean1,mean1,mean1,1]#,1.2,1.3,1.4,1.5,1.6]
# amp2d = np.multiply(0.3,[mean2,mean2,mean2,mean2,mean2,mean2,mean2,mean2,mean2,mean2])
amp2d = [0.8,0.7,0.6,0.3,0.2]#,0.1,0.2,0.2,0.2,0.2]
angleMode1e = np.hstack((angleMode1c,angleMode1d))
angleMode2e = np.hstack((angleMode2c,angleMode2d))
amp1e = np.hstack((amp1c,amp1d))
amp2e = np.hstack((amp2c,amp2d))

angleMode1f = [100,100,100,100,100,100,100]
angleMode2f = [ 60, 80,100,120,140,160,180]
amp1f = [0.7,0.7,0.7,0.7,0.7,0.7,0.7]#,1.2,1.3,1.4,1.5,1.6]
amp2f = [mean2,1.05,1.2,1.1,1.1,1.1,1.1]#,0.1,0.2,0.2,0.2,0.2]
angleMode1g = np.hstack((angleMode1e,angleMode1f))
angleMode2g = np.hstack((angleMode2e,angleMode2f))
amp1g = np.hstack((amp1e,amp1f))
amp2g = np.hstack((amp2e,amp2f))

angleMode1h = [120,140,160,180,200,220,240,260,280]#,120,140,160,180,200]
angleMode2h = [ 180, 180,180,180,180,180,180,180,180]#,40,40,40,40,40]
amp1h = [0.6,0.9,1.2,1.1,1.1,1.1,1.1,1.1,.9]#,1.2,1.3,1.4,1.5,1.6]
amp2h = [1,.8,.6,.5,.5,.5,.6,.7,.8]#,0.1,0.2,0.2,0.2,0.2]
angleMode1i = np.hstack((angleMode1g,angleMode1h))
angleMode2i = np.hstack((angleMode2g,angleMode2h))
amp1i = np.hstack((amp1g,amp1h))
amp2i = np.hstack((amp2g,amp2h))

angleMode1j = [280,280,280,280]#,280,280,280]#,120,140,160,180,200]
angleMode2j = [200,220,240,260]#,280,300,320]#,40,40,40,40,40]
amp1j = [0.9,0.9,0.9,0.9]#,0.9,0.9,0.9]#,1.2,1.3,1.4,1.5,1.6]
amp2j = [0.9,0.9,0.9,0.9]#,0.9,0.9,0.9]#,0.1,0.2,0.2,0.2,0.2]
angleMode1k = np.hstack((angleMode1i,angleMode1j))
angleMode2k = np.hstack((angleMode2i,angleMode2j))
amp1k = np.hstack((amp1i,amp1j))
amp2k = np.hstack((amp2i,amp2j))

angleMode1l = [300,320,340,360]#,380,400]#,280,280,280]#,120,140,160,180,200]
angleMode2l = [260,260,260,260]#,260,400]#,280,300,320]#,40,40,40,40,40]
amp1l = [0.9,0.7,0.5,0.9]#,0.9,0.9]#,1.2,1.3,1.4,1.5,1.6]
amp2l = [0.9,0.99,1.1,0.9]#,0.9,0.9]#,0.1,0.2,0.2,0.2,0.2]
angleMode1m = np.hstack((angleMode1k,angleMode1l))
angleMode2m = np.hstack((angleMode2k,angleMode2l))
amp1m = np.hstack((amp1k,amp1l))
amp2m = np.hstack((amp2k,amp2l))

angleMode1n = [360,360,360,360]#,380,400]#,280,280,280]#,120,140,160,180,200]
angleMode2n = [280,300,320,340]#,260,400]#,280,300,320]#,40,40,40,40,40]
amp1n = [0.9,1.1,0.99,0.9]#,0.9,0.9]#,1.2,1.3,1.4,1.5,1.6]
amp2n = [0.9,0.99,1.1,0.9]#,0.9,0.9]#,0.1,0.2,0.2,0.2,0.2]
angleMode1o = np.hstack((angleMode1m,angleMode1n))
angleMode2o = np.hstack((angleMode2m,angleMode2n))
amp1o = np.hstack((amp1m,amp1n))
amp2o = np.hstack((amp2m,amp2n))

angleMode1p = [380,400,420,440,460,480]#,380,400]#,280,280,280]#,120,140,160,180,200]
angleMode2p = [340,340,340,340,340,340]#,260,400]#,280,300,320]#,40,40,40,40,40]
amp1p = [0.9,1.1,0.99,0.9,0.9,0.8,0.7]#,0.9,0.9]#,1.2,1.3,1.4,1.5,1.6]
amp2p = [0.9,0.99,1.1,0.9,0.8,0.8,0.8]#,0.9,0.9]#,0.1,0.2,0.2,0.2,0.2]
angleMode1 = np.hstack((angleMode1o,angleMode1p))
angleMode2 = np.hstack((angleMode2o,angleMode2p))
amp1 = np.hstack((amp1o,amp1p))
amp2 = np.hstack((amp2o,amp2p))

modeColors = cm.rainbow(np.linspace(0, 1, len(angleMode1)))

for ii in range(len(angleMode1)):
    plt.figure(figsize=(11, 5))
    if ii == 0:
        axLeft = plt.subplot2grid((2, 4), (0, 0), rowspan=2, colspan=2)
        axLeft.plot(angleMode1[0], angleMode2[0], '-')
        axLeft.scatter(angleMode1[0],angleMode2[0],15,color=modeColors[0,:],zorder=3)
        axLeft.set_xlim([-100, 500])
        axLeft.set_ylim([-100, 500])
        axLeft.set_xlabel('Offshore Migration (phase) $\Longrightarrow$')
        axLeft.set_ylabel('Onshore Migration (phase) $\Longrightarrow$')
        axRight = plt.subplot2grid((2, 4), (0, 2), rowspan=2, colspan=2)
        ceofMode1 = SS[:, 0] * amp1[ii] * np.cos(angleMode1[ii] * np.pi / 180 - thetaS[:, 0])
        ceofMode2 = SS[:, 1] * amp2[ii] * np.cos(angleMode2[ii] * np.pi / 180 - thetaS[:, 1])
        combined = np.mean(alllinesS, axis=0) + ceofMode1 + ceofMode2
        axRight.plot(xinterp,combined,color=modeColors[ii,:],label=r'$\phi_1$ = {}'.format(angleMode1[ii]))
    else:
        axLeft = plt.subplot2grid((2, 4), (0, 0), rowspan=2, colspan=2)
        axLeft.plot(angleMode1[0:ii+1], angleMode2[0:ii+1], '-')
        axLeft.scatter(angleMode1[0:ii+1],angleMode2[0:ii+1],15,color=modeColors[0:ii+1,:],zorder=3)
        axLeft.set_xlim([-100, 500])
        axLeft.set_ylim([-100, 500])
        axLeft.set_xlabel('Offshore Migration (phase) $\Longrightarrow$')
        axLeft.set_ylabel('Onshore Migration (phase) $\Longrightarrow$')
        axRight = plt.subplot2grid((2, 4), (0, 2), rowspan=2, colspan=2)
        for ff in range(ii+1):
            if ff < ii-2:
                ceofMode1 = SS[:, 0] * amp1[ff] * np.cos(angleMode1[ff] * np.pi / 180 - thetaS[:, 0])
                ceofMode2 = SS[:, 1] * amp2[ff] * np.cos(angleMode2[ff] * np.pi / 180 - thetaS[:, 1])
                combined = np.mean(alllinesS, axis=0) + ceofMode1 + ceofMode2
                axRight.plot(xinterp,combined,color=[0.3,0.3,0.3],alpha=0.4,label=r'$\phi_1$ = {}'.format(angleMode1[ff]))
            elif ff == ii-2:
                ceofMode1 = SS[:, 0] * amp1[ff] * np.cos(angleMode1[ff] * np.pi / 180 - thetaS[:, 0])
                ceofMode2 = SS[:, 1] * amp2[ff] * np.cos(angleMode2[ff] * np.pi / 180 - thetaS[:, 1])
                combined = np.mean(alllinesS, axis=0) + ceofMode1 + ceofMode2
                axRight.plot(xinterp,combined,color=modeColors[ff,:],alpha=0.4,label=r'$\phi_1$ = {}'.format(angleMode1[ff]))
            elif ff == ii-1:
                ceofMode1 = SS[:, 0] * amp1[ff] * np.cos(angleMode1[ff] * np.pi / 180 - thetaS[:, 0])
                ceofMode2 = SS[:, 1] * amp2[ff] * np.cos(angleMode2[ff] * np.pi / 180 - thetaS[:, 1])
                combined = np.mean(alllinesS, axis=0) + ceofMode1 + ceofMode2
                axRight.plot(xinterp,combined,color=modeColors[ff,:],alpha=0.7,label=r'$\phi_1$ = {}'.format(angleMode1[ff]))
            else:
                ceofMode1 = SS[:, 0] * amp1[ff] * np.cos(angleMode1[ff] * np.pi / 180 - thetaS[:, 0])
                ceofMode2 = SS[:, 1] * amp2[ff] * np.cos(angleMode2[ff] * np.pi / 180 - thetaS[:, 1])
                combined = np.mean(alllinesS, axis=0) + ceofMode1 + ceofMode2
                axRight.plot(xinterp,combined,color=modeColors[ff,:],label=r'$\phi_1$ = {}'.format(angleMode1[ff]))

    axRight.set_xlabel('cross-shore (m)')
    axRight.set_ylabel('depth (m)')
    plt.tight_layout()
    if ii < 10:
        plt.savefig('/home/dylananderson/projects/duckGeomorph/migrationPlots/figure0{}.png'.format(ii),dpi=300)
    else:
        plt.savefig('/home/dylananderson/projects/duckGeomorph/migrationPlots/figure{}.png'.format(ii),dpi=300)
    plt.close()




plt.style.use('dark_background')
plt.figure(figsize=(10,6))
thetaMode2 =60
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])
thetaMode1 = np.arange(25,85,50)
ax1c = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))
#ax1c.legend()
# ax1c.set_title(r'3b. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))
ax1c.set_ylabel('Depth (MHHW, m)')
ax1c.set_xlabel('Cross-shore (m)')
plt.show()



plt.style.use('dark_background')
plt.figure(figsize=(10,6))
thetaMode2 =30
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])
thetaMode1 = np.arange(220,360,15)
ax1c = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))
#ax1c.legend()
# ax1c.set_title(r'3b. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))
ax1c.set_ylabel('Depth (MHHW, m)')
ax1c.set_xlabel('Cross-shore (m)')
plt.show()






plt.style.use('dark_background')
plt.figure(figsize=(8,5))
thetaMode2 =20
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])
thetaMode1 = np.arange(45,105,40)
ax1c = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)+2))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[ii+2,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]),linewidth=2)
#ax1c.legend()
# ax1c.set_title(r'3b. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))
# ax1c.set_ylabel('Depth (MHHW, m)')
# ax1c.set_xlabel('Cross-shore (m)')
ax1c.set_xlim([0,250])
ax1c.set_ylim([-5, 0])
plt.show()



plt.style.use('dark_background')
plt.figure(figsize=(8,5))
thetaMode2 =-70
eofPred2 = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,0])
thetaMode1 = np.arange(150,170,30)
ax1c = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)+2))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 1])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax1c.plot(xinterp,combined,color=colors[(len(colors)-2-ii),:],label=r'$\phi_1$ = {}'.format(thetaMode1[-ii]),linewidth=2)
#ax1c.legend()
# ax1c.set_title(r'3b. Stationary $\phi_2$ = {} + Propagating Offshore Signal'.format(thetaMode2))
# ax1c.set_ylabel('Depth (MHHW, m)')
# ax1c.set_xlabel('Cross-shore (m)')
ax1c.set_xlim([0,500])
ax1c.set_ylim([-7, 0])
plt.show()





























thetasMode1 = np.arange(0,1200,5)
thetasMode2 = np.arange(0,1200,5)

eofHypo = np.nan * np.ones((len(thetasMode1),len(xinterpS)))
eofHypo2 = np.nan * np.ones((len(thetasMode1),len(xinterpS)))

for tt in range(len(thetasMode1)):
   eofHypo[tt,:] = SS[:, 0] * np.nanmean(RS[:, 0]) * np.cos(thetasMode1[tt] * np.pi / 180 - thetaS[:, 0])
   eofHypo2[tt,:] = SS[:, 1] * np.nanmean(RS[:, 1]) * np.cos(thetasMode1[tt] * np.pi / 180 - thetaS[:, 1])


thetag, xintg = np.meshgrid(thetasMode1, xinterpS)

fig50 = plt.figure(figsize=(14,10))

ax1 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
ax1b = ax1.twiny()
plt0 = ax1.pcolor(thetag,xintg,eofHypo.T, vmin=-0.75, vmax=0.75, cmap='bwr')
plt0b = ax1b.plot(range(1200),np.ones(1200))
ax1b.cla()
ax1b.set_xlim(0,810)
# ax1.set_ylim(0,810)
ax1b.set_xticks([0,98,196,294,392,490,589,687,785])#,883,981,1079])
ax1b.set_xticklabels(['0','','1 year','','2 years','','3 years','','4 years'])
ax1b.set_xlabel('Time elapsed (Average Migration Speed)',fontsize=14)
ax1.set_title(r'Offshore Migration $\Longrightarrow$',fontsize=18)
ax1.set_ylabel('Cross-shore (m, distance from shoreline)',fontsize=14)
ax1.set_xlabel('Phase of Offshore Propagation (degrees)',fontsize=14)
# fig.colorbar(plt0, ax=ax1, orientation='horizontal')

ax2 = plt.subplot2grid((1,2),(0,1),rowspan=1,colspan=1)
ax2b = ax2.twiny()
plt1 = ax2.pcolor(thetag,xintg,eofHypo2.T, vmin=-0.75, vmax=0.75, cmap='bwr')
plt2b = ax2b.plot(range(810),np.ones(810))
ax2b.cla()
ax2b.set_xticks([0,132,263,394,525,656,788])
ax2b.set_xticklabels(['0','','1 year','','2 years','','3 years'])
ax2b.set_xlabel('Time elapsed (Average Migration Speed)',fontsize=14)
ax2.set_xlabel('Phase of Onshore Propagation (degrees)',fontsize=14)
ax2.set_title(r'Onshore Migration $\Longrightarrow$',fontsize=18)
ax2.set_ylim([0,700])
ax2b.set_ylim(0,700)

cb1 = fig.colorbar(plt1, ax=ax2)#, orientation='horizontal')
cb1.set_label('Difference from Mean Profile (m)',fontsize=14)
plt.show()


t1 = 0
t2 = -1
fig = plt.figure(figsize=(14,10))
#plt.set_cmap('RdBu')#bwr')
# plt.set_cmap('RdBu_r')

tg, xg = np.meshgrid(time, xinterp)
ax2 = plt.subplot2grid((7,3),(0,0),rowspan=7,colspan=1)
plt1 = ax2.pcolor(xg,tg,eofPredog.T, vmin=-.75, vmax=0.75,cmap='bwr')
ax2.set_ylim([time[t1], time[t2]])
#fig.colorbar(plt1, ax=ax2, orientation='horizontal')
ax2.set_title('Observed Offshore CEOF: {:.2f}%'.format(percentVS[0]))
# ax2.get_yaxis().set_ticks([])
ax2.set_xlabel('Cross-shore (m)')

ax3 = plt.subplot2grid((7,3),(0,1),rowspan=7,colspan=1)
plt2 = ax3.pcolor(xg,tg,eofPred2og.T, vmin=-.85, vmax=.85,cmap='bwr')
ax3.set_ylim([time[t1], time[t2]])
#fig.colorbar(plt2, ax=ax3, orientation='horizontal')
ax3.set_title('Observed Onshore CEOF: {:.2f}%'.format(percentVS[1]))
# ax3.get_yaxis().set_ticks([])
ax3.set_xlabel('Cross-shore (m)')



ax1 = plt.subplot2grid((7,3),(0,2),rowspan=3,colspan=1)
ax1b = ax1.twinx()
plt0 = ax1.pcolor(xintg,thetag,eofHypo.T, vmin=-0.75, vmax=0.75, cmap='bwr')
plt0b = ax1b.plot(range(1200),np.ones(1200))
ax1b.cla()
ax1b.set_ylim(0,810)
ax1.set_ylim(0,810)
ax1b.set_yticks([0,98,196,294,392,490,589,687,785])
ax1b.set_yticklabels(['0','','1 year','','2 years','','3 years','','4 years'])
ax1b.set_ylabel('Average Time Elapsed',fontsize=12)
ax1.set_title(r'Offshore Migration $\Longrightarrow$',fontsize=12)
# ax1.set_ylabel('Cross-shore (m, distance from shoreline)',fontsize=12)
ax1.set_ylabel('Phase (deg)',fontsize=12)
ax1.set_xlim(0,510)
# fig.colorbar(plt0, ax=ax1, orientation='horizontal')

ax22 = plt.subplot2grid((7,3),(3,2),rowspan=3,colspan=1)
ax22b = ax22.twinx()
plt12 = ax22.pcolor(xintg,thetag,eofHypo2.T, vmin=-0.75, vmax=0.75, cmap='bwr')
plt22b = ax22b.plot(range(1200),np.ones(1200))
ax22b.cla()
ax22b.set_yticks([0,132,263,394,525,656,788,920,1052])
ax22b.set_yticklabels(['0','','1 year','','2 years','','3 years','','4 years'])
ax22b.set_ylabel('Average Time Elapsed',fontsize=12)
ax22.set_ylabel('Phase (deg)',fontsize=12)
ax22.set_title(r'Onshore Migration $\Longrightarrow$',fontsize=12)
ax22.set_xlabel('Cross-shore (m, distance from shoreline)',fontsize=12)

ax22.set_ylim([0,1100])
ax22b.set_ylim(0,1100)
ax22.set_xlim([0,510])

#ax4 = plt.subplot2grid((7,3),(6,2),rowspan=1,colspan=1)
#cb = fig.colorbar(plt12,ax=ax4)

cbar_ax = fig.add_axes([0.71, 0.1, 0.2, 0.03])
cb1c = plt.colorbar(plt12, cax=cbar_ax, orientation='horizontal')
cb1c.set_label('Difference from Mean (m)')
plt.tight_layout()
plt.show()















# thetasMode1 = np.arange(-180,180,2)
# thetasMode2 = np.arange(-100,260,2)

thetasMode1 = np.arange(-270,1000,2)
thetasMode2 = np.arange(-270,1800,2)
innerBarX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
outerBarX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerBarZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
outerBarZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterBarX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterBarZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))

outerTroughZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
outerTroughX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerTroughZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerTroughX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterTroughZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterTroughX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))

innerBarDZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
outerBarDZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterBarDZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))

# innerBar = np.nan * np.ones((alllines.shape[0],))
# outerBar = np.nan * np.ones((alllines.shape[0],))
#
# zOuterBar = np.nan * np.ones((alllines.shape[0],))
# zInnerBar = np.nan * np.ones((alllines.shape[0],))
#
# innerTrough = np.nan * np.ones((alllines.shape[0],))
# outerTrough = np.nan * np.ones((alllines.shape[0],))
# zOuterTrough = np.nan * np.ones((alllines.shape[0],))
# zInnerTrough = np.nan * np.ones((alllines.shape[0],))
#

for tt in range(len(thetasMode1)):
   eofPred = SS[:, 0] * np.nanmean(RS[:, 0]) * np.cos(thetasMode1[tt] * np.pi / 180 - thetaS[:, 0])
   for yy in range(len(thetasMode2)):
      eofPred2 = SS[:, 1] * np.nanmean(RS[:, 1]) * np.cos(thetasMode2[yy] * np.pi / 180 - thetaS[:, 1])
      combined = np.real(np.mean(alllinesS,axis=0)+eofPred+eofPred2)
      fname = "/home/dylananderson/projects/duckGeomorph/sandBarTool/Mode1_{}_Mode2_{}.png".format(thetasMode1[tt],thetasMode2[yy])
      xFRFbar, xFRFtrough = findSandBarAndTrough1D(xinterpS, combined, plotFname=None, smoothLengthScale=5, profileTrend=np.mean(alllines,axis=0))
      if xFRFbar is not None:
            #for sandbarX in xFRFbar:
            #    ax[1].plot(sandbarX, time[tt], 'ro', label='bar')

         if len(xFRFbar) > 1:
            outerBarX[yy,tt] = np.max(xFRFbar)
            zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == xinterpS))
            outerBarZ[yy,tt] = combined[zOuterInd[0]]
            outerTroughZ[yy,tt] = np.max(xFRFtrough)
            outerTroughZInd = np.where((np.round(2*np.max(xFRFtrough))/2 == xinterpS))
            outerTroughZ[yy,tt] = combined[outerTroughZInd[0]]

            innerBarX[yy,tt] = np.min(xFRFbar)
            zInnerInd = np.where((np.round(2*np.min(xFRFbar))/2 == xinterpS))
            innerBarZ[yy,tt] = combined[zInnerInd[0]]
            innerTroughZ[yy,tt] = np.max(xFRFtrough)
            innerTroughZInd = np.where((np.round(2*np.min(xFRFtrough))/2 == xinterpS))
            innerTroughZ[yy,tt] = combined[innerTroughZInd[0]]

            outerBarDZ[yy,tt] = outerBarZ[yy,tt]-outerTroughZ[yy,tt]
            innerBarDZ[yy,tt] = innerBarZ[yy,tt]-innerTroughZ[yy,tt]

            # innerOuterBarX[yy,tt] = np.max(xFRFbar)
            # zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == xinterpS))
            # innerOuterBarZ[yy,tt] = combined[zOuterInd[0]]

         else:
            # innerBarX[yy,tt] = xFRFbar
            # zInnerInd = np.where((np.round(2*xFRFbar)/2 == xinterpS))
            # innerBarZ[yy,tt] = combined[zInnerInd[0]]
            #
            # outerBarX[yy,tt] = xFRFbar
            # zOuterInd = np.where((np.round(2*xFRFbar)/2 == xinterpS))
            # outerBarZ[yy,tt] = combined[zOuterInd[0]]

            innerOuterBarX[yy,tt] = xFRFbar
            zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == xinterpS))
            innerOuterBarZ[yy,tt] = combined[zOuterInd[0]]
            innerOuterTroughZ[yy, tt] = xFRFtrough
            innerOuterTroughZInd = np.where((np.round(2 * np.max(xFRFtrough)) / 2 == xinterpS))
            innerOuterTroughZ[yy, tt] = combined[innerOuterTroughZInd[0]]
            innerOuterBarDZ[yy,tt] = innerOuterBarZ[yy,tt]-innerOuterTroughZ[yy,tt]


      # if xFRFtrough is not None:
        #     #for troughX in xFRFtrough:
        #     #    ax[1].plot(troughX, time[tt], 'bd', label='trough')
        #     if len(xFRFtrough) > 1:
        #         outerTrough[tt] = np.max(xFRFtrough)
        #         zOuterInd = np.where((np.round(2*np.max(xFRFtrough))/2 == bathyX))
        #         zOuterTrough[tt] = bathy[zOuterInd[0]]
        #
        #         innerTrough[tt] = np.min(xFRFtrough)
        #         zInnerInd = np.where((np.round(2*np.min(xFRFtrough))/2 == bathyX))
        #         zInnerTrough[tt] = bathy[zInnerInd[0]]
        #     else:
        #         outerTrough[tt] = xFRFtrough
        #         zOuterInd = np.where((np.round(2*xFRFtrough)/2 == bathyX))
        #         zOuterTrough[tt] = bathy[zOuterInd[0]]
        #barPatch = mpatches.Patch(color='red', label='sandbar')
        #troughPatch = mpatches.Patch(color='blue', label='trough')

meshTheta1, meshTheta2 = np.meshgrid(thetasMode1,thetasMode2)

# fig2 = plt.figure(figsize=(12,10))
# ax1 = plt.subplot2grid((2,2),(0,0),rowspan=1,colspan=1)
# p1 = ax1.pcolor(meshTheta1,meshTheta2,innerBarX,vmin=30,vmax=300)
# cb1 = plt.colorbar(p1,ax=ax1)
# cb1.set_label('xFRF')
# ax2 = plt.subplot2grid((2,2),(0,1),rowspan=1,colspan=1)
# p2 = ax2.pcolor(meshTheta1,meshTheta2,outerBarX,vmin=30,vmax=300)
# ax2.set_xlabel('Mode 1 Phase')
# ax1.set_xlabel('Mode 1 Phase')
# ax1.set_ylabel('Mode 2 Phase')
# ax1.set_title('Inner Bar Cross-shore Location')
# ax2.set_title('Outer Bar Cross-shore Location')
# cb2 = plt.colorbar(p2,ax=ax2)
# cb2.set_label('xFRF')
#
#
# ax3 = plt.subplot2grid((2,2),(1,0),rowspan=1,colspan=1)
# p3 = ax3.pcolor(meshTheta1,meshTheta2,innerBarDZ)#,vmin=-6,vmax=-1.5)
# cb3 = plt.colorbar(p3,ax=ax3)
# cb3.set_label('Depth (m)')
#
# ax4 = plt.subplot2grid((2,2),(1,1),rowspan=1,colspan=1)
# p4 = ax4.pcolor(meshTheta1,meshTheta2,outerBarDZ)#,vmin=-6,vmax=-1.5)
# cb4 = plt.colorbar(p4,ax=ax4)
# cb4.set_label('Depth (m)')
#
# ax4.set_xlabel('Mode 1 Phase')
# ax3.set_xlabel('Mode 1 Phase')
# ax3.set_ylabel('Mode 2 Phase')
# ax3.set_title('Inner Bar Depth')
# ax4.set_title('Outer Bar Depth')
# plt.tight_layout()
#
# plt.show()

# index1 = np.where((timeS > DT.datetime(1984,9,1)) & (timeS < DT.datetime(1988,9,1)))
# indexedTime = timeS[index1]
# indexedOrdinal = [i.toordinal() for i in indexedTime]
# plt.figure()
# mode1prop = np.unwrap(phitS[index1[0],0])
# mode2prop = np.unwrap(phitS[index1[0],1])
# plt.plot(mode1prop,mode2prop)




fig2 = plt.figure(figsize=(14,8))
ax1 = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
p1 = ax1.pcolor(meshTheta1,meshTheta2,innerBarX,vmin=30,vmax=300)
cb1 = plt.colorbar(p1,ax=ax1)
cb1.set_label('xFRF')
ax2 = plt.subplot2grid((2,3),(0,1),rowspan=1,colspan=1)
p2 = ax2.pcolor(meshTheta1,meshTheta2,outerBarX,vmin=30,vmax=300)
ax2.set_xlabel('Mode 1 Phase')
ax1.set_xlabel('Mode 1 Phase')
# ax3.set_xlabel('Mode 1 Phase')
ax1.set_ylabel('Mode 2 Phase')
ax1.set_title('Inner Bar Cross-shore Location')
ax2.set_title('Outer Bar Cross-shore Location')

cb2 = plt.colorbar(p2,ax=ax2)
cb2.set_label('xFRF')
ax1c = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
p1c = ax1c.pcolor(meshTheta1,meshTheta2,innerOuterBarX,vmin=30,vmax=300)
cb1c = plt.colorbar(p1c,ax=ax1c)
cb1c.set_label('xFRF')
ax1c.set_title('Single Bar Cross-shore Location')

ax3 = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
p3 = ax3.pcolor(meshTheta1,meshTheta2,innerBarDZ,vmin=0,vmax=1)
cb3 = plt.colorbar(p3,ax=ax3)
cb3.set_label('Depth (m)')

ax4 = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
p4 = ax4.pcolor(meshTheta1,meshTheta2,outerBarDZ,vmin=0,vmax=1)
cb4 = plt.colorbar(p4,ax=ax4)
cb4.set_label('Depth (m)')

ax4c = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
p4c = ax4c.pcolor(meshTheta1,meshTheta2,innerOuterBarDZ,vmin=0,vmax=1)
cb4c = plt.colorbar(p4c,ax=ax4c)
cb4c.set_label('Depth (m)')

ax4c.set_xlabel('Mode 1 Phase')
ax4.set_xlabel('Mode 1 Phase')
ax3.set_xlabel('Mode 1 Phase')
ax3.set_ylabel('Mode 2 Phase')
ax3.set_title('Inner Bar Magnitude')
ax4.set_title('Outer Bar Magnitude')
ax4c.set_title('Single Bar Magnitude')

plt.tight_layout()


plt.figure(figsize=(14,10))
ax101 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
p101 = ax101.pcolor(meshTheta1,meshTheta2,innerBarZ-outerBarZ)#vmin=0,vmax=1)
cb101 = plt.colorbar(p101,ax=ax101)
ax101.set_xlabel(r'Mode 2 Phase (offshore $\longrightarrow$)')
ax101.set_ylabel(r'Mode 2 Phase (onshore $\longrightarrow$)')
cb101.set_label('Difference in Bar Crest Depth')
ax101.set_title('Onshore Propagation Separates the Crest''s Depths')
ax101.set_xlim([-180,180])
ax101.set_ylim([-100,260])

ax102 = plt.subplot2grid((1,2),(0,1),rowspan=1,colspan=1)
p102 = ax102.pcolor(meshTheta1,meshTheta2,outerBarX-innerBarX)#vmin=0,vmax=1)
cb102 = plt.colorbar(p102,ax=ax102)
ax102.set_xlabel(r'Mode 2 Phase (offshore $\longrightarrow$)')
ax102.set_ylabel(r'Mode 2 Phase (onshore $\longrightarrow$)')
cb102.set_label('Difference in Bar Crest Cross-shore Location')
ax102.set_title('Offshore Propagation Converges the Crest''s Horizontal Distance')
ax102.set_xlim([-180,180])
ax102.set_ylim([-100,260])

plt.show()

# import matplotlib
# params = {
#     'text.latex.preamble': ['\\usepackage{gensymb}'],
#     'image.origin': 'lower',
#     'image.interpolation': 'nearest',
#     'image.cmap': 'gray',
#     'axes.grid': False,
#     'savefig.dpi': 150,  # to adjust notebook inline plot size
#     'axes.labelsize': 8, # fontsize for x and y labels (was 10)
#     'axes.titlesize': 8,
#     'font.size': 8, # was 10
#     'legend.fontsize': 6, # was 10
#     'xtick.labelsize': 8,
#     'ytick.labelsize': 8,
#     'text.usetex': True,
#     'figure.figsize': [3.39, 2.10],
#     'font.family': 'serif',
# }
# matplotlib.rcParams.update('default')
#

import copy
doublebar = copy.deepcopy(outerBarX)
doubleIndex = np.where((doublebar > 0))
doublebar[doubleIndex] = 0.5

singlebar = copy.deepcopy(innerOuterBarX)
singleIndex = np.where((singlebar > 0))
doublebar[singleIndex] = 0.75


plt.style.use('dark_background')
fig = plt.figure(figsize=(12,16))
ax1 = plt.subplot2grid((4,3),(0,0),rowspan=1,colspan=1)
p1 = ax1.pcolor(meshTheta1,meshTheta2,innerBarX,vmin=30,vmax=300,cmap='plasma')
ax1.set_xlim([-180,180])
ax1.set_ylim([-100,260])

# cb1 = plt.colorbar(p1,ax=ax1)
# cb1.set_label('xFRF')
# ax1.text(50,200,'Location')

ax2 = plt.subplot2grid((4,3),(0,1),rowspan=1,colspan=1)
p2 = ax2.pcolor(meshTheta1,meshTheta2,outerBarX,vmin=30,vmax=300,cmap='plasma')
ax2.set_xlim([-180,180])
ax2.set_ylim([-100,260])

# ax2.text(50,200,'Location')

# ax2.set_xlabel('Mode 1 Phase')
# ax1.set_xlabel('Mode 1 Phase')
# ax3.set_xlabel('Mode 1 Phase')
# ax1.set_ylabel('Mode 2 Phase')
ax1.set_ylabel(r'Onshore Propagation (deg) $\Longrightarrow$')

ax1.set_title('Inner Bar',fontsize=16)
ax2.set_title('Outer Bar',fontsize=16)

# cb2 = plt.colorbar(p2,ax=ax2)
# cb2.set_label('xFRF')
ax1c = plt.subplot2grid((4,3),(0,2),rowspan=1,colspan=1)
p1c = ax1c.pcolor(meshTheta1,meshTheta2,innerOuterBarX,vmin=30,vmax=300,cmap='plasma')
ax1c.set_xlim([-180,180])
ax1c.set_ylim([-100,260])


# cb1c = plt.colorbar(p1c,ax=ax1c)
# cb1c.set_label('xFRF')
ax1c.set_title('Single Bar',fontsize=16)

ax3 = plt.subplot2grid((4,3),(1,0),rowspan=1,colspan=1)
p3 = ax3.pcolor(meshTheta1,meshTheta2,innerBarDZ,vmin=0,vmax=1)
ax3.set_xlim([-180,180])
ax3.set_ylim([-100,260])

# cb3 = plt.colorbar(p3,ax=ax3)
# cb3.set_label('Depth (m)')
# ax3.text(50,200,'Magnitude')

ax4 = plt.subplot2grid((4,3),(1,1),rowspan=1,colspan=1)
p4 = ax4.pcolor(meshTheta1,meshTheta2,outerBarDZ,vmin=0,vmax=1)
ax4.set_xlim([-180,180])
ax4.set_ylim([-100,260])

# cb4 = plt.colorbar(p4,ax=ax4)
# cb4.set_label('Depth (m)')
# ax4.text(50,200,'Magnitude')

ax4c = plt.subplot2grid((4,3),(1,2),rowspan=1,colspan=1)
p4c = ax4c.pcolor(meshTheta1,meshTheta2,innerOuterBarDZ,vmin=0,vmax=1)
ax4c.set_xlim([-180,180])
ax4c.set_ylim([-100,260])

# cb4c = plt.colorbar(p4c,ax=ax4c)
# cb4c.set_label('Depth (m)')
#
# ax101 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=1)
# p101 = ax101.pcolor(meshTheta1,meshTheta2,innerBarZ-outerBarZ,cmap='inferno')#vmin=0,vmax=1)
# # cb101 = plt.colorbar(p101,ax=ax101)
#
# ax101 = plt.subplot2grid((3,3),(2,1),rowspan=1,colspan=1)
# p101 = ax101.pcolor(meshTheta1,meshTheta2,outerBarX-innerBarX,cmap='cividis')#vmin=0,vmax=1)
# # cb101 = plt.colorbar(p101,ax=ax101)


ax4c.set_xlabel(r'Offshore Propagation (deg) $\Longrightarrow$')
ax4.set_xlabel(r'Offshore Propagation (deg) $\Longrightarrow$')
ax3.set_xlabel(r'Offshore Propagation (deg) $\Longrightarrow$')
ax3.set_ylabel(r'Onshore Propagation (deg) $\Longrightarrow$')
# ax3.set_title('Inner Bar Magnitude')
# ax4.set_title('Outer Bar Magnitude')
# ax4c.set_title('Single Bar Magnitude')

plt.subplots_adjust(right=0.87)
cbar_ax = fig.add_axes([0.91, 0.72, 0.02, 0.15])
cb1c = plt.colorbar(p1c,cax=cbar_ax)
cb1c.set_label('Distance from shoreline (m)')

cbar_ax2 = fig.add_axes([0.91, 0.52, 0.02, 0.15])
cb1c2 = plt.colorbar(p4c,cax=cbar_ax2)
cb1c2.set_label('Bar Magnitude (m)')

conceptbar_ax = fig.add_axes([0.3, 0.05, 0.4, 0.4])
con = conceptbar_ax.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
conceptbar_ax.set_xlabel(r'Offshore Propagation (deg) $\Longrightarrow$')
conceptbar_ax.set_ylabel(r'Onshore Propagation (deg) $\Longrightarrow$')
conceptbar_ax.set_xlim([-180,180])
conceptbar_ax.set_ylim([-100,260])

# cb4c = plt.colorbar(p4c,ax=ax4c)
# cb4c.set_label('Depth (m)')

# plt.tight_layout()

plt.show()






conceptualFig = plt.figure(figsize=(10,12))
plt.style.use('default')
axLeft = plt.subplot2grid((7,6),(4,0),rowspan=3,colspan=3)
con = axLeft.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
axLeft.set_xlabel(r'Offshore Propagation (deg) $\Longrightarrow$')
axLeft.set_ylabel(r'Onshore Propagation (deg) $\Longrightarrow$')
axLeft.set_xlim([-180,180])
axLeft.set_ylim([-100,260])

axRight = plt.subplot2grid((7,6),(4,3),rowspan=3,colspan=3)
con = axRight.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
axRight.set_xlabel(r'Offshore Propagation (deg) $\Longrightarrow$')
axRight.set_ylabel(r'Onshore Propagation (deg) $\Longrightarrow$')
axRight.set_xlim([-180,180])
axRight.set_ylim([-100,260])
x = np.arange(-110,30,1)
axRight.plot(x ,45*np.ones(len(x)),'g--',label='Outer Bar Decay')
y = np.arange(145,290,1)
axRight.plot(-90*np.ones(len(y)),y,'r--')
axRight.plot(-90*np.ones(len(y)),y-360,'r--')
y2 = np.arange(0,180,1)
axRight.plot(-80*np.ones(len(y2)),y2,'k--',label='2-bar to low-tide terrace')
x2 = np.arange(-0,135,1)
axRight.plot(x2 ,25*np.ones(len(x2)),'m--',label='Outer Bar Decay')

ax1 = plt.subplot2grid((7,6),(0,0),rowspan=2,colspan=2)
p1 = ax1.pcolor(meshTheta1,meshTheta2,innerBarX,vmin=30,vmax=300,cmap='plasma')
ax1.set_xlim([-180,180])
ax1.set_ylim([-100,260])

ax2 = plt.subplot2grid((7,6),(0,2),rowspan=2,colspan=2)
p2 = ax2.pcolor(meshTheta1,meshTheta2,outerBarX,vmin=30,vmax=300,cmap='plasma')
ax2.set_xlim([-180,180])
ax2.set_ylim([-100,260])

ax1.set_ylabel(r'Onshore Propagation (deg) $\Longrightarrow$')

ax1.set_title('Inner Bar',fontsize=16)
ax2.set_title('Outer Bar',fontsize=16)

ax1c = plt.subplot2grid((7,6),(0,4),rowspan=2,colspan=2)
p1c = ax1c.pcolor(meshTheta1,meshTheta2,innerOuterBarX,vmin=30,vmax=300,cmap='plasma')
ax1c.set_xlim([-180,180])
ax1c.set_ylim([-100,260])

ax1c.set_title('Single Bar',fontsize=16)

ax3 = plt.subplot2grid((7,6),(2,0),rowspan=2,colspan=2)
p3 = ax3.pcolor(meshTheta1,meshTheta2,innerBarDZ,vmin=0,vmax=1)
ax3.set_xlim([-180,180])
ax3.set_ylim([-100,260])
ax4 = plt.subplot2grid((7,6),(2,2),rowspan=2,colspan=2)
p4 = ax4.pcolor(meshTheta1,meshTheta2,outerBarDZ,vmin=0,vmax=1)
ax4.set_xlim([-180,180])
ax4.set_ylim([-100,260])
ax4c = plt.subplot2grid((7,6),(2,4),rowspan=2,colspan=2)
p4c = ax4c.pcolor(meshTheta1,meshTheta2,innerOuterBarDZ,vmin=0,vmax=1)

ax4c.set_xlim([-180,180])
ax4c.set_ylim([-100,260])
ax4c.set_xlabel(r'Offshore Propagation (deg) $\Longrightarrow$')
ax4.set_xlabel(r'Offshore Propagation (deg) $\Longrightarrow$')
ax3.set_xlabel(r'Offshore Propagation (deg) $\Longrightarrow$')
ax3.set_ylabel(r'Onshore Propagation (deg) $\Longrightarrow$')
plt.tight_layout()
plt.show()















# for ii in range(11): #range(numClusters):
#     finder = orderC[ii]
#     #index = np.where(yc == finder)
#     if flatarrayC[ii,ii] < 0.55:
#         sizeMarker = 50
#     elif flatarrayC[ii,ii] < 0.70 and flatarrayC[ii,ii]>0.549:
#         sizeMarker = 125
#     else:
#         sizeMarker = 200
#
#     if ii == 5:
#         ax4c.scatter(orderAngleC[finder],orderAngle2C[finder]+360,sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#     elif ii == 6:
#         ax4c.scatter(orderAngleC[finder],orderAngle2C[finder],sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#         ax4c.scatter(orderAngleC[finder]+360,orderAngle2C[finder]+360,sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#     elif ii == 7:
#         ax4c.scatter(orderAngleC[finder],orderAngle2C[finder],sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#         ax4c.scatter(orderAngleC[finder]+360,orderAngle2C[finder]+360,sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#     elif ii == 3:
#         ax4c.scatter(orderAngleC[finder],orderAngle2C[finder],sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#         ax4c.scatter(orderAngleC[finder]-360,orderAngle2C[finder],sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#     elif ii == 1:
#         ax4c.scatter(orderAngleC[finder],orderAngle2C[finder],sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#         ax4c.scatter(orderAngleC[finder],orderAngle2C[finder]+360,sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#     else:
#         ax4c.scatter(orderAngleC[finder],orderAngle2C[finder],sizeMarker,color=colorsC[ii,:],edgecolor='k',zorder=3)
#
#
# plt.show()
#plt.plot(time,innerBarX,'bo')
#plt.plot(time,outerBarX,'ro')




# outerBarSurveyX = np.nan*np.ones((len(alllinesS),))
# outerBarSurveyZ = np.nan*np.ones((len(alllinesS),))
# innerBarSurveyX = np.nan*np.ones((len(alllinesS),))
# innerBarSurveyZ = np.nan*np.ones((len(alllinesS),))
# outerInnerBarSurveyX = np.nan*np.ones((len(alllinesS),))
# outerInnerBarSurveyZ = np.nan*np.ones((len(alllinesS),))
#
# outerTroughSurveyX = np.nan*np.ones((len(alllinesS),))
# outerTroughSurveyZ = np.nan*np.ones((len(alllinesS),))
# innerTroughSurveyX = np.nan*np.ones((len(alllinesS),))
# innerTroughSurveyZ = np.nan*np.ones((len(alllinesS),))
# outerInnerTroughSurveyX = np.nan*np.ones((len(alllinesS),))
# outerInnerTroughSurveyZ = np.nan*np.ones((len(alllinesS),))
#
#
# for tt in range(alllinesS.shape[0]):
#     bathy = alllinesS[tt,:]
#     bathyX = xinterpS
#     if bathy[~np.isnan(bathy)].any():
#
#         fname = "/home/dylananderson/projects/duckGeomorph/sandBarTool/Line1_{}.png".format(tt)
#
#         #xFRFbar, xFRFtrough = mL.findSandBarAndTrough1D(bathyX, bathy, plotFname=fname, smoothLengthScale=10, profileTrend=np.mean(alllines,axis=0))
#         xFRFbar, xFRFtrough = mL.findSandBarAndTrough1D(bathyX, bathy, plotFname=None, smoothLengthScale=2, profileTrend=np.mean(alllinesS,axis=0))
#
#         if xFRFbar is not None:
#             #for sandbarX in xFRFbar:
#             #    ax[1].plot(sandbarX, time[tt], 'ro', label='bar')
#
#             if len(xFRFbar) > 1:
#                 outerBarSurveyX[tt] = np.max(xFRFbar)
#                 zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == bathyX))
#                 outerBarSurveyZ[tt] = bathy[zOuterInd[0]]
#
#                 innerBarSurveyX[tt] = np.min(xFRFbar)
#                 zInnerInd = np.where((np.round(2*np.min(xFRFbar))/2 == bathyX))
#                 innerBarSurveyZ[tt] = bathy[zInnerInd[0]]
#
#             else:
#                 outerInnerBarSurveyX[tt] = xFRFbar
#                 zOuterInd = np.where((np.round(2*xFRFbar)/2 == bathyX))
#                 outerInnerBarSurveyZ[tt] = bathy[zOuterInd[0]]
#
#         if xFRFtrough is not None:
#             #for troughX in xFRFtrough:
#             #    ax[1].plot(troughX, time[tt], 'bd', label='trough')
#             if len(xFRFtrough) > 1:
#                 outerTroughSurveyX[tt] = np.max(xFRFtrough)
#                 zOuterInd = np.where((np.round(2*np.max(xFRFtrough))/2 == bathyX))
#                 outerTroughSurveyZ[tt] = bathy[zOuterInd[0]]
#
#                 innerTroughSurveyX[tt] = np.min(xFRFtrough)
#                 zInnerInd = np.where((np.round(2*np.min(xFRFtrough))/2 == bathyX))
#                 innerTroughSurveyZ[tt] = bathy[zInnerInd[0]]
#             else:
#                 outerInnerTroughSurveyX[tt] = xFRFtrough
#                 zOuterInd = np.where((np.round(2*xFRFtrough)/2 == bathyX))
#                 outerInnerTroughSurveyZ[tt] = bathy[zOuterInd[0]]
#         #barPatch = mpatches.Patch(color='red', label='sandbar')
#         #troughPatch = mpatches.Patch(color='blue', label='trough')
#
#
# fig2 = plt.figure(figsize=(8,8))
# plt.plot(time,innerBarSurveyX,'bo')
# plt.plot(time,outerBarSurveyX,'ro')
# plt.show()













mode = 0

fig = plt.figure(figsize=(10,5))
ax1 = plt.subplot2grid((3,2),(0,0),rowspan=1,colspan=2)
xg,tg = np.meshgrid(xinterpS,timeS)
plt0 = ax1.pcolor(tg,xg,(alllinesS-np.mean(alllinesS, axis=0)), vmin=-1.8, vmax=1.8,cmap='bwr')

ax = plt.subplot2grid((3,2),(1,0),rowspan=1,colspan=2,sharex=ax1)
temp = RS[:,0]*10
temp1 = phit2S[:,0]
ax.scatter(timeS,temp1,temp,temp,cmap='Reds')

ax2 = plt.subplot2grid((3,2),(2,0),rowspan=1,colspan=2, sharex=ax1)
temp = RS[:,1]*10
temp1 = phit2S[:,1]
ax2.scatter(timeS,temp1,temp,temp,cmap='Blues')

ax.set_xlim([timeS[0], timeS[-1]])
ax1.set_xlim([timeS[0], timeS[-1]])
ax2.set_xlim([timeS[0], timeS[-1]])



ordinalDates = [i.toordinal() for i in timeS]
from scipy import stats
def myfunc(x):
    return slope * x + intercept

# firstSlopeInd = np.where((timeS > DT.datetime(2016,11,1)) & (timeS < DT.datetime(2019,11,1)))
firstSlopeInd = np.where((timeS > DT.datetime(1981,2,1)) & (timeS < DT.datetime(1983,12,1)))
x = np.array(ordinalDates)[firstSlopeInd]
y = phit2S[firstSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))

secondSlopeInd = np.where((timeS > DT.datetime(1987,12,1)) & (timeS < DT.datetime(1990,12,1)))
x = np.array(ordinalDates)[secondSlopeInd]
y = phit2S[secondSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))

thirdSlopeInd = np.where((timeS > DT.datetime(1990,12,1)) & (timeS < DT.datetime(1992,12,1)))
x = np.array(ordinalDates)[thirdSlopeInd]
y = phit2S[thirdSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))

fourthSlopeInd = np.where((timeS > DT.datetime(1992,12,1)) & (timeS < DT.datetime(1997,10,1)))
x = np.array(ordinalDates)[fourthSlopeInd]
y = phit2S[fourthSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))


fifthSlopeInd = np.where((timeS > DT.datetime(1997,12,15)) & (timeS < DT.datetime(1998,12,25)))
x = np.array(ordinalDates)[fifthSlopeInd]
y = phit2S[fifthSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))


sixthSlopeInd = np.where((timeS > DT.datetime(1984,9,15)) & (timeS < DT.datetime(1987,12,1)))
x = np.array(ordinalDates)[sixthSlopeInd]
y = phit2S[sixthSlopeInd,0]
z = np.real(RS[sixthSlopeInd,0])
largerMags = np.where((z > 0.4))

slope, intercept, r, p, std_err = stats.linregress(x[largerMags[1]], y[0,largerMags[1]])
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))

seventhSlopeInd = np.where((timeS > DT.datetime(2015,12,15)) & (timeS < DT.datetime(2018,12,1)))
x = np.array(ordinalDates)[seventhSlopeInd]
y = phit2S[seventhSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))

eigthSlopeInd = np.where((timeS > DT.datetime(2011,2,15)) & (timeS < DT.datetime(2015,1,1)))
x = np.array(ordinalDates)[eigthSlopeInd]
# y = np.unwrap(phit2S[eigthSlopeInd,0]*np.pi/180)*180/np.pi
y = phit2S[eigthSlopeInd,0]
z = np.real(RS[eigthSlopeInd,0])
largerMags = np.where((z > 0.4))

slope, intercept, r, p, std_err = stats.linregress(x[largerMags[1]], y[0,largerMags[1]])
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))


ninthSlopeInd = np.where((timeS > DT.datetime(2002,9,15)) & (timeS < DT.datetime(2004,1,1)))
x = np.array(ordinalDates)[ninthSlopeInd]
y = phit2S[ninthSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))

tenthSlopeInd = np.where((timeS > DT.datetime(2019,9,15)) & (timeS < DT.datetime(2020,8,1)))
x = np.array(ordinalDates)[tenthSlopeInd]
y = phit2S[tenthSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))

elevethSlopeInd = np.where((timeS > DT.datetime(1999,3,1)) & (timeS < DT.datetime(2001,12,1)))
x = np.array(ordinalDates)[elevethSlopeInd]
y = phit2S[elevethSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))

twelveSlopeInd = np.where((timeS > DT.datetime(2006,3,1)) & (timeS < DT.datetime(2009,6,1)))
x = np.array(ordinalDates)[twelveSlopeInd]
y = phit2S[twelveSlopeInd,0]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],350,'{}'.format(np.round(slope*365)))

# thirteenSlopeInd = np.where((timeS > DT.datetime(200,3,1)) & (timeS < DT.datetime(2009,6,1)))
# x = np.array(ordinalDates)[thirteenSlopeInd]
# y = phit2S[thirteenSlopeInd,0]
# slope, intercept, r, p, std_err = stats.linregress(x, y)
# mymodel = list(map(myfunc, x))
# ax.plot(x, mymodel ,'k')

plt.show()


# plt.figure()
# # plt.plot(np.cos(phit2S[:,0]*np.pi/180),np.cos(phit2S[:,1]*np.pi/180),'.')
# plt.polar(phit2S[:,0],RS[:,0])


def moving_average(a, n=21):
   ret = np.cumsum(a, dtype=float)
   ret[n:] = ret[n:] - ret[:-n]
   return ret[n - 1:] / n



def weightedMovingAverage(a,b,n=3):
    cut = np.floor(n/2)
    index = np.arange(int(cut),int(len(a)-cut),1)
    output = np.nan * np.ones((len(index),))
    counter = 0
    for ff in index:
        subset = a[int(ff-cut):int(ff+cut+1)]
        weights = (b[int(ff-cut):int(ff+cut+1)])
        output[counter] = np.average(subset,weights=weights)
        counter = counter+1
    return output



# mode1Phase = moving_average(np.unwrap(temp1*(np.pi/180)),n=5)
mode1Phase = moving_average(np.unwrap(phit2S[:,0]*(np.pi/180)),n=5)
mode2Phase = moving_average(np.unwrap(phit2S[:,1]*(np.pi/180)),n=5)
mode3Phase = moving_average(np.unwrap(phit2S[:,2]*(np.pi/180)),n=5)

# a = np.unwrap(phit2S[:,0]*(np.pi/180))
# b = RS[:,0]
# n = 5
# mode1PhaseTest = weightedMovingAverage(a,b,n)


mode1Mag = moving_average(RS[:,0],n=5)
mode2Mag = moving_average(RS[:,1],n=5)
mode3Mag = moving_average(RS[:,2],n=5)

mode1Wrap = (mode1Phase + np.pi) % (2 * np.pi) - np.pi
mode2Wrap = (mode2Phase + np.pi) % (2 * np.pi) - np.pi
mode3Wrap = (mode3Phase + np.pi) % (2 * np.pi) - np.pi

index = np.where((mode1Wrap<0))
mode1Wrap[index[0]] = mode1Wrap[index[0]]+np.pi*2
# phases =


plt.figure()
plt.plot(timeS[2:-2],mode1Phase)









mode = 0

fig = plt.figure(figsize=(10,5))
ax1 = plt.subplot2grid((3,2),(0,0),rowspan=1,colspan=2)
xg,tg = np.meshgrid(xinterpS,timeS)
plt0 = ax1.pcolor(tg,xg,(alllinesS-np.mean(alllinesS, axis=0)), vmin=-1.8, vmax=1.8,cmap='bwr')

ax = plt.subplot2grid((3,2),(1,0),rowspan=1,colspan=2,sharex=ax1)
# temp = RS[:,mode]*10
# temp1 = phit2S[:,mode]
temp = mode1Mag*10
temp1 = mode1Wrap
# index = np.where((temp1<0))
# temp1[index[0]] = temp1[index[0]]+360
# ax.scatter(timeS,np.unwrap(phit2S[:,mode]*(np.pi/180)),temp,phit2S[:,mode])
ax.scatter(timeS[2:-2],temp1,temp,mode1Mag,cmap='Reds')

mode = 1
ax2 = plt.subplot2grid((3,2),(2,0),rowspan=1,colspan=2, sharex=ax1)
# temp = RS[:,1]*10
# temp1 = phit2S[:,mode]
temp = mode2Mag*10
temp1 = mode2Wrap
# index = np.where((temp1<0))
# temp1[index[0]] = temp1[index[0]]+360
ax2.scatter(timeS[2:-2],temp1,temp,mode2Mag,cmap='Blues')

ax.set_xlim([timeS[0], timeS[-1]])
ax1.set_xlim([timeS[0], timeS[-1]])
ax2.set_xlim([timeS[0], timeS[-1]])
plt.show()

timeS2 = timeS[2:-2]
ordinalDates = [i.toordinal() for i in timeS2]
from scipy import stats
def myfunc(x):
    return slope * x + intercept

# firstSlopeInd = np.where((timeS > DT.datetime(2016,11,1)) & (timeS < DT.datetime(2019,11,1)))
firstSlopeInd = np.where((timeS2 > DT.datetime(1981,2,1)) & (timeS2 < DT.datetime(1983,12,1)))
x = np.array(ordinalDates)[firstSlopeInd]
y = mode1Wrap[firstSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

secondSlopeInd = np.where((timeS2 > DT.datetime(1987,12,1)) & (timeS2 < DT.datetime(1990,9,1)))
x = np.array(ordinalDates)[secondSlopeInd]
y = mode1Wrap[secondSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

thirdSlopeInd = np.where((timeS2 > DT.datetime(1990,9,1)) & (timeS2 < DT.datetime(1992,11,1)))
x = np.array(ordinalDates)[thirdSlopeInd]
y = mode1Wrap[thirdSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

fourthSlopeInd = np.where((timeS2 > DT.datetime(1992,10,15)) & (timeS2 < DT.datetime(1993,10,15)))
x = np.array(ordinalDates)[fourthSlopeInd]
y = mode1Wrap[fourthSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365)))

fourthSlopeIndb = np.where((timeS2 > DT.datetime(1993,10,20)) & (timeS2 < DT.datetime(1997,10,15)))
x = np.array(ordinalDates)[fourthSlopeIndb]
y = mode1Wrap[fourthSlopeIndb]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

fifthSlopeInd = np.where((timeS2 > DT.datetime(1998,1,1)) & (timeS2 < DT.datetime(1998,12,25)))
x = np.array(ordinalDates)[fifthSlopeInd]
y = mode1Wrap[fifthSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))


sixthSlopeInd = np.where((timeS2 > DT.datetime(1984,1,1)) & (timeS2 < DT.datetime(1987,12,1)))
x = np.array(ordinalDates)[sixthSlopeInd]
y = mode1Wrap[sixthSlopeInd]
z = np.real(mode1Mag[sixthSlopeInd])
largerMags = np.where((z > 0.4))
slope, intercept, r, p, std_err = stats.linregress(x[largerMags[0]], y[largerMags[0]])

# slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

seventhSlopeInd = np.where((timeS2 > DT.datetime(2015,12,15)) & (timeS2 < DT.datetime(2018,12,1)))
x = np.array(ordinalDates)[seventhSlopeInd]
y = mode1Wrap[seventhSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

eigthSlopeInd = np.where((timeS2 > DT.datetime(2011,2,15)) & (timeS2 < DT.datetime(2015,12,1)))
x = np.array(ordinalDates)[eigthSlopeInd]
# y = np.unwrap(phit2S[eigthSlopeInd,0]*np.pi/180)*180/np.pi
y = mode1Wrap[eigthSlopeInd]
z = np.real(mode1Mag[eigthSlopeInd])
largerMags = np.where((z > 0.45))
slope, intercept, r, p, std_err = stats.linregress(x[largerMags[0]], y[largerMags[0]])

# slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))


ninthSlopeInd = np.where((timeS2 > DT.datetime(2002,9,15)) & (timeS2 < DT.datetime(2004,1,1)))
x = np.array(ordinalDates)[ninthSlopeInd]
y = mode1Wrap[ninthSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

tenthSlopeInd = np.where((timeS2 > DT.datetime(2019,9,15)) & (timeS2 < DT.datetime(2020,8,1)))
x = np.array(ordinalDates)[tenthSlopeInd]
y = mode1Wrap[tenthSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

elevethSlopeInd = np.where((timeS2 > DT.datetime(1999,3,1)) & (timeS2 < DT.datetime(2001,12,1)))
x = np.array(ordinalDates)[elevethSlopeInd]
y = mode1Wrap[elevethSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

twelveSlopeInd = np.where((timeS2 > DT.datetime(2006,3,1)) & (timeS2 < DT.datetime(2009,6,1)))
x = np.array(ordinalDates)[twelveSlopeInd]
y = mode1Wrap[twelveSlopeInd]
slope, intercept, r, p, std_err = stats.linregress(x, y)
mymodel = list(map(myfunc, x))
ax.plot(x, mymodel ,'k')
ax.text(x[0],6,'{}'.format(np.round(slope*365*100)/100))

# thirteenSlopeInd = np.where((timeS > DT.datetime(200,3,1)) & (timeS < DT.datetime(2009,6,1)))
# x = np.array(ordinalDates)[thirteenSlopeInd]
# y = phit2S[thirteenSlopeInd,0]
# slope, intercept, r, p, std_err = stats.linregress(x, y)
# mymodel = list(map(myfunc, x))
# ax.plot(x, mymodel ,'k')

plt.show()


#
# def weighted_moving_average(x,y,step_size=0.05,width=1):
#     bin_centers  = np.arange(np.min(x),np.max(x)-0.5*step_size,step_size)+0.5*step_size
#     bin_avg = np.zeros(len(bin_centers))
#
#     #We're going to weight with a Gaussian function
#     #def gaussian(x,amp=1,mean=0,sigma=1):
#     #    return amp*np.exp(-(x-mean)**2/(2*sigma**2))
#
#     for index in range(0,len(bin_centers)):
#         bin_center = bin_centers[index]
#         weights = gaussian(x,mean=bin_center,sigma=width)
#         bin_avg[index] = np.average(y,weights=weights)
#
#     return (bin_centers,bin_avg)



mode1Phase3WA = weightedMovingAverage(np.unwrap(phit2S[:,0]*(np.pi/180)),RS[:,0],n=3)
mode2Phase3WA = weightedMovingAverage(np.unwrap(phit2S[:,1]*(np.pi/180)),RS[:,1],n=3)
mode3Phase3WA = weightedMovingAverage(np.unwrap(phit2S[:,2]*(np.pi/180)),RS[:,2],n=3)

mode1Mag3WA = weightedMovingAverage(RS[:,0],RS[:,0],n=3)
mode2Mag3WA = weightedMovingAverage(RS[:,1],RS[:,1],n=3)
mode3Mag3WA = weightedMovingAverage(RS[:,2],RS[:,2],n=3)

mode1Wrap3WA = (mode1Phase3WA + np.pi) % (2 * np.pi) - np.pi
mode2Wrap3WA = (mode2Phase3WA + np.pi) % (2 * np.pi) - np.pi
mode3Wrap3WA = (mode3Phase3WA + np.pi) % (2 * np.pi) - np.pi

mode1Phase5WA = weightedMovingAverage(np.unwrap(phit2S[:,0]*(np.pi/180)),RS[:,0],n=5)
mode2Phase5WA = weightedMovingAverage(np.unwrap(phit2S[:,1]*(np.pi/180)),RS[:,1],n=5)
mode3Phase5WA = weightedMovingAverage(np.unwrap(phit2S[:,2]*(np.pi/180)),RS[:,2],n=5)

mode1Mag5WA = weightedMovingAverage(RS[:,0],RS[:,0],n=5)
mode2Mag5WA = weightedMovingAverage(RS[:,1],RS[:,1],n=5)
mode3Mag5WA = weightedMovingAverage(RS[:,2],RS[:,2],n=5)

mode1Wrap5WA = (mode1Phase5WA + np.pi) % (2 * np.pi) - np.pi
mode2Wrap5WA = (mode2Phase5WA + np.pi) % (2 * np.pi) - np.pi
mode3Wrap5WA = (mode3Phase5WA + np.pi) % (2 * np.pi) - np.pi


mode1Phase3 = moving_average(np.unwrap(phit2S[:,0]*(np.pi/180)),n=3)
mode2Phase3 = moving_average(np.unwrap(phit2S[:,1]*(np.pi/180)),n=3)
mode3Phase3 = moving_average(np.unwrap(phit2S[:,2]*(np.pi/180)),n=3)

mode1Mag3 = moving_average(RS[:,0],n=3)
mode2Mag3 = moving_average(RS[:,1],n=3)
mode3Mag3 = moving_average(RS[:,2],n=3)

mode1Wrap3 = (mode1Phase3 + np.pi) % (2 * np.pi) - np.pi
mode2Wrap3 = (mode2Phase3 + np.pi) % (2 * np.pi) - np.pi
mode3Wrap3 = (mode3Phase3 + np.pi) % (2 * np.pi) - np.pi
days3 = moving_average(days2,n=3)


mode1Phase9 = moving_average(np.unwrap(phit2S[:,0]*(np.pi/180)),n=9)
mode2Phase9 = moving_average(np.unwrap(phit2S[:,1]*(np.pi/180)),n=9)
mode3Phase9 = moving_average(np.unwrap(phit2S[:,2]*(np.pi/180)),n=9)

mode1Mag9 = moving_average(RS[:,0],n=9)
mode2Mag9 = moving_average(RS[:,1],n=9)
mode3Mag9 = moving_average(RS[:,2],n=9)

mode1Wrap9 = (mode1Phase9 + np.pi) % (2 * np.pi) - np.pi
mode2Wrap9 = (mode2Phase9 + np.pi) % (2 * np.pi) - np.pi
mode3Wrap9 = (mode3Phase9 + np.pi) % (2 * np.pi) - np.pi
days9 = moving_average(days2,n=9)

mode1Phase7 = moving_average(np.unwrap(phit2S[:,0]*(np.pi/180)),n=7)
mode2Phase7 = moving_average(np.unwrap(phit2S[:,1]*(np.pi/180)),n=7)
mode3Phase7 = moving_average(np.unwrap(phit2S[:,2]*(np.pi/180)),n=7)

mode1Mag7 = moving_average(RS[:,0],n=7)
mode2Mag7 = moving_average(RS[:,1],n=7)
mode3Mag7 = moving_average(RS[:,2],n=7)

mode1Wrap7 = (mode1Phase7 + np.pi) % (2 * np.pi) - np.pi
mode2Wrap7 = (mode2Phase7 + np.pi) % (2 * np.pi) - np.pi
mode3Wrap7 = (mode3Phase7 + np.pi) % (2 * np.pi) - np.pi
days7 = moving_average(days2,n=7)


mode1Phase5 = moving_average(np.unwrap(phit2S[:,0]*(np.pi/180)),n=5)
mode2Phase5 = moving_average(np.unwrap(phit2S[:,1]*(np.pi/180)),n=5)
mode3Phase5 = moving_average(np.unwrap(phit2S[:,2]*(np.pi/180)),n=5)

mode1Mag5 = moving_average(RS[:,0],n=5)
mode2Mag5 = moving_average(RS[:,1],n=5)
mode3Mag5 = moving_average(RS[:,2],n=5)

mode1Wrap5 = (mode1Phase5 + np.pi) % (2 * np.pi) - np.pi
mode2Wrap5 = (mode2Phase5 + np.pi) % (2 * np.pi) - np.pi
mode3Wrap5 = (mode3Phase5 + np.pi) % (2 * np.pi) - np.pi
days5 = moving_average(days2,n=5)

#
# index = np.where((mode1Wrap5<0))
# mode1Wrap5[index[0]] = mode1Wrap5[index[0]]+np.pi*2
indexShortEnough3 = np.where((days3 < 125))
indexShortEnough5 = np.where((days5 < 125))
indexShortEnough7 = np.where((days7 < 125))
indexShortEnough9 = np.where((days9 < 125))

time5 = timeS[2:-2]
months5 = [i.month for i in time5]
time7 = timeS[3:-3]
months7 = [i.month for i in time7]


seasons = np.nan * np.ones((np.shape(np.array(months7))))
for ff in range(len(months7)):
    if months7[ff] == 12:
        seasons[ff] = 1
    elif months7[ff] >=1 and months7[ff] <= 2:
        seasons[ff] = 1
    elif months7[ff] >=3 and months7[ff] <= 5:
        seasons[ff] = 2
    elif months7[ff] >=6 and months7[ff] <= 8:
        seasons[ff] = 3
    elif months7[ff] >=9 and months7[ff] <= 11:
        seasons[ff] = 4

plt.figure(figsize=(12,12))
ax200 = plt.subplot2grid((3,3),(0,1),rowspan=2,colspan=2)

xs = np.diff(mode1Phase7)[indexShortEnough7]/days7[indexShortEnough7]*180/np.pi
ys = np.diff(mode2Phase7)[indexShortEnough7]/days7[indexShortEnough7]*180/np.pi
zs = seasons[indexShortEnough7]
cols = cm.rainbow(np.linspace(0, 1, 4))
for ff in range(4):
    tempInd = np.where((zs == ff+1))
    if ff == 0:
        scat1 = ax200.scatter(xs[tempInd],ys[tempInd],20,color=cols[ff,:],label='winter',alpha=0.7,edgecolor='none')
    elif ff == 1:
        scat1 = ax200.scatter(xs[tempInd],ys[tempInd],20,color=cols[ff,:],label='spring',alpha=0.7,edgecolor='none')
    elif ff == 2:
        scat1 = ax200.scatter(xs[tempInd],ys[tempInd],20,color=cols[ff,:],label='summer',alpha=0.7,edgecolor='none')
    elif ff == 3:
        scat1 = ax200.scatter(xs[tempInd],ys[tempInd],20,color=cols[ff,:],label='fall',alpha=0.7,edgecolor='none')
plt.legend()
ax200.set_xlim([-1.2,4])
ax200.set_ylim([-1.75, 5.25])
ax200.set_title('Between Survey CEOF Progression')
ax200.text(0,-1.55,r'Offshore Migration (deg/day) $\Longrightarrow$',fontsize=16)
ax200.text(-1.1,0,r'Onshore Migration (deg/day) $\Longrightarrow$',rotation=90,fontsize=16)

ax201 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=1)



bins = np.array((-1.5,-0.75, 0, 0.75, 1.5, 2.25, 3, 3.75, 4.5, 5.25))
ydataSum = np.nan * np.ones((len(bins)-1,))
ydataFall = np.nan * np.ones((len(bins)-1,))
ydataSpr = np.nan * np.ones((len(bins)-1,))
ydataWin = np.nan * np.ones((len(bins)-1,))

for xx in range(len(bins)-1):
    tempo = np.where((ys > bins[xx]) & (ys < bins[xx+1]))
    subset = seasons[tempo]
    for yy in range(4):
        counting = np.where((subset == (yy+1)))
        if yy == 0:
            ydataWin[xx] = len(counting[0])/len(tempo[0])
        elif yy == 1:
            ydataSpr[xx] = len(counting[0])/len(tempo[0])
        elif yy == 2:
            ydataSum[xx] = len(counting[0])/len(tempo[0])
        elif yy == 3:
            ydataFall[xx] = len(counting[0])/len(tempo[0])
ax201.barh(bins[0:-1]+0.375,ydataWin,color=cols[0,:],height=.75)
ax201.barh(bins[0:-1]+0.375,ydataSpr,left=ydataWin,color=cols[1,:],height=.75)
ax201.barh(bins[0:-1]+0.375,ydataSum,left=ydataWin+ydataSpr,color=cols[2,:],height=.75)
ax201.barh(bins[0:-1]+0.375,ydataFall,left=ydataWin+ydataSpr+ydataSum,color=cols[3,:],height=.75)
ax201.set_ylim([-1.75, 5.25])
ax201.set_ylabel(r'Onshore Migration (deg/day) $\Longrightarrow$')
ax201.set_xlabel('Percent of Surveys')

ax202 = plt.subplot2grid((3,3),(2,1),rowspan=1,colspan=2)

bins = np.array((-1, -0.4, 0.2, 0.8, 1.4, 2, 2.6, 3.2, 4))
ydataSum = np.nan * np.ones((len(bins)-1,))
ydataFall = np.nan * np.ones((len(bins)-1,))
ydataSpr = np.nan * np.ones((len(bins)-1,))
ydataWin = np.nan * np.ones((len(bins)-1,))

for xx in range(len(bins)-1):
    tempo = np.where((xs > bins[xx]) & (xs < bins[xx+1]))
    subset = seasons[tempo]
    for yy in range(4):
        counting = np.where((subset == (yy+1)))
        if yy == 0:
            ydataWin[xx] = len(counting[0])/len(tempo[0])
        elif yy == 1:
            ydataSpr[xx] = len(counting[0])/len(tempo[0])
        elif yy == 2:
            ydataSum[xx] = len(counting[0])/len(tempo[0])
        elif yy == 3:
            ydataFall[xx] = len(counting[0])/len(tempo[0])
ax202.bar(bins[0:-1]+0.3,ydataWin,color=cols[0,:],width=.6)
ax202.bar(bins[0:-1]+0.3,ydataSpr,bottom=ydataWin,color=cols[1,:],width=.6)
ax202.bar(bins[0:-1]+0.3,ydataSum,bottom=ydataWin+ydataSpr,color=cols[2,:],width=.6)
ax202.bar(bins[0:-1]+0.3,ydataFall,bottom=ydataWin+ydataSpr+ydataSum,color=cols[3,:],width=.6)
ax202.set_xlim([-1.2,4])
ax202.set_xlabel(r'Offshore Migration (deg/day) $\Longrightarrow$')
ax202.set_ylabel('Percent of Surveys')
plt.tight_layout()
# ax202.hist(x, 5, density=True, histtype='bar', stacked=True)
# ax202.set_title('stacked bar')
plt.show()




#
# index = np.where((mode1Wrap5<0))
# mode1Wrap5[index[0]] = mode1Wrap5[index[0]]+np.pi*2
indexShortEnough3 = np.where((days3 < 125))
indexShortEnough5 = np.where((days5 < 125))
indexShortEnough7 = np.where((days7 < 125))
indexShortEnough9 = np.where((days9 < 125))

time3 = timeS[1:-1]
months3 = [i.month for i in time3]
time5 = timeS[2:-2]
months5 = [i.month for i in time5]
time7 = timeS[3:-3]
months7 = [i.month for i in time7]


seasons = np.nan * np.ones((np.shape(np.array(months5))))
for ff in range(len(months5)):
    if months5[ff] == 12:
        seasons[ff] = 1
    elif months5[ff] >=1 and months5[ff] <= 2:
        seasons[ff] = 1
    elif months5[ff] >=3 and months5[ff] <= 5:
        seasons[ff] = 2
    elif months5[ff] >=6 and months5[ff] <= 8:
        seasons[ff] = 3
    elif months5[ff] >=9 and months5[ff] <= 11:
        seasons[ff] = 4

plt.figure(figsize=(12,12))
ax200 = plt.subplot2grid((3,3),(0,1),rowspan=2,colspan=2)

xs = np.diff(mode1Phase5)[indexShortEnough5]/days5[indexShortEnough5]*180/np.pi
ys = np.diff(mode2Phase5)[indexShortEnough5]/days5[indexShortEnough5]*180/np.pi
zs = seasons[indexShortEnough5]
cols = cm.rainbow(np.linspace(0, 1, 4))
for ff in range(4):
    tempInd = np.where((zs == ff+1))
    if ff == 0:
        scat1 = ax200.scatter(xs[tempInd],ys[tempInd],20,color=cols[ff,:],label='winter',alpha=0.7,edgecolor='none')
    elif ff == 1:
        scat1 = ax200.scatter(xs[tempInd],ys[tempInd],20,color=cols[ff,:],label='spring',alpha=0.7,edgecolor='none')
    elif ff == 2:
        scat1 = ax200.scatter(xs[tempInd],ys[tempInd],20,color=cols[ff,:],label='summer',alpha=0.7,edgecolor='none')
    elif ff == 3:
        scat1 = ax200.scatter(xs[tempInd],ys[tempInd],20,color=cols[ff,:],label='fall',alpha=0.7,edgecolor='none')
plt.legend()
ax200.set_xlim([-1.2,5])
ax200.set_ylim([-2, 6])
ax200.set_title('Between Survey CEOF Progression')
ax200.text(0,-1.75,r'Offshore Migration (deg/day) $\Longrightarrow$',fontsize=16)
ax200.text(-1.1,0,r'Onshore Migration (deg/day) $\Longrightarrow$',rotation=90,fontsize=16)

ax201 = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=1)

bins = np.array((-2,-1, 0, 1, 2, 3, 4, 5, 6))
ydataSum = np.nan * np.ones((len(bins)-1,))
ydataFall = np.nan * np.ones((len(bins)-1,))
ydataSpr = np.nan * np.ones((len(bins)-1,))
ydataWin = np.nan * np.ones((len(bins)-1,))

for xx in range(len(bins)-1):
    tempo = np.where((ys > bins[xx]) & (ys < bins[xx+1]))
    subset = seasons[tempo]
    for yy in range(4):
        counting = np.where((subset == (yy+1)))
        if yy == 0:
            ydataWin[xx] = len(counting[0])/len(tempo[0])
        elif yy == 1:
            ydataSpr[xx] = len(counting[0])/len(tempo[0])
        elif yy == 2:
            ydataSum[xx] = len(counting[0])/len(tempo[0])
        elif yy == 3:
            ydataFall[xx] = len(counting[0])/len(tempo[0])
ax201.barh(bins[0:-1]+0.5,ydataWin,color=cols[0,:],height=1)
ax201.barh(bins[0:-1]+0.5,ydataSpr,left=ydataWin,color=cols[1,:],height=1)
ax201.barh(bins[0:-1]+0.5,ydataSum,left=ydataWin+ydataSpr,color=cols[2,:],height=1)
ax201.barh(bins[0:-1]+0.5,ydataFall,left=ydataWin+ydataSpr+ydataSum,color=cols[3,:],height=1)
ax201.set_ylim([-2, 6])
ax201.set_ylabel(r'Onshore Migration (deg/day) $\Longrightarrow$')
ax201.set_xlabel('Percent of Surveys')

ax202 = plt.subplot2grid((3,3),(2,1),rowspan=1,colspan=2)

bins = np.array((-1, -0, 1, 2, 3, 4, 5))
ydataSum = np.nan * np.ones((len(bins)-1,))
ydataFall = np.nan * np.ones((len(bins)-1,))
ydataSpr = np.nan * np.ones((len(bins)-1,))
ydataWin = np.nan * np.ones((len(bins)-1,))

for xx in range(len(bins)-1):
    tempo = np.where((xs > bins[xx]) & (xs < bins[xx+1]))
    subset = seasons[tempo]
    for yy in range(4):
        counting = np.where((subset == (yy+1)))
        if yy == 0:
            ydataWin[xx] = len(counting[0])/len(tempo[0])
        elif yy == 1:
            ydataSpr[xx] = len(counting[0])/len(tempo[0])
        elif yy == 2:
            ydataSum[xx] = len(counting[0])/len(tempo[0])
        elif yy == 3:
            ydataFall[xx] = len(counting[0])/len(tempo[0])
ax202.bar(bins[0:-1]+0.5,ydataWin,color=cols[0,:],width=1)
ax202.bar(bins[0:-1]+0.5,ydataSpr,bottom=ydataWin,color=cols[1,:],width=1)
ax202.bar(bins[0:-1]+0.5,ydataSum,bottom=ydataWin+ydataSpr,color=cols[2,:],width=1)
ax202.bar(bins[0:-1]+0.5,ydataFall,bottom=ydataWin+ydataSpr+ydataSum,color=cols[3,:],width=1)
ax202.set_xlim([-1.2,5])
ax202.set_xlabel(r'Offshore Migration (deg/day) $\Longrightarrow$')
ax202.set_ylabel('Percent of Surveys')
plt.tight_layout()
# ax202.hist(x, 5, density=True, histtype='bar', stacked=True)
# ax202.set_title('stacked bar')
plt.show()




bigShifts = np.where((xs > 3))



#
#
#
# fig = plt.figure(figsize=(10,5))
# ax1 = plt.subplot2grid((3,4),(0,0),rowspan=1,colspan=2)
# xg,tg = np.meshgrid(xinterpS,timeS)
# plt0 = ax1.pcolor(tg,xg,(alllinesS-np.mean(alllinesS, axis=0)), vmin=-1.8, vmax=1.8,cmap='bwr')
#
# ax = plt.subplot2grid((3,4),(1,0),rowspan=1,colspan=2,sharex=ax1)
# temp = mode1Mag*15
# temp1 = mode1Wrap
# # index = np.where((temp1<0))
# # temp1[index[0]] = temp1[index[0]]+360
# ax.scatter(timeS[1:-1],temp1,temp,mode1Mag,cmap='Reds')
#
# ax2 = plt.subplot2grid((3,4),(2,0),rowspan=1,colspan=2, sharex=ax1)
# temp = mode2Mag*10
# temp1 = mode2Wrap
# ax2.scatter(timeS[1:-1],temp1,temp,mode2Mag,cmap='Blues')
# ax1.set_xlim([timeS[0], timeS[110]])
# # ax1.set_xlim([timeS[0], timeS[-1]])
# # ax2.set_xlim([timeS[0], timeS[-1]])
#
# ax1b = plt.subplot2grid((3,4),(0,2),rowspan=1,colspan=2)
# xg,tg = np.meshgrid(xinterpS,timeS)
# plt0b = ax1b.pcolor(tg,xg,(alllinesS-np.mean(alllinesS, axis=0)), vmin=-1.8, vmax=1.8,cmap='bwr')
# axb = plt.subplot2grid((3,4),(1,2),rowspan=1,colspan=2,sharex=ax1b)
# temp = mode1Mag*15
# temp1 = mode1Wrap
# # index = np.where((temp1<0))
# # temp1[index[0]] = temp1[index[0]]+360
# axb.scatter(timeS[1:-1],temp1,temp,mode1Mag,cmap='Reds')
#
# ax2b = plt.subplot2grid((3,4),(2,2),rowspan=1,colspan=2, sharex=ax1b)
# temp = mode2Mag*10
# temp1 = mode2Wrap
# ax2b.scatter(timeS[1:-1],temp1,temp,mode2Mag,cmap='Blues')
# ax1b.set_xlim([timeS[-80], timeS[-1]])
#
#
# plt.show()







fig = plt.figure(figsize=(8,8))
ax1 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
# p1 = ax1.scatter(phit2S[:,0],phit2S[:,1],RS[:,1]*100,RS[:,0],vmin=0,vmax=1.75,cmap='Reds')
# p1 = ax1.scatter(mode1Wrap,mode2Wrap,mode1Mag*50,mode2Mag,vmin=0,vmax=1.5,cmap='Reds',alpha=0.5,edgecolor='none')
p1 = ax1.scatter(mode1Wrap,mode2Wrap,mode2Mag*150,'r',edgecolor='none')

ax1.set_xlabel('Mode 1 Phase')
ax1.set_ylabel('Mode 2 Phase')
fig.colorbar(p1,ax=ax1)



import datetime as DT





index1 = np.where((time5 > DT.datetime(1984,11,1)) & (time5 < DT.datetime(1987,5,1)))
indexedTime = time5[index1]
indexedOrdinal = [i.toordinal() for i in indexedTime]
indexedTime1Month = [i.month for i in indexedTime]
mode1prop = mode1Phase5[index1[0]]-6*np.pi#np.unwrap(phitS[index1[0],0])
mode2prop = mode2Phase5[index1[0]]-6*np.pi#np.unwrap(phitS[index1[0],1])

# index2 = np.where((time7 > DT.datetime(1981,8,1)) & (time7 < DT.datetime(1987,9,1)))
# indexedTime2 = time7[index2]
# indexedOrdinal2 = [i.toordinal() for i in indexedTime2]
# mode1prop2 = mode1Phase7[index2[0]]-2*np.pi#np.unwrap(phitS[index1[0],0])
# mode2prop2 = mode2Phase7[index2[0]]-2*np.pi#np.unwrap(phitS[index1[0],1])

index2 = np.where((time5 > DT.datetime(1991,9,1)) & (time5 < DT.datetime(1995,10,1)))
indexedTime2 = time5[index2]
indexedOrdinal2 = [i.toordinal() for i in indexedTime2]
indexedTime2Month = [i.month for i in indexedTime2]
mode1prop2 = mode1Phase5[index2[0]]-14*np.pi#np.unwrap(phitS[index1[0],0])
mode2prop2 = mode2Phase5[index2[0]]-14*np.pi#np.unwrap(phitS[index1[0],1])

index3 = np.where((time5 > DT.datetime(1995,9,1)) & (time5 < DT.datetime(1999,7,1)))
indexedTime3 = time5[index3]
indexedOrdinal3 = [i.toordinal() for i in indexedTime3]
indexedTime3Month = [i.month for i in indexedTime3]
mode1prop3 = mode1Phase5[index3[0]]-18*np.pi#np.unwrap(phitS[index1[0],0])
mode2prop3 = mode2Phase5[index3[0]]-18*np.pi#np.unwrap(phitS[index1[0],1])

index4 = np.where((time5 > DT.datetime(2012,9,1)) & (time5 < DT.datetime(2016,11,1)))
indexedTime4 = time5[index4]
indexedOrdinal4 = [i.toordinal() for i in indexedTime4]
indexedTime4Month = [i.month for i in indexedTime4]
mode1prop4 = mode1Phase5[index4[0]]-38*np.pi#np.unwrap(phitS[index1[0],0])
mode2prop4 = mode2Phase5[index4[0]]-42*np.pi#np.unwrap(phitS[index1[0],1])

# index4 = np.where((time5 > DT.datetime(2016,10,1)) & (time5 < DT.datetime(2019,11,1)))
# indexedTime4 = time5[index4]
# indexedOrdinal4 = [i.toordinal() for i in indexedTime4]
# indexedTime4Month = [i.month for i in indexedTime4]
# mode1prop4 = mode1Phase5[index4[0]]-38*np.pi#np.unwrap(phitS[index1[0],0])
# mode2prop4 = mode2Phase5[index4[0]]-42*np.pi#np.unwrap(phitS[index1[0],1])

# plt.figure()
# plt.plot(mode1prop*180/np.pi,mode2prop*180/np.pi,'-o')
# plt.plot(mode1prop2*180/np.pi,mode2prop2*180/np.pi,'o-')
# plt.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'o-')
#
# plt.show()




fig2 = plt.figure(figsize=(14,8))
ax1 = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
p1 = ax1.pcolor(meshTheta1,meshTheta2,innerBarX,vmin=30,vmax=300)
plt.plot(mode1prop*180/np.pi,mode2prop*180/np.pi,'-o')
plt.plot(mode1prop2*180/np.pi,mode2prop2*180/np.pi,'o-')
plt.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'o-')
plt.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'o-')

cb1 = plt.colorbar(p1,ax=ax1)
cb1.set_label('xFRF')
ax2 = plt.subplot2grid((2,3),(0,1),rowspan=1,colspan=1)
p2 = ax2.pcolor(meshTheta1,meshTheta2,outerBarX,vmin=30,vmax=300)
plt.plot(mode1prop*180/np.pi,mode2prop*180/np.pi,'-o')
plt.plot(mode1prop2*180/np.pi,mode2prop2*180/np.pi,'o-')
plt.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'o-')
plt.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'o-')

ax2.set_xlabel('Mode 1 Phase')
ax1.set_xlabel('Mode 1 Phase')
# ax3.set_xlabel('Mode 1 Phase')
ax1.set_ylabel('Mode 2 Phase')
ax1.set_title('Inner Bar Cross-shore Location')
ax2.set_title('Outer Bar Cross-shore Location')

cb2 = plt.colorbar(p2,ax=ax2)
cb2.set_label('xFRF')
ax1c = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
p1c = ax1c.pcolor(meshTheta1,meshTheta2,innerOuterBarX,vmin=30,vmax=300)
plt.plot(mode1prop*180/np.pi,mode2prop*180/np.pi,'-o')
plt.plot(mode1prop2*180/np.pi,mode2prop2*180/np.pi,'o-')
plt.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'o-')
plt.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'o-')


cb1c = plt.colorbar(p1c,ax=ax1c)
cb1c.set_label('xFRF')
ax1c.set_title('Single Bar Cross-shore Location')

ax3 = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
p3 = ax3.pcolor(meshTheta1,meshTheta2,innerBarDZ,vmin=0,vmax=1)
plt.plot(mode1prop*180/np.pi,mode2prop*180/np.pi,'-o')
plt.plot(mode1prop2*180/np.pi,mode2prop2*180/np.pi,'o-')
plt.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'o-')
plt.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'o-')

cb3 = plt.colorbar(p3,ax=ax3)
cb3.set_label('Depth (m)')

ax4 = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
p4 = ax4.pcolor(meshTheta1,meshTheta2,outerBarDZ,vmin=0,vmax=1)
plt.plot(mode1prop*180/np.pi,mode2prop*180/np.pi,'-o')
plt.plot(mode1prop2*180/np.pi,mode2prop2*180/np.pi,'o-')
plt.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'o-')
plt.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'o-')

cb4 = plt.colorbar(p4,ax=ax4)
cb4.set_label('Depth (m)')

ax4c = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
p4c = ax4c.pcolor(meshTheta1,meshTheta2,innerOuterBarDZ,vmin=0,vmax=1)
plt.plot(mode1prop*180/np.pi,mode2prop*180/np.pi,'-o')
plt.plot(mode1prop2*180/np.pi,mode2prop2*180/np.pi,'o-')
plt.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'o-')
plt.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'o-')

cb4c = plt.colorbar(p4c,ax=ax4c)
cb4c.set_label('Depth (m)')

ax4c.set_xlabel('Mode 1 Phase')
ax4.set_xlabel('Mode 1 Phase')
ax3.set_xlabel('Mode 1 Phase')
ax3.set_ylabel('Mode 2 Phase')
ax3.set_title('Inner Bar Magnitude')
ax4.set_title('Outer Bar Magnitude')
ax4c.set_title('Single Bar Magnitude')

plt.tight_layout()
plt.show()






mode1PhaseWeighted5 = weightedMovingAverage(np.unwrap(phit2S[:,0]*(np.pi/180)),RS[:,0],5)
mode2PhaseWeighted5 = weightedMovingAverage(np.unwrap(phit2S[:,1]*(np.pi/180)),RS[:,1],5)
# mode1PhaseWeighted5 = weightedMovingAverage(np.unwrap(phit2S[:,0]*(np.pi/180)),RS[:,0],3)
# mode2PhaseWeighted5 = weightedMovingAverage(np.unwrap(phit2S[:,1]*(np.pi/180)),RS[:,1],3)

mode1PhaseWeighted5 = np.unwrap(phit2S[:,0]*(np.pi/180))
mode2PhaseWeighted5 = np.unwrap(phit2S[:,1]*(np.pi/180))


# index1 = np.where((time3 > DT.datetime(1984,11,1)) & (time3 < DT.datetime(1987,5,1)))
index1 = np.where((time3 > DT.datetime(1984,10,25)) & (time3 < DT.datetime(1988,2,20)))

indexedTime = time3[index1]
indexedOrdinal = [i.toordinal() for i in indexedTime]
indexedTime1Month = [i.month for i in indexedTime]
mode1prop = mode1PhaseWeighted5[index1[0]]-6*np.pi+2*np.pi#np.unwrap(phitS[index1[0],0])
mode2prop = mode2PhaseWeighted5[index1[0]]-6*np.pi#np.unwrap(phitS[index1[0],1])

# index2 = np.where((time7 > DT.datetime(1981,8,1)) & (time7 < DT.datetime(1987,9,1)))
# indexedTime2 = time7[index2]
# indexedOrdinal2 = [i.toordinal() for i in indexedTime2]
# mode1prop2 = mode1Phase7[index2[0]]-2*np.pi#np.unwrap(phitS[index1[0],0])
# mode2prop2 = mode2Phase7[index2[0]]-2*np.pi#np.unwrap(phitS[index1[0],1])

index2 = np.where((time3 > DT.datetime(1991,9,15)) & (time3 < DT.datetime(1995,10,1)))
indexedTime2 = time3[index2]
indexedOrdinal2 = [i.toordinal() for i in indexedTime2]
indexedTime2Month = [i.month for i in indexedTime2]
mode1prop2 = mode1PhaseWeighted5[index2[0]]-14*np.pi+2*np.pi#np.unwrap(phitS[index1[0],0])
mode2prop2 = mode2PhaseWeighted5[index2[0]]-14*np.pi#np.unwrap(phitS[index1[0],1])

# index3 = np.where((time3 > DT.datetime(1995,8,1)) & (time3 < DT.datetime(1999,7,1)))
index3 = np.where((time3 > DT.datetime(1995,11,1)) & (time3 < DT.datetime(2002,3,10)))

indexedTime3 = time3[index3]
indexedOrdinal3 = [i.toordinal() for i in indexedTime3]
indexedTime3Month = [i.month for i in indexedTime3]
mode1prop3 = mode1PhaseWeighted5[index3[0]]-18*np.pi+2*np.pi#np.unwrap(phitS[index1[0],0])
mode2prop3 = mode2PhaseWeighted5[index3[0]]-18*np.pi#np.unwrap(phitS[index1[0],1])

index4 = np.where((time3 > DT.datetime(2012,9,1)) & (time3 < DT.datetime(2016,11,1)))
indexedTime4 = time3[index4]
indexedOrdinal4 = [i.toordinal() for i in indexedTime4]
indexedTime4Month = [i.month for i in indexedTime4]
mode1prop4 = mode1PhaseWeighted5[index4[0]]-38*np.pi+2*np.pi#np.unwrap(phitS[index1[0],0])
mode2prop4 = mode2PhaseWeighted5[index4[0]]-42*np.pi#np.unwrap(phitS[index1[0],1])



mode1propDeg = mode1prop*180/np.pi
mode1propDegAlt = np.ones((len(mode1propDeg),))
mode2propDeg = mode2prop*180/np.pi
mode2propDegAlt = np.ones((len(mode1propDeg),))

index1RS = RS[index1, 0]
mode1propDegAlt[0] = mode1propDeg[0]
for ii in range(len(mode1propDeg)):
    if np.real(index1RS[0][ii]) > 0.55:
        mode1propDegAlt[ii] = mode1propDeg[ii]
    else:
        print('found one, initial degree {}'.format(mode1propDeg[ii]))
        mode1propDeg[ii:] = mode1propDeg[ii:] - np.subtract(mode1propDeg[ii],mode1propDeg[ii-1])
        mode1propDegAlt[ii] = mode1propDeg[ii]
        print('processed degree {}'.format(mode1propDeg[ii]))



mode1propDeg3 = mode1prop3*180/np.pi
mode1propDegAlt3 = np.ones((len(mode1propDeg3),))
mode2propDeg3 = mode2prop3*180/np.pi
mode2propDegAlt3 = np.ones((len(mode1propDeg3),))

index3RS = RS[index3, 0]
mode1propDegAlt3[0] = mode1propDeg3[0]
for ii in range(len(mode1propDeg3)):
    if np.real(index3RS[0][ii]) > 0.55:
        mode1propDegAlt3[ii] = mode1propDeg3[ii]
    else:
        print('found one, initial degree {}'.format(mode1propDeg3[ii]))
        mode1propDeg3[ii:] = mode1propDeg3[ii:] - np.subtract(mode1propDeg3[ii],mode1propDeg3[ii-1])
        mode1propDegAlt3[ii] = mode1propDeg3[ii]
        print('processed degree {}'.format(mode1propDeg3[ii]))


mode1propDeg2 = mode1prop2*180/np.pi
mode1propDegAlt2 = np.ones((len(mode1propDeg2),))
mode2propDeg2 = mode2prop2*180/np.pi
mode2propDegAlt2 = np.ones((len(mode1propDeg2),))

index2RS = RS[index2, 0]
mode1propDegAlt2[0] = mode1propDeg2[0]
for ii in range(len(mode1propDeg2)):
    if np.real(index2RS[0][ii]) > 0.55:
        mode1propDegAlt2[ii] = mode1propDeg2[ii]
    else:
        print('found one, initial degree {}'.format(mode1propDeg2[ii]))
        mode1propDeg2[ii:] = mode1propDeg2[ii:] - np.subtract(mode1propDeg2[ii],mode1propDeg2[ii-1])
        mode1propDegAlt2[ii] = mode1propDeg2[ii]
        print('processed degree {}'.format(mode1propDeg2[ii]))




fig = plt.figure(figsize=(12,12))
ax100 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)

con1 = ax100.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
lines1 = ax100.plot(mode1prop*180/np.pi,mode2prop*180/np.pi,'-',label='2.5 years')
# # scatter1 = ax100.scatter(mode1prop*180/np.pi,mode2prop*180/np.pi,35,indexedTime1Month,cmap='twilight_shifted',zorder=3)
# scatter1 = ax100.scatter(mode1prop*180/np.pi,mode2prop*180/np.pi,35,RS[index1,0],cmap='bwr',zorder=3)
# lines1 = ax100.plot(mode1propDegAlt,mode2prop*180/np.pi,'-',label='2.5 years')
scatter1 = ax100.scatter(mode1prop*180/np.pi,mode2prop*180/np.pi,35,indexedTime1Month,cmap='twilight_shifted',zorder=3)
# scatter1 = ax100.scatter(mode1propDegAlt,mode2prop*180/np.pi,35,RS[index1,0],cmap='bwr',zorder=3)


lines2 = ax100.plot(mode1prop2*180/np.pi,mode2prop2*180/np.pi,'-',label='4.1 years')
scatter2 = ax100.scatter(mode1prop2*180/np.pi,mode2prop2*180/np.pi,35,indexedTime2Month,cmap='twilight_shifted',zorder=3)
# lines2 = ax100.plot(mode1propDegAlt2,mode2prop2*180/np.pi,'-',label='4.1 years')
# scatter1 = ax100.scatter(mode1propDegAlt2,mode2prop2*180/np.pi,35,RS[index2,0],cmap='bwr',zorder=3,vmin=0.25)

# lines3 = ax100.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'-',label='3.8 years')
# # scatter3 = ax100.scatter(mode1prop3*180/np.pi,mode2prop3*180/np.pi,35,indexedTime3Month,cmap='twilight_shifted',zorder=3)
# scatter1 = ax100.scatter(mode1prop3*180/np.pi,mode2prop3*180/np.pi,35,RS[index3,0],cmap='bwr',zorder=3)

lines3 = ax100.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'-',label='3.8 years')
scatter3 = ax100.scatter(mode1prop3*180/np.pi,mode2prop3*180/np.pi,35,indexedTime3Month,cmap='twilight_shifted',zorder=3)
# lines3 = ax100.plot(mode1propDegAlt3,mode2prop3*180/np.pi,'-',label='3.8 years')
# scatter1 = ax100.scatter(mode1propDegAlt3,mode2prop3*180/np.pi,35,RS[index3,0],cmap='bwr',zorder=3)


lines4 = ax100.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'-',label='4.2 years')
scatter4 = ax100.scatter(mode1prop4*180/np.pi,mode2prop4*180/np.pi,35,indexedTime4Month,cmap='twilight_shifted',zorder=3)
# scatter1 = ax100.scatter(mode1prop4*180/np.pi,mode2prop4*180/np.pi,35,RS[index4,0],cmap='bwr',zorder=3)

# plt.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'o-')
# plt.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'o-')
plt.legend()
ax100.set_xlim([-300, 1000])
ax100.set_ylim([-290, 1600])
cb = fig.colorbar(scatter1,ax=ax100)
plt.show()






fig = plt.figure(figsize=(12,12))
ax100 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)

con1 = ax100.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
# lines1 = ax100.plot(mode1prop*180/np.pi,mode2prop*180/np.pi,'-',label='2.5 years')
# # scatter1 = ax100.scatter(mode1prop*180/np.pi,mode2prop*180/np.pi,35,indexedTime1Month,cmap='twilight_shifted',zorder=3)
# scatter1 = ax100.scatter(mode1prop*180/np.pi,mode2prop*180/np.pi,35,RS[index1,0],cmap='bwr',zorder=3)
lines1 = ax100.plot(mode1propDegAlt,mode2prop*180/np.pi,'-',label='2.5 years')
# scatter1 = ax100.scatter(mode1prop*180/np.pi,mode2prop*180/np.pi,35,indexedTime1Month,cmap='twilight_shifted',zorder=3)
scatter1 = ax100.scatter(mode1propDegAlt,mode2prop*180/np.pi,35,RS[index1,0],cmap='bwr',zorder=3)


lines2 = ax100.plot(mode1propDegAlt2,mode2prop2*180/np.pi,'-',label='4.1 years')
# scatter2 = ax100.scatter(mode1prop2*180/np.pi,mode2prop2*180/np.pi,35,indexedTime2Month,cmap='twilight_shifted',zorder=3)
scatter1 = ax100.scatter(mode1propDegAlt2,mode2prop2*180/np.pi,35,RS[index2,0],cmap='bwr',zorder=3,vmin=0.25)

# lines3 = ax100.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'-',label='3.8 years')
# # scatter3 = ax100.scatter(mode1prop3*180/np.pi,mode2prop3*180/np.pi,35,indexedTime3Month,cmap='twilight_shifted',zorder=3)
# scatter1 = ax100.scatter(mode1prop3*180/np.pi,mode2prop3*180/np.pi,35,RS[index3,0],cmap='bwr',zorder=3)

lines3 = ax100.plot(mode1propDegAlt3,mode2prop3*180/np.pi,'-',label='3.8 years')
# scatter3 = ax100.scatter(mode1prop3*180/np.pi,mode2prop3*180/np.pi,35,indexedTime3Month,cmap='twilight_shifted',zorder=3)
scatter1 = ax100.scatter(mode1propDegAlt3,mode2prop3*180/np.pi,35,RS[index3,0],cmap='bwr',zorder=3)


lines4 = ax100.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'-',label='4.2 years')
# scatter4 = ax100.scatter(mode1prop4*180/np.pi,mode2prop4*180/np.pi,35,indexedTime4Month,cmap='twilight_shifted',zorder=3)
scatter1 = ax100.scatter(mode1prop4*180/np.pi,mode2prop4*180/np.pi,35,RS[index4,0],cmap='bwr',zorder=3)

# plt.plot(mode1prop3*180/np.pi,mode2prop3*180/np.pi,'o-')
# plt.plot(mode1prop4*180/np.pi,mode2prop4*180/np.pi,'o-')
plt.legend()
ax100.set_xlim([-300, 1000])
ax100.set_ylim([-290, 1600])
cb = fig.colorbar(scatter1,ax=ax100)
plt.show()





fig = plt.figure(figsize=(12,12))
ax100 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
con1 = ax100.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')

ax100.set_xlim([-270, 630])
ax100.set_ylim([-270, 630])
plt.show()







fig2 = plt.figure(figsize=(16,12))


# index1 = np.where((timeS > DT.datetime(1981,7,1)) & (timeS < DT.datetime(1984,9,1)))
index1 = np.where((timeS > DT.datetime(2016,11,1)) & (timeS < DT.datetime(2019,11,1)))

newTime = [date - timeS[index1[0][0]] for date in timeS[index1[0]]]
yearTime = [daysPassed.days/365 for daysPassed in newTime]
phiTemp = phit2S[index1,mode]
phiTemp[0][15:] = phiTemp[0][15:]+360
subsetOfAllines = alllines[index1[0],:]
ax1 = plt.subplot2grid((3,4),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.Reds(np.linspace(0, 1, len(subsetOfAllines)+5))
subsetAgain = np.arange(0,15)+5
for i in subsetAgain:
    ax1.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax2 = plt.subplot2grid((3,4),(0,1),rowspan=1,colspan=1)
subsetAgain = np.arange(16,38)+5
for i in subsetAgain:
    ax2.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax3 = plt.subplot2grid((3,4),(0,2),rowspan=2,colspan=2)

ax3.scatter(yearTime,phiTemp,RS[index1,mode]*25,colors[4:-1,:],label='2016-2019')



index1 = np.where((timeS > DT.datetime(1995,8,15)) & (timeS < DT.datetime(1998,1,20)))

newTime = [date - timeS[index1[0][0]] for date in timeS[index1[0]]]
yearTime = [daysPassed.days/365 for daysPassed in newTime]
phiTemp = phit2S[index1,mode]
phiTemp[0][7:] = phiTemp[0][7:]+360
subsetOfAllines = alllines[index1[0],:]
ax1b = plt.subplot2grid((3,4),(1,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.Blues(np.linspace(0, 1, len(subsetOfAllines)+5))
subsetAgain = np.arange(0,12)+5
for i in subsetAgain:
    ax1b.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax2b = plt.subplot2grid((3,4),(1,1),rowspan=1,colspan=1)
subsetAgain = np.arange(12,39)+5
for i in subsetAgain:
    ax2b.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

#ax3b = plt.subplot2grid((3,4),(1,2),rowspan=1,colspan=1)
ax3.scatter(yearTime,phiTemp,RS[index1,mode]*25,colors[4:-1],label='1995-1998')




# index1 = np.where((timeS > DT.datetime(2011,11,15)) & (timeS < DT.datetime(2016,8,10)))
index1 = np.where((timeS > DT.datetime(2012,6,15)) & (timeS < DT.datetime(2016,8,10)))

newTime = [date - timeS[index1[0][0]] for date in timeS[index1[0]]]
yearTime = [daysPassed.days/365 for daysPassed in newTime]
phiTemp = phit2S[index1,mode]
phiTemp[0][5:] = phiTemp[0][5:]+360
subsetOfAllines = alllines[index1[0],:]
ax1c = plt.subplot2grid((3,4),(2,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.Greens(np.linspace(0, 1, len(subsetOfAllines)+5))
subsetAgain = np.arange(0,16)+5
for i in subsetAgain:
    ax1c.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax2c = plt.subplot2grid((3,4),(2,1),rowspan=1,colspan=1)
subsetAgain = np.arange(16,41)+5
for i in subsetAgain:
    ax2c.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

#ax3c = plt.subplot2grid((4,3),(2,2),rowspan=1,colspan=1)
ax3.scatter(yearTime,phiTemp,RS[index1,mode]*25,colors[4:-1],label='2012-2016')
# ax3.set_ylim([60,490])
# ax3.set_xlim([0,4.5])


# index1 = np.where((timeS > DT.datetime(2011,11,15)) & (timeS < DT.datetime(2016,8,10)))
index1 = np.where((timeS > DT.datetime(1988,6,15)) & (timeS < DT.datetime(1991,6,10)))

newTime = [date - timeS[index1[0][0]] for date in timeS[index1[0]]]
yearTime = [daysPassed.days/365 for daysPassed in newTime]
phiTemp = phit2S[index1,mode]
phiTemp[0][9:] = phiTemp[0][9:]+360
subsetOfAllines = alllines[index1[0],:]
ax1c = plt.subplot2grid((3,4),(2,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.Purples(np.linspace(0, 1, len(subsetOfAllines)+5))
subsetAgain = np.arange(0,16)+5
for i in subsetAgain:
    ax1c.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax2c = plt.subplot2grid((3,4),(2,3),rowspan=1,colspan=1)
subsetAgain = np.arange(16,56)+5
for i in subsetAgain:
    ax2c.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

#ax3c = plt.subplot2grid((4,3),(2,2),rowspan=1,colspan=1)
ax3.scatter(yearTime,phiTemp,RS[index1,mode]*25,colors[4:-1],label='1988-1991')
# ax3.set_ylim([60,490])
# ax3.set_xlim([0,4.5])
ax3.legend()










mode = 0


fig3 = plt.figure(figsize=(16,12))


index1 = np.where((timeS > DT.datetime(1981,7,1)) & (timeS < DT.datetime(1984,9,1)))
# index1 = np.where((timeS > DT.datetime(1993,8,1)) & (timeS < DT.datetime(1996,8,1)))

newTime = [date - timeS[index1[0][0]] for date in timeS[index1[0]]]
yearTime = [daysPassed.days/365 for daysPassed in newTime]
phiTemp = phit2S[index1,mode]
phiTemp[0][15:] = phiTemp[0][15:]+360
subsetOfAllines = alllines[index1[0],:]
ax1 = plt.subplot2grid((3,4),(0,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.Reds(np.linspace(0, 1, len(subsetOfAllines)+5))
subsetAgain = np.arange(0,15)+5
for i in subsetAgain:
    ax1.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax2 = plt.subplot2grid((3,4),(0,1),rowspan=1,colspan=1)
subsetAgain = np.arange(16,53)+5
for i in subsetAgain:
    ax2.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax3 = plt.subplot2grid((3,4),(0,2),rowspan=2,colspan=2)

ax3.scatter(yearTime,phiTemp,RS[index1,mode]*25,colors[4:-1,:])#,label='2016-2019')



index1 = np.where((timeS > DT.datetime(1995,8,15)) & (timeS < DT.datetime(1998,1,20)))

newTime = [date - timeS[index1[0][0]] for date in timeS[index1[0]]]
yearTime = [daysPassed.days/365 for daysPassed in newTime]
phiTemp = phit2S[index1,mode]
phiTemp[0][7:] = phiTemp[0][7:]+360
subsetOfAllines = alllines[index1[0],:]
ax1b = plt.subplot2grid((3,4),(1,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.Blues(np.linspace(0, 1, len(subsetOfAllines)+5))
subsetAgain = np.arange(0,12)+5
for i in subsetAgain:
    ax1b.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax2b = plt.subplot2grid((3,4),(1,1),rowspan=1,colspan=1)
subsetAgain = np.arange(12,39)+5
for i in subsetAgain:
    ax2b.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

#ax3b = plt.subplot2grid((3,4),(1,2),rowspan=1,colspan=1)
ax3.scatter(yearTime,phiTemp,RS[index1,mode]*25,colors[4:-1],label='1995-1998')




# index1 = np.where((timeS > DT.datetime(2011,11,15)) & (timeS < DT.datetime(2016,8,10)))
index1 = np.where((timeS > DT.datetime(2012,6,15)) & (timeS < DT.datetime(2016,8,10)))

newTime = [date - timeS[index1[0][0]] for date in timeS[index1[0]]]
yearTime = [daysPassed.days/365 for daysPassed in newTime]
phiTemp = phit2S[index1,mode]
phiTemp[0][5:] = phiTemp[0][5:]+360
subsetOfAllines = alllines[index1[0],:]
ax1c = plt.subplot2grid((3,4),(2,0),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.Greens(np.linspace(0, 1, len(subsetOfAllines)+5))
subsetAgain = np.arange(0,16)+5
for i in subsetAgain:
    ax1c.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax2c = plt.subplot2grid((3,4),(2,1),rowspan=1,colspan=1)
subsetAgain = np.arange(16,41)+5
for i in subsetAgain:
    ax2c.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

#ax3c = plt.subplot2grid((4,3),(2,2),rowspan=1,colspan=1)
ax3.scatter(yearTime,phiTemp,RS[index1,mode]*25,colors[4:-1],label='2012-2016')
# ax3.set_ylim([60,490])
# ax3.set_xlim([0,4.5])


# index1 = np.where((timeS > DT.datetime(2011,11,15)) & (timeS < DT.datetime(2016,8,10)))
index1 = np.where((timeS > DT.datetime(1988,6,15)) & (timeS < DT.datetime(1991,6,10)))

newTime = [date - timeS[index1[0][0]] for date in timeS[index1[0]]]
yearTime = [daysPassed.days/365 for daysPassed in newTime]
phiTemp = phit2S[index1,mode]
phiTemp[0][9:] = phiTemp[0][9:]+360
subsetOfAllines = alllines[index1[0],:]
ax1c = plt.subplot2grid((3,4),(2,2),rowspan=1,colspan=1)
import matplotlib.cm as cm
colors = cm.Purples(np.linspace(0, 1, len(subsetOfAllines)+5))
subsetAgain = np.arange(0,16)+5
for i in subsetAgain:
    ax1c.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

ax2c = plt.subplot2grid((3,4),(2,3),rowspan=1,colspan=1)
subsetAgain = np.arange(16,56)+5
for i in subsetAgain:
    ax2c.plot(xinterpS, subsetOfAllines[i-5,:],color=colors[i,:])#, label=time[i])

#ax3c = plt.subplot2grid((4,3),(2,2),rowspan=1,colspan=1)
ax3.scatter(yearTime,phiTemp,RS[index1,mode]*25,colors[4:-1],label='1988-1991')
# ax3.set_ylim([60,490])
# ax3.set_xlim([0,4.5])






#
# index2 = np.where((timeS > DT.datetime(1992,8,1)) & (timeS < DT.datetime(1998,1,1)))
# newTime2 = [date - timeS[index2[0][0]] for date in timeS[index2[0]]]
# yearTime2 = [daysPassed.days/365 for daysPassed in newTime2]
# phiTemp2 = phit2S[index2,mode]
# # phiIndex2 = np.where((phiTemp2[0]<75))
# # phiTemp2[0][phiIndex2[0]] = phiTemp2[0][phiIndex2[0]] + 360
# # ax2.scatter(yearTime,np.unwrap(phit2S[index1,mode]),temp[index1],phit2S[index1,mode])
# # ax2.scatter(yearTime2,phiTemp2,RS[index2,mode]*10,phit2S[index2,mode])
# ax2.scatter(yearTime2,phiTemp2,RS[index2,mode]*10,'r')
#
#
#
# index2 = np.where((timeS > DT.datetime(1990,9,1)) & (timeS < DT.datetime(1994,9,1)))
# newTime2 = [date - timeS[index2[0][0]] for date in timeS[index2[0]]]
# yearTime2 = [daysPassed.days/365 for daysPassed in newTime2]
# phiTemp2 = phit2S[index2,mode]
# # phiIndex2 = np.where((phiTemp2[0]<75))
# # phiTemp2[0][phiIndex2[0]] = phiTemp2[0][phiIndex2[0]] + 360
# # ax2.scatter(yearTime,np.unwrap(phit2S[index1,mode]),temp[index1],phit2S[index1,mode])
# ax2.scatter(yearTime2,phiTemp2,RS[index2,mode]*10,phit2S[index2,mode])




mode = 1

fig, ax = plt.subplots(2,2)

ax[0,0].scatter(xinterpS, SS[:,mode],12,theta2S[:,mode],cmap=cmap)
ax[0,0].set_ylabel('Spatial Magnitude (m)')
ax[0,0].set_xlabel('Cross-shore (m)')
ax[1,0].scatter(xinterpS, theta2S[:,mode],12,theta2S[:,mode],cmap=cmap)
ax[1,0].set_ylabel('Spatial Phase (deg)')
ax[1,0].set_xlabel('Cross-shore (m)')
ax[0,1].scatter(timeS,RS[:,mode],12,phit2S[:,mode])
ax[0,1].set_ylabel('Temporal Magnitude (m)')
ax[0,1].set_xlabel('Time')
ax[1,1].scatter(timeS,phit2S[:,mode],12,phit2S[:,mode])
ax[1,1].set_ylabel('Temporal Phase (deg)')
ax[1,1].set_xlabel('Time')


mode = 2

fig, ax = plt.subplots(2,2)

ax[0,0].scatter(xinterpS, SS[:,mode],12,theta2S[:,mode],cmap=cmap)
ax[0,0].set_ylabel('Spatial Magnitude (m)')
ax[0,0].set_xlabel('Cross-shore (m)')
ax[1,0].scatter(xinterpS, theta2S[:,mode],12,theta2S[:,mode],cmap=cmap)
ax[1,0].set_ylabel('Spatial Phase (deg)')
ax[1,0].set_xlabel('Cross-shore (m)')
ax[0,1].scatter(timeS,RS[:,mode],12,phit2S[:,mode])
ax[0,1].set_ylabel('Temporal Magnitude (m)')
ax[0,1].set_xlabel('Time')
ax[1,1].scatter(timeS,phit2S[:,mode],12,phit2S[:,mode])
ax[1,1].set_ylabel('Temporal Phase (deg)')
ax[1,1].set_xlabel('Time')


#
# mode = 0
# southUnwrapped = np.unwrap(phitS[:,0],discont=2.*np.pi/2)
# deltaPhi = np.diff(southUnwrapped)
# plt.figure()
# # plt.plot(timeS[1:],deltaPhi,'.',color='r')
# # plt.scatter(timeS[1:],deltaPhi,15,RS[1:,mode])
# longInts = np.where((days2<90))
# hsSub = hs90[longInts]
# weSub = we90[longInts]
# lwpSub = lwp90[longInts]
# cumWeSub = cumWe[longInts]
# deltaPhiSub = deltaPhi[longInts]
# normWeSub = normWe[longInts]
# mags = RS[1:,mode]
# magsSub = mags[longInts]
# daysSub = days2[longInts]
# fvSub = fv50[longInts]
# hlSub = hl50[longInts]
# hsSub2 = hs50[longInts]
# x = hsSub
# y = deltaPhiSub
# plt.scatter(x,y,15,fvSub)
# plt.ylim([-1.5,2])
# plt.ylabel('d(phi)/dt (radians)')
# plt.xlabel('Hs (90th percentile)')
# coef = np.polyfit(x,y,1)
# poly1d_fn = np.poly1d(coef)
# # poly1d_fn is now a function which takes in x and returns an estimate for y
# plt.plot(x, poly1d_fn(x), '--k')
#




def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )


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
import datetime as DT
hsCombined = np.append(Hs,hs26m)
hsSmooth = moving_average(hsCombined,3)
tpCombined = np.append(Tp,tp26m)
dmCombined = np.append(Dm,dm26m)
tWave = [DT.datetime.fromtimestamp(x) for x in timeWave]
tC = np.append(np.array(tWave),tWave26m)

badDirs = np.where((dmCombined > 360))
dmCombined[badDirs] = dmCombined[badDirs]*np.nan
badtp = np.where((tpCombined < 1))
tpCombined[badtp] = tpCombined[badtp]*np.nan

waveNorm = dmCombined - 72
neg = np.where((waveNorm > 180))
waveNorm[neg[0]] = waveNorm[neg[0]]-360
offpos = np.where((waveNorm>90))
offneg = np.where((waveNorm<-90))
waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
waveNorm[offneg[0]] = waveNorm[offneg[0]]*0

Lo = (9.81/(2*np.pi)) * np.square(tpCombined)
Ir = 0.122/(np.sqrt((hsCombined/Lo)))
HoverL = hsCombined/Lo
lwpC = 1025*np.square(hsCombined)*tpCombined*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
weC = np.square(hsCombined)*tpCombined
ws = (2.65-1)*9.81*np.square(0.00015)/(18*0.000001)
fV = hsCombined/(ws*tpCombined)



hsList = []
tpList = []
dmList = []
lwpList = []
avgLWP = []
weList = []
avgWE = []
fvList = []
irList = []
hlList = []
tList = []
hList = []
dList = []
avgFV = []
t2List = []
avgHs = []
for xx in range(len(timeS)-1):
    t1 = timeS[xx]
    t2 = timeS[xx+1]
    tempWave = np.where((tC < t2) & (tC > t1))

    t2List.append(t2.month)
    dList.append(len(tempWave[0])/24)
    tList.append(tC[tempWave])
    hList.append(len(tempWave[0]))
    hsList.append(hsCombined[tempWave])
    avgHs.append(np.nanmean(hsCombined[tempWave]))
    tpList.append(tpCombined[tempWave])
    dmList.append(dmCombined[tempWave])
    lwpList.append(lwpC[tempWave])
    avgLWP.append(np.nanmean(np.abs(lwpC[tempWave])))
    weList.append(weC[tempWave])
    # avgWE.append(np.nansum(weC[tempWave])/len(tempWave[0]))
    avgWE.append(np.nanmean(weC[tempWave]))

    fvList.append(fV[tempWave])
    avgFV.append(np.nanmean(fV[tempWave]))
    irList.append(Ir[tempWave])
    hlList.append(HoverL[tempWave])





mode1Phase = moving_average(np.unwrap(phit2S[:,0]*(np.pi/180)),n=5)
mode2Phase = moving_average(np.unwrap(phit2S[:,1]*(np.pi/180)),n=5)
mode3Phase = moving_average(np.unwrap(phit2S[:,2]*(np.pi/180)),n=5)
mode1Phase = np.hstack((mode1Phase[0],np.hstack((mode1Phase[0],np.hstack((mode1Phase[:],np.hstack((mode1Phase[-1],mode1Phase[-1]))))))))
mode2Phase = np.hstack((mode2Phase[0],np.hstack((mode2Phase[0],np.hstack((mode2Phase[:],np.hstack((mode2Phase[-1],mode2Phase[-1]))))))))
mode3Phase = np.hstack((mode3Phase[0],np.hstack((mode3Phase[0],np.hstack((mode3Phase[:],np.hstack((mode3Phase[-1],mode3Phase[-1]))))))))

mode1Mag = moving_average(RS[:,0],n=5)
mode2Mag = moving_average(RS[:,1],n=5)
mode3Mag = moving_average(RS[:,2],n=5)
mode1Mag = np.hstack((mode1Mag[0],np.hstack((mode1Mag[0],np.hstack((mode1Mag[:],np.hstack((mode1Mag[-1],mode1Mag[-1]))))))))
mode2Mag = np.hstack((mode2Mag[0],np.hstack((mode2Mag[0],np.hstack((mode2Mag[:],np.hstack((mode2Mag[-1],mode2Mag[-1]))))))))
mode3Mag = np.hstack((mode3Mag[0],np.hstack((mode3Mag[0],np.hstack((mode3Mag[:],np.hstack((mode3Mag[-1],mode3Mag[-1]))))))))

smoothingHS = moving_average(np.array(avgHs),n=5)
smoothingHS = np.hstack((smoothingHS[0],np.hstack((smoothingHS[0],np.hstack((smoothingHS[:],np.hstack((smoothingHS[-1],smoothingHS[-1]))))))))


mode1Wrap = (mode1Phase + np.pi) % (2 * np.pi) - np.pi
mode2Wrap = (mode2Phase + np.pi) % (2 * np.pi) - np.pi
mode3Wrap = (mode3Phase + np.pi) % (2 * np.pi) - np.pi

mode1avgMag = np.nanmean(np.vstack((mode1Mag[0:-1],mode1Mag[1:])),axis=0)
deltaPhi = np.diff(mode1Phase)
mode2avgMag = np.nanmean(np.vstack((mode2Mag[0:-1],mode2Mag[1:])),axis=0)
deltaPhi2 = np.diff(mode2Phase)

dArray = np.array(dList)
shortEnough = np.where((dArray<60))

mode1PhiChange = deltaPhi[shortEnough]/np.array(hList)[shortEnough]
mode2PhiChange = deltaPhi2[shortEnough]/np.array(hList)[shortEnough]

mode1PhaseShort = mode1Wrap[shortEnough]
mode2PhaseShort = mode2Wrap[shortEnough]
mode1MagAvg = mode1avgMag[shortEnough]
mode1LWP = np.array(avgLWP)[shortEnough]
mode1WE = np.array(avgWE)[shortEnough]
mode1FV = np.array(avgFV)[shortEnough]
mode1HS = smoothingHS[shortEnough]
fig10 = plt.figure(figsize=(12,12))
ax100 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
sc1 = ax100.scatter(mode1PhiChange,mode2PhiChange,10*(mode1PhaseShort+np.pi),mode2PhaseShort+np.pi,cmap='rainbow')#,vmin=0,vmax=500)
fig10.colorbar(sc1,ax=ax100)
plt.show()

hsOrdinaltime = np.array([i.toordinal() for i in tC])
time5Ordinal = np.array([i.toordinal() for i in time3])
timeOrdinal = np.array([i.toordinal() for i in time])
avgMonthBeforeHs = np.nan * np.ones((len(time3),))
avgMonthBeforeWE = np.nan * np.ones((len(time3),))
avgMonthBeforeFV = np.nan * np.ones((len(time3),))
avgMonthBeforeLWP = np.nan * np.ones((len(time3),))

for ii in range(len(time3)):
    # beforeIndex = np.where((tC > DT.datetime(time5[ii].year,time5[ii].month,time5[ii].day)) & (tC < time5[ii]))
    beforeIndex = np.where((hsOrdinaltime > (time5Ordinal[ii]-30)) & (hsOrdinaltime < (time5Ordinal[ii])))

    tempHs = hsCombined[beforeIndex]
    avgMonthBeforeHs[ii] = np.nanmean(tempHs)
    tempWE = weC[beforeIndex]
    avgMonthBeforeWE[ii] = np.nanmean(tempWE)
    tempLWP = lwpC[beforeIndex]
    avgMonthBeforeLWP[ii] = np.nanmean(np.abs(tempLWP))
    tempFV = fV[beforeIndex]
    avgMonthBeforeFV[ii] = np.nanmean(tempFV)

maxMonthBeforeWE = np.nan * np.ones((len(time),))
for ii in range(len(time)):
    # beforeIndex = np.where((tC > DT.datetime(time5[ii].year,time5[ii].month,time5[ii].day)) & (tC < time5[ii]))
    beforeIndex = np.where((hsOrdinaltime > (timeOrdinal[ii]-30)) & (hsOrdinaltime < (timeOrdinal[ii])))
    tempWE = weC[beforeIndex]
    # tempWE = hsCombined[beforeIndex]

    maxMonthBeforeWE[ii] = np.nanmax(tempWE)


plt.figure()
plt.plot(time,maxMonthBeforeWE)
plt.show()


mode1Phase = np.unwrap(phit2S[:,0]*(np.pi/180))
mode2Phase = np.unwrap(phit2S[:,1]*(np.pi/180))

# # Wave Height
# lowBreak = 0.65
# medLowBreaka = 0.6
# medLowBreakb = 0.9
# medHighBreaka = 0.85
# medHighBreakb = 1.15
# highBreaka = 1.1
# highBreakb = 1.35
# highesta = 1.3
lowBreak = 6
medLowBreaka = 6
medLowBreakb = 7.6
medHighBreaka = 7.6
medHighBreakb = 9.1
highBreaka = 9.1
highBreakb = 17
highesta = 18

binClassHs = np.nan * np.ones((len(time3),))
for ii in range(len(time3)):
    hs1 = avgMonthBeforeFV[ii]
    if hs1 > medLowBreaka and hs1 < medLowBreakb:
        binClassHs[ii] = 1
    elif hs1 >= medHighBreaka and hs1 < medHighBreakb:
        binClassHs[ii] = 2
    elif hs1 >= highBreaka and hs1 < highBreakb:
        binClassHs[ii] = 3
    elif hs1 >= highesta:
        binClassHs[ii] = 4
    elif hs1 < lowBreak:
        binClassHs[ii] = 0


lowBreak = 100
medLowBreaka = 100
medLowBreakb = 200
medHighBreaka = 200
medHighBreakb = 350
highBreaka = 350
highBreakb = 600
highesta = 802

binClassLWP = np.nan * np.ones((len(time7),))
for ii in range(len(time7)):
    hs1 = avgMonthBeforeLWP[ii]
    if hs1 > medLowBreaka and hs1 < medLowBreakb:
        binClassLWP[ii] = 1
    elif hs1 >= medHighBreaka and hs1 < medHighBreakb:
        binClassLWP[ii] = 2
    elif hs1 >= highBreaka and hs1 < highBreakb:
        binClassLWP[ii] = 3
    elif hs1 >= highesta:
        binClassLWP[ii] = 4
    elif hs1 < lowBreak:
        binClassLWP[ii] = 0


lowBreak = 50
medLowBreaka = 50
medLowBreakb = 100
medHighBreaka = 100
medHighBreakb = 150
highBreaka = 150
highBreakb = 250
highesta = 250
# lowBreak = 1
# medLowBreaka = 1
# medLowBreakb = 1.5
# medHighBreaka = 1.5
# medHighBreakb = 2
# highBreaka = 2
# highBreakb = 4.5
# highesta = 4.5

binClassWE = np.nan * np.ones((len(time),))
for ii in range(len(time)):
    hs1 = maxMonthBeforeWE[ii]
    if hs1 > medLowBreaka and hs1 < medLowBreakb:
        binClassWE[ii] = 1
    elif hs1 >= medHighBreaka and hs1 < medHighBreakb:
        binClassWE[ii] = 2
    elif hs1 >= highBreaka and hs1 < highBreakb:
        binClassWE[ii] = 3
    elif hs1 >= highesta:
        binClassWE[ii] = 4
    elif hs1 < lowBreak:
        binClassWE[ii] = 0



from itertools import groupby
grouped_L = [(k, sum(1 for i in g)) for k,g in groupby(binClassHs)]

def repeatingNumbers(numList):
    i = 0
    output = []
    while i < len(numList) - 1:
        n = numList[i]
        startIndex = i
        while i < len(numList) - 1 and numList[i] == numList[i + 1]:
            i = i + 1

        endIndex = i

        output.append([n, startIndex, endIndex, (endIndex-startIndex + 1)])
        print("{0} >> {1}".format(n, [startIndex, endIndex]))
        i = i + 1


    return output

groupedOutput = repeatingNumbers(binClassHs)
lowWaves = []
medLowWaves = []
medHighWaves = []
bigWaves = []
highestWaves = []

for tt in range(len(groupedOutput)):
    if groupedOutput[tt][3] > 1:
        if groupedOutput[tt][0] == 0:
            lowWaves.append([groupedOutput[tt][1],groupedOutput[tt][2]])
        if groupedOutput[tt][0] == 1:
            medLowWaves.append([groupedOutput[tt][1], groupedOutput[tt][2]])
        if groupedOutput[tt][0] == 2:
            medHighWaves.append([groupedOutput[tt][1],groupedOutput[tt][2]])
        if groupedOutput[tt][0] == 3:
            bigWaves.append([groupedOutput[tt][1], groupedOutput[tt][2]])

groupedOutputLWP = repeatingNumbers(binClassLWP)
lowLWP = []
medLowLWP = []
medHighLWP = []
bigLWP  = []
highLWP = []
highestLWP = []

for tt in range(len(groupedOutputLWP)):
    if groupedOutputLWP[tt][3] > 1:
        if groupedOutputLWP[tt][0] == 0:
            lowLWP.append([groupedOutputLWP[tt][1],groupedOutputLWP[tt][2]])
        if groupedOutputLWP[tt][0] == 1:
            medLowLWP.append([groupedOutputLWP[tt][1], groupedOutputLWP[tt][2]])
        if groupedOutputLWP[tt][0] == 2:
            medHighLWP.append([groupedOutputLWP[tt][1],groupedOutputLWP[tt][2]])
        if groupedOutputLWP[tt][0] == 3:
            highLWP.append([groupedOutputLWP[tt][1], groupedOutputLWP[tt][2]])
        if groupedOutputLWP[tt][0] == 4:
            highestLWP.append([groupedOutputLWP[tt][1], groupedOutputLWP[tt][2]])


groupedOutputWE = repeatingNumbers(binClassWE)
lowWE = []
medLowWE = []
medHighWE = []
bigWE = []
highestWE = []
highWE = []

for tt in range(len(groupedOutputWE)):
    if groupedOutputWE[tt][3] >= 1:
        if groupedOutputWE[tt][0] == 0:
            lowWE.append([groupedOutputWE[tt][1],groupedOutputWE[tt][2]])
        if groupedOutputWE[tt][0] == 1:
            medLowWE.append([groupedOutputWE[tt][1], groupedOutputWE[tt][2]])
        if groupedOutputWE[tt][0] == 2:
            medHighWE.append([groupedOutputWE[tt][1],groupedOutputWE[tt][2]])
        if groupedOutputWE[tt][0] == 3:
            highWE.append([groupedOutputWE[tt][1], groupedOutputWE[tt][2]])
        if groupedOutputWE[tt][0] == 4:
            highestWE.append([groupedOutputWE[tt][1], groupedOutputWE[tt][2]])



from scipy.interpolate import interp1d


fig = plt.figure(figsize=(10,5))
ax1 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=2)
xg,tg = np.meshgrid(xinterpS,timeS)
plt0 = ax1.pcolor(tg,xg,(alllinesS-np.mean(alllinesS, axis=0)), vmin=-1.8, vmax=1.8,cmap='bwr')
ax1.set_ylabel('Cross-shore (m)',fontsize=24)
ax1.set_title('Sand Bar Migrations at the FRF',fontsize=24)
fig.set_size_inches(22, 8)
plt.show()
#plt.savefig("sample.png", dpi=400)

mode = 0

fig = plt.figure(figsize=(10,5))
ax1 = plt.subplot2grid((4,2),(0,0),rowspan=1,colspan=2)
xg,tg = np.meshgrid(xinterpS,timeS)
plt0 = ax1.pcolor(tg,xg,(alllinesS-np.mean(alllinesS, axis=0)), vmin=-1.8, vmax=1.8,cmap='bwr')

ax = plt.subplot2grid((4,2),(1,0),rowspan=1,colspan=2,sharex=ax1)
# temp = RS[:,mode]*10
# temp1 = phit2S[:,mode]
temp = mode1Mag*10
temp1 = mode1Wrap
# index = np.where((temp1<0))
# temp1[index[0]] = temp1[index[0]]+360
# ax.scatter(timeS,np.unwrap(phit2S[:,mode]*(np.pi/180)),temp,phit2S[:,mode])
ax.scatter(time,temp1,temp,mode1Mag,cmap='Reds')

mode = 1
ax2 = plt.subplot2grid((4,2),(2,0),rowspan=1,colspan=2, sharex=ax1)
# temp = RS[:,1]*10
# temp1 = phit2S[:,mode]
temp = mode2Mag*10
temp1 = mode2Wrap
# index = np.where((temp1<0))
# temp1[index[0]] = temp1[index[0]]+360
ax2.scatter(time,temp1,temp,mode2Mag,cmap='Blues')

ax3 = plt.subplot2grid((4,2),(3,0),rowspan=1,colspan=2, sharex=ax1)
dayAvg = 60
hsSmoother = moving_average(hsCombined,24*dayAvg)
ax3.plot(tC[int(24*dayAvg/2-1):-int(24*dayAvg/2)],hsSmoother)
dayAvg = 90
hsSmoother = moving_average(hsCombined,24*dayAvg)
ax3.plot(tC[int(24*dayAvg/2-1):-int(24*dayAvg/2)],hsSmoother)
ax3.plot(time7,avgMonthBeforeHs,'o-')
ax.set_xlim([timeS[2], timeS[-2]])
ax1.set_xlim([timeS[2], timeS[-2]])
ax2.set_xlim([timeS[2], timeS[-2]])
ax3.set_xlim([timeS[2], timeS[-2]])

plt.show()




import matplotlib.colors as mcolors


fig = plt.figure(figsize=(12,12))
ax100 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
con1 = ax100.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
coloredlines = cm.rainbow(np.linspace(0,1,len(lowWaves)))

for ff in range(len(lowWaves)):
    si = lowWaves[ff][0]
    ei = lowWaves[ff][1]
    # si = medLowWaves[ff][0]
    # ei = medLowWaves[ff][1]
    # si = medHighWaves[ff][0]
    # ei = medHighWaves[ff][1]
    # si = bigWaves[ff][0]
    # ei = bigWaves[ff][1]
    indexedWaves = avgMonthBeforeFV[si-1:ei+2]
    indexedTime = time7[si-1:ei+2]

    indexedOrdinalTime = [i.toordinal() for i in indexedTime]
    indexedTimeMonth = [i.month for i in indexedTime]
    mode1prop = mode1Phase7[si-1:ei+2]
    mode2prop = mode2Phase7[si-1:ei+2]
    mode1propWrap = np.unwrap((mode1prop + np.pi) % (2 * np.pi) - np.pi)
    mode2propWrap = np.unwrap((mode2prop + np.pi) % (2 * np.pi) - np.pi)
    x1 = np.sort(mode1propWrap*180/np.pi)
    y1 = mode2propWrap*180/np.pi
    mags = mode2Mag7[si-1:ei+2]

    if y1[-1] < -20:
        y1 = y1 + 360
    newx = np.arange(x1[0],x1[-1],1)
    f = interp1d(x1,y1,kind='linear')
    newy = f(newx)

    if x1[-1] < 60 and y1[0]>-100:
        lines1 = ax100.plot(newx,newy,'-',linewidth=2,color='red',label='{}'.format(indexedTime[0].year))
        scater = ax100.scatter(x1,y1,mags*40,'k')
        add_arrow(lines1[0], color='r', size=25)
    elif x1[-1] > 60:
        if indexedTime[0].year < 2014 and indexedTime[0].year > 1981:
            if x1[0] < 73 and y1[0] < -65:
                print('skipping')
            else:
                lines1 = ax100.plot(newx,newy,'-',linewidth=2,color='blue',label='{}'.format(indexedTime[0].year))
                scater = ax100.scatter(x1, y1, mags * 40, 'k')

                add_arrow(lines1[0], color='blue', size=25)
        #
            # if indexedTime[0].year == 1993:
            #
            # else:
            #
        # elif indexedTime[0].year == 1993:
        #     lines1 = ax100.plot(newx,newy,'-',linewidth=2,color='green',label='{}'.format(indexedTime[0].year))
        #     add_arrow(lines1[0], color='blue', size=25)
            # lines1 = ax100.plot(newx,newy,'-',linewidth=2,color=coloredlines[ff,:],label='{}'.format(indexedTime[0].year))
            #scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedWaves,cmap='rainbow',zorder=3,vmin=3.5,vmax=6.5)
            # scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedTimeMonth,cmap='twilight_shifted',zorder=3,vmin=1,vmax=12)

ax100.set_xlim([-180, 340])
ax100.set_ylim([-160, 360])
# cb = fig.colorbar(scatter1,ax=ax100)
ax100.set_xlabel(r'Offshore Propagation $\Longrightarrow$',fontsize=16)
ax100.set_ylabel(r'Onshore Propagation $\Longrightarrow$',fontsize=16)
ax100.set_title('Non-dimensional Fall Velocity < 6')
# plt.legend()
plt.show()




fig = plt.figure(figsize=(18,12))

# ax100 = plt.subplot2grid((1,3),(0,0),rowspan=1,colspan=1)
ax100 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)

con1 = ax100.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
coloredlines = cm.rainbow(np.linspace(0,1,len(lowWaves)))

for ff in range(len(lowWaves)):
    si = lowWaves[ff][0]
    ei = lowWaves[ff][1]
    # si = medLowWaves[ff][0]
    # ei = medLowWaves[ff][1]
    # si = medHighWaves[ff][0]
    # ei = medHighWaves[ff][1]
    # si = bigWaves[ff][0]
    # ei = bigWaves[ff][1]
    #indexedWaves = avgMonthBeforeFV[si-1:ei+2]
    indexedTime = time3[si-2:ei]

    indexedOrdinalTime = [i.toordinal() for i in indexedTime]
    indexedTimeMonth = [i.month for i in indexedTime]
    mode1prop = mode1Phase3WA[si-2:ei]
    mode2prop = mode2Phase3WA[si-2:ei]
    mode1propWrap = np.unwrap((mode1prop + np.pi) % (2 * np.pi) - np.pi)
    mode2propWrap = np.unwrap((mode2prop + np.pi) % (2 * np.pi) - np.pi)
    x1 = np.sort(mode1propWrap*180/np.pi)
    y1 = mode2propWrap*180/np.pi
    mags = mode1Mag[si:ei+2]

    y1 = np.sort(y1)
    rise = np.abs(y1[-1]-y1[0])
    run = np.abs(x1[-1]-x1[0])

    slope = rise/run

    colormap = cm.coolwarm
    # colormap = cm.jet
    #normalize = mcolors.Normalize(vmin=0.33, vmax=3)
    normalize = mcolors.LogNorm(vmin=0.5, vmax=3)

    s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
    colorOfLine = colormap(normalize(slope))

    if y1[-1] < -20:
        y1 = y1 + 360
    newx = np.arange(x1[0],x1[-1],1)
    f = interp1d(x1,y1,kind='linear')
    newy = f(newx)

    if indexedTime[0].year < 2005 or indexedTime[0].year > 2008:
        if indexedTime[0].year < 2018:
            if np.nanmean(mags) > 0.4:

                if x1[-1] < 60 and y1[0]>-100:

                    #else:
                    lines1 = ax100.plot(newx,newy,'-',linewidth=2,color=colorOfLine,label='{}'.format(indexedTime[0].year))
                        # lines1 = ax100.plot(newx,newy,'-',linewidth=2,color='red',label='{}'.format(indexedTime[0].year))
                    #scater = ax100.scatter(x1,y1,mags*50,color=colorOfLine)
                    add_arrow(lines1[0], color='k', size=25)
                elif x1[-1] > 60:
                    if x1[-1] > 200 and y1[-1]>200:
                        lines1 = ax100.plot(newx-360,newy,'-',linewidth=2,color=colorOfLine,label='{}'.format(indexedTime[0].year))
                        ## lines1 = ax100.plot(newx,newy,'-',linewidth=2,color='red',label='{}'.format(indexedTime[0].year))
                        #scater = ax100.scatter(x1-360,y1,mags*50,color=colorOfLine)
                        add_arrow(lines1[0], color='k', size=25)
                    else:
                        if indexedTime[0].year < 2014 and indexedTime[0].year > 1981:
                            if x1[0] < 73 and y1[0] < -65:
                                print('skipping')
                            else:
                                lines1 = ax100.plot(newx,newy,'-',linewidth=2,color=colorOfLine,label='{}'.format(indexedTime[0].year))
                                ##scater = ax100.scatter(x1, y1, mags * 40, 'k')
                                #scater = ax100.scatter(x1, y1, mags * 50, color=colorOfLine)

                            add_arrow(lines1[0], color='k', size=25)
                        #
                    # if indexedTime[0].year == 1993:
                #
                # else:
                #
        # elif indexedTime[0].year == 1993:
        #     lines1 = ax100.plot(newx,newy,'-',linewidth=2,color='green',label='{}'.format(indexedTime[0].year))
        #     add_arrow(lines1[0], color='blue', size=25)
            # lines1 = ax100.plot(newx,newy,'-',linewidth=2,color=coloredlines[ff,:],label='{}'.format(indexedTime[0].year))
            #scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedWaves,cmap='rainbow',zorder=3,vmin=3.5,vmax=6.5)
            # scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedTimeMonth,cmap='twilight_shifted',zorder=3,vmin=1,vmax=12)

ax100.set_xlim([-180, 280])
ax100.set_ylim([-160, 320])
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# axins1 = inset_axes(ax100,
#                     width="50%",  # width = 50% of parent_bbox width
#                     height="5%",  # height : 5%
#                     loc='upper right')
# # cb = fig.colorbar(s_map, cax=axins1, orientation='horizontal')#,yticks=[0.5,1,2])#,yticklabels=['0.5','1','2'])
# # cb.ax.set_yticklabels(['0.5','1','2'], fontsize=16, weight='bold')
# # axins1.xaxis.set_ticks_position("bottom")
# # axins1.set_xticks([0.5,1,2])
# # axins1.set_xticklabels(['0.5','1','2'], fontsize=14)#, weight='bold')
# # cb = fig.colorbar(lines1[0],ax=ax100)
# #for l in axins1.xaxis.get_ticklabels():
# #    l.set_weight("bold")
# #    l.set_fontsize(12)
ax100.set_xlabel(r'Offshore Propagation $\Longrightarrow$',fontsize=12)
ax100.set_ylabel(r'Onshore Propagation $\Longrightarrow$',fontsize=12)
ax100.set_title('Non-dimensional Fall Velocity < 6')
# plt.legend()
plt.show()






# fig = plt.figure(figsize=(12,12))
ax100b = plt.subplot2grid((1,3),(0,1),rowspan=1,colspan=1)
con1 = ax100b.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')

for ff in range(len(bigWaves)):
    # si = lowWaves[ff][0]
    # ei = lowWaves[ff][1]
    # si = medLowWaves[ff][0]
    # ei = medLowWaves[ff][1]
    # si = medHighWaves[ff][0]
    # ei = medHighWaves[ff][1]
    si = bigWaves[ff][0]
    ei = bigWaves[ff][1]
    if ff == 0:
        indexedWaves = avgMonthBeforeFV[si:ei+1]
        indexedTime = time3[si:ei+1]
        indexedOrdinalTime = [i.toordinal() for i in indexedTime]
        indexedTimeMonth = [i.month for i in indexedTime]
        mode1prop = mode1Phase3WA[si:ei+1]
        mode2prop = mode2Phase3WA[si:ei+1]
        mode1propWrap = np.unwrap((mode1prop + np.pi) % (2 * np.pi) - np.pi)
        mode2propWrap = np.unwrap((mode2prop + np.pi) % (2 * np.pi) - np.pi)
        x1 = np.sort(mode1propWrap*180/np.pi)
        y1 = mode2propWrap*180/np.pi
        mags = mode2Mag[si:ei+1]

    else:
        indexedWaves = avgMonthBeforeFV[si-1:ei]
        indexedTime = time3[si-1:ei]
        indexedOrdinalTime = [i.toordinal() for i in indexedTime]
        indexedTimeMonth = [i.month for i in indexedTime]
        mode1prop = mode1Phase3WA[si-1:ei]
        mode2prop = mode2Phase3WA[si-1:ei]
        mode1propWrap = np.unwrap((mode1prop + np.pi) % (2 * np.pi) - np.pi)
        mode2propWrap = np.unwrap((mode2prop + np.pi) % (2 * np.pi) - np.pi)
        x1 = np.sort(mode1propWrap*180/np.pi)
        y1 = mode2propWrap*180/np.pi
        mags = mode2Mag[si-1:ei]
    y1 = np.sort(y1)
    rise = np.abs(y1[-1]-y1[0])
    run = np.abs(x1[-1]-x1[0])

    slope = rise/run

    colormap = cm.coolwarm
    # colormap = cm.jet
    #normalize = mcolors.Normalize(vmin=0.33, vmax=3)
    normalize = mcolors.LogNorm(vmin=0.5, vmax=3)

    s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
    colorOfLine = colormap(normalize(slope))
    if indexedTime[0].year < 2005 or indexedTime[0].year > 2008:
        if indexedTime[0].year < 2018:

            if y1[0] < -110:
               y1 = y1 + 360
            #newx = np.arange(x1[0],x1[-1],1)
            #f = interp1d(x1,y1,kind='linear')
            #newy = f(newx)
            if x[-1] > 155 and y[-1] > 100:
                x1 = x1-360

            # if x1[0] < -100:
            #     lines1 = ax100.plot(x1+360, y1, '-', linewidth=2, color=colorOfLine,
            #                         label='{}'.format(indexedTime[0].year))
            #
            #     # lines1 = ax100.plot(newx+360, newy, 'k-', linewidth=2)
            #     scater = ax100.scatter(x1+360, y1, mags * 50, color=colorOfLine)
            #
            # else:
            # lines1 = ax100.plot(x1,y1,'k-')
                #lines1 = ax100.plot(newx,newy,'k-',linewidth=2)
            lines1 = ax100b.plot(x1, y1, '-', linewidth=2, color=colorOfLine, label='{}'.format(indexedTime[0].year))

            #scater = ax100b.scatter(x1, y1, mags * 50, color=colorOfLine)

            #scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedWaves,cmap='rainbow',zorder=3,vmin=3.5,vmax=6.5)
            #add_arrow(lines1[0], color='k', size=25)
            # scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedTimeMonth,cmap='twilight_shifted',zorder=3,vmin=1,vmax=12)



ax100b.set_xlim([-180, 280])
ax100b.set_ylim([-160, 320])
# cb = fig.colorbar(scatter1,ax=ax100)
ax100b.set_xlabel(r'Offshore Propagation $\Longrightarrow$',fontsize=12)
# ax100b.set_ylabel(r'Onshore Propagation $\Longrightarrow$')
ax100b.set_title('Non-dimensional Fall Velocity > 9')
# ax100b.legend()
# plt.show()








colorparams =  avgMonthBeforeHs
# Choose a colormap
colormap = cm.rainbow
#colormap = cm.jet
normalize = mcolors.Normalize(vmin=1.3, vmax=1.8)
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
# s_map.set_array(colorparams)

# fig = plt.figure(figsize=(12,12))
ax100c = plt.subplot2grid((1,3),(0,2),rowspan=1,colspan=1)
con1 = ax100c.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
coloredlines = cm.rainbow(np.linspace(0,1,len(medHighLWP)))

for ff in range(len(medHighLWP)):
    # si = lowLWP[ff][0]
    # ei = lowLWP[ff][1]
    # si = medLowLWP[ff][0]
    # ei = medLowLWP[ff][1]
    si = medHighLWP[ff][0]
    ei = medHighLWP[ff][1]
    # si = bigWaves[ff][0]
    # ei = bigWaves[ff][1]
    if si == 0:
        si = 1
    indexedWaves = avgMonthBeforeFV[si-1:ei+1]
    indexedTime = time5[si-1:ei+1]
    indexedOrdinalTime = [i.toordinal() for i in indexedTime]
    indexedTimeMonth = [i.month for i in indexedTime]
    mode1prop = mode1Phase5WA[si-1:ei+1]
    mode2prop = mode2Phase5WA[si-1:ei+1]
    mode1propWrap = np.unwrap((mode1prop + np.pi) % (2 * np.pi) - np.pi)
    mode2propWrap = np.unwrap((mode2prop + np.pi) % (2 * np.pi) - np.pi)
    x1 = np.sort(mode1propWrap*180/np.pi)
    y1 = mode2propWrap*180/np.pi

    mags = mode2Mag[si - 1:ei+1]

    y1 = np.sort(y1)
    rise = np.abs(y1[-1]-y1[0])
    run = np.abs(x1[-1]-x1[0])

    slope = rise/run

    colormap = cm.bwr
    # colormap = cm.jet
    #normalize = mcolors.Normalize(vmin=0.33, vmax=3)
    normalize = mcolors.LogNorm(vmin=0.5, vmax=3)

    s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
    colorOfLine = colormap(normalize(slope))
    color = colormap(normalize(np.max(avgMonthBeforeHs[si-1:ei+1])))
    if y1[-1] < -20:
        y1 = y1 + 360
    newx = np.arange(x1[0],x1[-1],1)
    f = interp1d(x1,y1,kind='linear')
    newy = f(newx)
    if indexedTime[0].year < 2004 or indexedTime[0].year > 2010:
        if y1[-1]-y1[0] > 200:
            print('skip')
        else:
            if newx[-1] > 220:
                print('i found it')
                lines1 = ax100c.plot(newx-360,newy,'-',linewidth=2,color='k')#color,label='{}'.format(np.max(avgMonthBeforeHs[si-1:ei+1])))
                # lines1 = ax100.plot(x1, y1, '-', linewidth=2, color='k', label='{}'.format(indexedTime[0].year))
                #scater = ax100.scatter(x1-360, y1, mags * 50, color=colorOfLine)

                add_arrow(lines1[0], color='k', size=25)
            else:
                if newy[-1] > 270:
                    lines1 = ax100c.plot(newx, newy-360, '-', linewidth=2, color='k')#,
                                        #label='{}'.format(np.max(avgMonthBeforeHs[si - 1:ei + 1])))
                    #scater = ax100.scatter(x1, y1-360, mags * 50, color=colorOfLine)

                    add_arrow(lines1[0], color='k', size=25)
                else:
                    lines1 = ax100c.plot(newx, newy, '-', linewidth=2, color='k')#color,
                                    #label='{}'.format(np.max(avgMonthBeforeHs[si - 1:ei + 1])))
                    #scater = ax100.scatter(x1, y1, mags * 50, color=colorOfLine)

                    # scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedWaves,cmap='rainbow',zorder=3,vmin=3.5,vmax=6.5)
                    add_arrow(lines1[0], color='k', size=25)
        # lines1 = ax100.plot(x1,y1,'k-')
        # scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedTimeMonth,cmap='twilight_shifted',zorder=3,vmin=1,vmax=12)

ax100c.set_xlim([-180, 280])
ax100c.set_ylim([-160, 320])
# cb = fig.colorbar(scatter1,ax=ax100)
ax100c.set_xlabel(r'Offshore Propagation $\Longrightarrow$',fontsize=12)
# ax100c.set_ylabel(r'Onshore Propagation $\Longrightarrow$')
ax100c.set_title('LWP > 350')
# plt.legend()
#cb1 = plt.colorbar(lines1,ax=ax100)
plt.show()








fig = plt.figure(figsize=(12,12))
ax100 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
con1 = ax100.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')

for ff in range(len(highestWE)):
    if ff > 0:
        # si = lowWaves[ff][0]
        # ei = lowWaves[ff][1]
        # si = medLowWaves[ff][0]
        # ei = medLowWaves[ff][1]
        si = highestWE[ff][0]-1
        ei = highestWE[ff][1]-0
        # si = bigWaves[ff][0]
        # ei = bigWaves[ff][1]
        indexedWaves = maxMonthBeforeWE[si-1:ei+1]
        indexedTime = time[si-1:ei+1]
        indexedOrdinalTime = [i.toordinal() for i in indexedTime]
        indexedTimeMonth = [i.month for i in indexedTime]
        mode1prop = mode1Phase[si-1:ei+1]
        mode2prop = mode2Phase[si-1:ei+1]
        mode1propWrap = np.unwrap((mode1prop + np.pi) % (2 * np.pi) - np.pi)
        mode2propWrap = np.unwrap((mode2prop + np.pi) % (2 * np.pi) - np.pi)
        x1 = np.sort(mode1propWrap*180/np.pi)
        y1 = mode2propWrap*180/np.pi

        #linewidthV = np.round(np.max(indexedWaves))
        # if y1[-1] < -20:
        #     y1 = y1 + 360
        #newx = np.arange(x1[0],x1[-1],1)
        #f = interp1d(x1,y1,kind='linear')
        #newy = f(newx)
        #if x1[0] < -100:
        #    lines1 = ax100.plot(newx+360, newy, 'k-', linewidth=2)
        #else:
        lines1 = ax100.plot(x1,y1,'-',label='{}/{}/{} - {}'.format(time[si].month,time[si].day,time[si].year,np.max(indexedWaves)))

        # lines1 = ax100.plot(x1,y1,'-',linewidth=linewidthV,label='{}/{}/{} - {}'.format(time[si].month,time[si].day,time[si].year,np.max(indexedWaves)))
        #    lines1 = ax100.plot(newx,newy,'k-',linewidth=2)
        #scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedWaves,cmap='rainbow',zorder=3,vmin=3.5,vmax=6.5)
        #add_arrow(lines1[0], color='k', size=25)
        # scatter1 = ax100.scatter(mode1propWrap*180/np.pi,mode2propWrap*180/np.pi,35,indexedTimeMonth,cmap='twilight_shifted',zorder=3,vmin=1,vmax=12)



ax100.set_xlim([-180, 340])
ax100.set_ylim([-160, 360])
# cb = fig.colorbar(scatter1,ax=ax100)
ax100.set_xlabel(r'Offshore Propagation $\Longrightarrow$')
ax100.set_ylabel(r'Onshore Propagation $\Longrightarrow$')
ax100.set_title('Non-dimensional Fall Velocity > 9')
plt.legend()
plt.show()

