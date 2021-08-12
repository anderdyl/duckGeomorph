

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


dbfile = open('sandbarsSouthernTransectRandomForestInputData13.pickle','rb')
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
#sortedbmus = data['sortedbmus']


def listToArray(data):

    param05 = [item for sublist in data['param05'] for item in sublist]
    param05Removed = [ele for ele in param05 if ele != []]
    param05Array = np.concatenate(param05Removed).ravel()

    param25 = [item for sublist in data['param25'] for item in sublist]
    param25Removed = [ele for ele in param25 if ele != []]
    param25Array = np.concatenate(param25Removed).ravel()

    param50 = [item for sublist in data['param50'] for item in sublist]
    param50Removed = [ele for ele in param50 if ele != []]
    param50Array = np.concatenate(param50Removed).ravel()

    param75 = [item for sublist in data['param75'] for item in sublist]
    param75Removed = [ele for ele in param75 if ele != []]
    param75Array = np.concatenate(param75Removed).ravel()

    param90 = [item for sublist in data['param90'] for item in sublist]
    param90Removed = [ele for ele in param90 if ele != []]
    param90Array = np.concatenate(param90Removed).ravel()

    param95 = [item for sublist in data['param95'] for item in sublist]
    param95Removed = [ele for ele in param95 if ele != []]
    param95Array = np.concatenate(param95Removed).ravel()

    output = dict()
    output['param05Array'] = param05Array
    output['param25Array'] = param25Array
    output['param50Array'] = param50Array
    output['param75Array'] = param75Array
    output['param90Array'] = param90Array
    output['param95Array'] = param95Array
    return output



hsArrays = listToArray(hsSplit)
hsTran05 = hsArrays['param05Array']
hsTran25 = hsArrays['param25Array']
hsTran50 = hsArrays['param50Array']
hsTran75 = hsArrays['param75Array']
hsTran90 = hsArrays['param90Array']
hsTran95 = hsArrays['param95Array']

tpArrays = listToArray(tpSplit)
tpTran05 = tpArrays['param05Array']
tpTran25 = tpArrays['param25Array']
tpTran50 = tpArrays['param50Array']
tpTran75 = tpArrays['param75Array']
tpTran90 = tpArrays['param90Array']
tpTran95 = tpArrays['param95Array']

dmArrays = listToArray(dmSplit)
dmTran05 = dmArrays['param05Array']
dmTran25 = dmArrays['param25Array']
dmTran50 = dmArrays['param50Array']
dmTran75 = dmArrays['param75Array']
dmTran90 = dmArrays['param90Array']
dmTran95 = dmArrays['param95Array']

lwpArrays = listToArray(lwpSplit)
lwpTran05 = lwpArrays['param05Array']
lwpTran25 = lwpArrays['param25Array']
lwpTran50 = lwpArrays['param50Array']
lwpTran75 = lwpArrays['param75Array']
lwpTran90 = lwpArrays['param90Array']
lwpTran95 = lwpArrays['param95Array']

weArrays = listToArray(weSplit)
weTran05 = weArrays['param05Array']
weTran25 = weArrays['param25Array']
weTran50 = weArrays['param50Array']
weTran75 = weArrays['param75Array']
weTran90 = weArrays['param90Array']
weTran95 = weArrays['param95Array']

hoverlArrays = listToArray(hoverlSplit)
hoverlTran05 = hoverlArrays['param05Array']
hoverlTran25 = hoverlArrays['param25Array']
hoverlTran50 = hoverlArrays['param50Array']
hoverlTran75 = hoverlArrays['param75Array']
hoverlTran90 = hoverlArrays['param90Array']
hoverlTran95 = hoverlArrays['param95Array']

irArrays = listToArray(irSplit)
irTran05 = irArrays['param05Array']
irTran25 = irArrays['param25Array']
irTran50 = irArrays['param50Array']
irTran75 = irArrays['param75Array']
irTran90 = irArrays['param90Array']
irTran95 = irArrays['param95Array']

fvArrays = listToArray(fvSplit)
fvTran05 = fvArrays['param05Array']
fvTran25 = fvArrays['param25Array']
fvTran50 = fvArrays['param50Array']
fvTran75 = fvArrays['param75Array']
fvTran90 = fvArrays['param90Array']
fvTran95 = fvArrays['param95Array']

windSpeedArrays = listToArray(windSpeedSplit)
windSpeedTran05 = windSpeedArrays['param05Array']
windSpeedTran25 = windSpeedArrays['param25Array']
windSpeedTran50 = windSpeedArrays['param50Array']
windSpeedTran75 = windSpeedArrays['param75Array']
windSpeedTran90 = windSpeedArrays['param90Array']
windSpeedTran95 = windSpeedArrays['param95Array']

windDirectionArrays = listToArray(windDirectionSplit)
windDirectionTran05 = windDirectionArrays['param05Array']
windDirectionTran25 = windDirectionArrays['param25Array']
windDirectionTran50 = windDirectionArrays['param50Array']
windDirectionTran75 = windDirectionArrays['param75Array']
windDirectionTran90 = windDirectionArrays['param90Array']
windDirectionTran95 = windDirectionArrays['param95Array']

resWaterLevelArrays = listToArray(resWaterLevelSplit)
resWaterLevelTran05 = resWaterLevelArrays['param05Array']
resWaterLevelTran25 = resWaterLevelArrays['param25Array']
resWaterLevelTran50 = resWaterLevelArrays['param50Array']
resWaterLevelTran75 = resWaterLevelArrays['param75Array']
resWaterLevelTran90 = resWaterLevelArrays['param90Array']
resWaterLevelTran95 = resWaterLevelArrays['param95Array']



# goodWindInd = np.where((windSpeedTran90 > -50))
#
# windDirectionTran05 = windDirectionTran05[goodWindInd]
# windDirectionTran25 = windDirectionTran25[goodWindInd]
# windDirectionTran50 = windDirectionTran50[goodWindInd]
# windDirectionTran75 = windDirectionTran75[goodWindInd]
# windDirectionTran90 = windDirectionTran90[goodWindInd]
# windDirectionTran95 = windDirectionTran95[goodWindInd]
#
# windSpeedTran05 = windSpeedTran05[goodWindInd]
# windSpeedTran25 = windSpeedTran25[goodWindInd]
# windSpeedTran50 = windSpeedTran50[goodWindInd]
# windSpeedTran75 = windSpeedTran75[goodWindInd]
# windSpeedTran90 = windSpeedTran90[goodWindInd]
# windSpeedTran95 = windSpeedTran95[goodWindInd]
#
# resWaterLevelTran05 = resWaterLevelTran05[goodWindInd]
# resWaterLevelTran25 = resWaterLevelTran25[goodWindInd]
# resWaterLevelTran50 = resWaterLevelTran50[goodWindInd]
# resWaterLevelTran75 = resWaterLevelTran75[goodWindInd]
# resWaterLevelTran90 = resWaterLevelTran90[goodWindInd]
# resWaterLevelTran95 = resWaterLevelTran95[goodWindInd]
#
# hsTran05 = hsTran05[goodWindInd]
# hsTran25 = hsTran25[goodWindInd]
# hsTran50 = hsTran50[goodWindInd]
# hsTran75 = hsTran75[goodWindInd]
# hsTran90 = hsTran90[goodWindInd]
# hsTran95 = hsTran95[goodWindInd]
#
# tpTran05 = tpTran05[goodWindInd]
# tpTran25 = tpTran25[goodWindInd]
# tpTran50 = tpTran50[goodWindInd]
# tpTran75 = tpTran75[goodWindInd]
# tpTran90 = tpTran90[goodWindInd]
# tpTran95 = tpTran95[goodWindInd]
#
# dmTran05 = dmTran05[goodWindInd]
# dmTran25 = dmTran25[goodWindInd]
# dmTran50 = dmTran50[goodWindInd]
# dmTran75 = dmTran75[goodWindInd]
# dmTran90 = dmTran90[goodWindInd]
# dmTran95 = dmTran95[goodWindInd]
#
# lwpTran05 = lwpTran05[goodWindInd]
# lwpTran25 = lwpTran25[goodWindInd]
# lwpTran50 = lwpTran50[goodWindInd]
# lwpTran75 = lwpTran75[goodWindInd]
# lwpTran90 = lwpTran90[goodWindInd]
# lwpTran95 = lwpTran95[goodWindInd]
#
# weTran05 = weTran05[goodWindInd]
# weTran25 = weTran25[goodWindInd]
# weTran50 = weTran50[goodWindInd]
# weTran75 = weTran75[goodWindInd]
# weTran90 = weTran90[goodWindInd]
# weTran95 = weTran95[goodWindInd]
#
# irTran05 = irTran05[goodWindInd]
# irTran25 = irTran25[goodWindInd]
# irTran50 = irTran50[goodWindInd]
# irTran75 = irTran75[goodWindInd]
# irTran90 = irTran90[goodWindInd]
# irTran95 = irTran95[goodWindInd]
#
# hoverlTran05 = hoverlTran05[goodWindInd]
# hoverlTran25 = hoverlTran25[goodWindInd]
# hoverlTran50 = hoverlTran50[goodWindInd]
# hoverlTran75 = hoverlTran75[goodWindInd]
# hoverlTran90 = hoverlTran90[goodWindInd]
# hoverlTran95 = hoverlTran95[goodWindInd]



# import matplotlib.cm as cm
# import matplotlib.colors as mcolors
# from scipy.stats.kde import gaussian_kde
#
# def getDistributionPercentages(data,numClusters):
#
#     #minData = np.nanmin(data)
#     #maxData = np.nanmax(data)
#     #dist_space = np.linspace(minData,maxData,numBins)
#     param05 = np.zeros((numClusters*numClusters,))
#     param25 = np.zeros((numClusters*numClusters,))
#     param50 = np.zeros((numClusters*numClusters,))
#     param75 = np.zeros((numClusters*numClusters,))
#     param90 = np.zeros((numClusters*numClusters,))
#     param95 = np.zeros((numClusters*numClusters,))
#
#     counter = 0
#     for xx in range(numClusters):
#         for yy in range(numClusters):
#
#             transitionData = data[xx][yy]
#             if len(data)>0:
#                 #kdes = gaussian_kde(transitionData)
#                 #param50[counter] = np.nanmean(transitionData)
#                 param05[counter] = np.nanpercentile(transitionData,5)
#                 param25[counter] = np.nanpercentile(transitionData,25)
#                 param50[counter] = np.nanpercentile(transitionData,50)
#                 param75[counter] = np.nanpercentile(transitionData,75)
#                 param90[counter] = np.nanpercentile(transitionData,90)
#                 param95[counter] = np.nanpercentile(transitionData,95)
#
#                 #kde_ci = np.array(np.percentile(kdes, (5, 95), axis=0))
#             counter = counter+1
#
#     output = dict()
#     output['param05'] = param05
#     output['param25'] = param25
#     output['param50'] = param50
#     output['param75'] = param75
#     output['param90'] = param90
#     output['param95'] = param95
#
#     return output



# transHs = getDistributionPercentages(hs,8)
# transTp = getDistributionPercentages(tp,8)
# transDm = getDistributionPercentages(dm,8)
# transWE = getDistributionPercentages(we,8)
# transHoverL = getDistributionPercentages(hoverl,8)
# transIr = getDistributionPercentages(ir,8)
# transLWP = getDistributionPercentages(lwp,8)
# transWindSpeed = getDistributionPercentages(windSpeed,8)
# transSusWindSpeed = getDistributionPercentages(sustWindSpeed,8)
# #transWindDirection = getDistributionPercentages(windDirection,8)
# transResWaterLevel = getDistributionPercentages(resWaterLevel,8)
# transWaterLevel = getDistributionPercentages(waterLevel,8)





# avgHsList = [item for sublist in avgHs for item in sublist]
# avgHsRemoved = [ele for ele in avgHsList if ele != []]
# avgHsArray = np.concatenate(avgHsRemoved).ravel()



#
#
#
#
# avgWEList = [item for sublist in avgWE for item in sublist]
# avgWERemoved = [ele for ele in avgWEList if ele != []]
# avgWEArray = np.concatenate(avgWERemoved).ravel()
#
avgDayList = [item for sublist in days for item in sublist]
avgDayRemoved = [ele for ele in avgDayList if ele != []]
avgDayArray = np.concatenate(avgDayRemoved).ravel()
#

#
# avgTpList = [item for sublist in avgTp for item in sublist]
# avgTpRemoved = [ele for ele in avgTpList if ele != []]
# avgTpArray = np.concatenate(avgTpRemoved).ravel()
#
# avgLWPList = [item for sublist in avgLWP for item in sublist]
# avgLWPRemoved = [ele for ele in avgLWPList if ele != []]
# avgLWPArray = np.concatenate(avgLWPRemoved).ravel()
#
prevBmuList = [item for sublist in prevBmu for item in sublist]
prevBmuRemoved = [ele for ele in prevBmuList if ele != []]
prevBmuArray = np.concatenate(prevBmuRemoved).ravel()

currentBmuList = [item for sublist in currentBmu for item in sublist]
currentBmuRemoved = [ele for ele in currentBmuList if ele != []]
currentBmuArray = np.concatenate(currentBmuRemoved).ravel()
#
# stormCumuWEList = [item for sublist in stormCumuWE for item in sublist]
# stormCumuWERemoved = [ele for ele in stormCumuWEList if ele != []]
# stormCumuWEArray = np.concatenate(stormCumuWERemoved).ravel()
#
# stormCumuLWPList = [item for sublist in stormCumuLWP for item in sublist]
# stormCumuLWPRemoved = [ele for ele in stormCumuLWPList if ele != []]
# stormCumuLWPArray = np.concatenate(stormCumuLWPRemoved).ravel()
#
# avgIrList = [item for sublist in avgIr for item in sublist]
# avgIrRemoved = [ele for ele in avgIrList if ele != []]
# avgIrArray = np.concatenate(avgIrRemoved).ravel()
#
# transitionType = np.empty((len(currentBmuArray),))
# for ii in range(len(currentBmuArray)):
#     c = currentBmuArray[ii]
#     p = prevBmuArray[ii]
#     if c == 0:
#         if p == 0:
#             transitionType[ii] = 0
#         elif p > 0 and p <= 4:
#             transitionType[ii] = 1
#         elif p > 4:
#             transitionType[ii] = 2
#     elif c == 1:
#         if p == 1:
#             transitionType[ii] = 0
#         elif p > 1 and p <= 5:
#             transitionType[ii] = 1
#         elif p > 5:
#             transitionType[ii] = 2
#         elif p < 1:
#             transitionType[ii] = 2
#     elif c == 2:
#         if p == 2:
#             transitionType[ii] = 0
#         elif p > 2 and p < 6:
#             transitionType[ii] = 1
#         elif p >= 6:
#             transitionType[ii] = 2
#         elif p < 2:
#             transitionType[ii] = 2
#     elif c == 3:
#         if p == 3:
#             transitionType[ii] = 0
#         elif p > 3 and p < 7:
#             transitionType[ii] = 1
#         elif p >= 7:
#             transitionType[ii] = 2
#         elif p < 3:
#             transitionType[ii] = 2
#     elif c == 4:
#         if p == 4:
#             transitionType[ii] = 0
#         elif p > 4 and p < 8:
#             transitionType[ii] = 1
#         elif p == 0:
#             transitionType[ii] = 1
#         elif p < 4 and p > 0:
#             transitionType[ii] = 2
#     elif c == 5:
#         if p == 5:
#             transitionType[ii] = 0
#         elif p > 5 and p < 9:
#             transitionType[ii] = 1
#         elif p >= 0 and p <= 1:
#             transitionType[ii] = 1
#         elif p > 1 and p < 5:
#             transitionType[ii] = 2
#     elif c == 6:
#         if p == 6:
#             transitionType[ii] = 0
#         elif p > 6 and p < 9:
#             transitionType[ii] = 1
#         elif p >= 0 and p <= 2:
#             transitionType[ii] = 1
#         elif p > 2 and p < 6:
#             transitionType[ii] = 2
#     elif c == 7:
#         if p == 7:
#             transitionType[ii] = 0
#         elif p >= 0 and p <= 3:
#             transitionType[ii] = 1
#         elif p > 3 and p < 7:
#             transitionType[ii] = 2


transitionType = np.empty((len(currentBmuArray),))
for ii in range(len(currentBmuArray)):
    c = currentBmuArray[ii]
    p = prevBmuArray[ii]
    if c == 0:
        if p == 0:
            transitionType[ii] = 0
        elif p > 0 and p <= 6:
            transitionType[ii] = 2
        elif p > 6:
            transitionType[ii] = 1
    elif c == 1:
        if p == 1:
            transitionType[ii] = 0
        elif p > 1 and p <= 6:
            transitionType[ii] = 2
        elif p > 6:
            transitionType[ii] = 1
        elif p == 0:
            transitionType[ii] = 1
    elif c == 2:
        if p == 2:
            transitionType[ii] = 0
        elif p > 2 and p <= 8:
            transitionType[ii] = 2
        elif p > 8:
            transitionType[ii] = 1
        elif p < 2:
            transitionType[ii] = 1
    elif c == 3:
        if p == 3:
            transitionType[ii] = 0
        elif p > 3 and p <= 9:
            transitionType[ii] = 2
        elif p > 9:
            transitionType[ii] = 1
        elif p < 3:
            transitionType[ii] = 1
    elif c == 4:
        if p == 4:
            transitionType[ii] = 0
        elif p > 4 and p <= 10:
            transitionType[ii] = 2
        elif p > 10:
            transitionType[ii] = 1
        elif p < 4:
            transitionType[ii] = 1
    elif c == 5:
        if p == 5:
            transitionType[ii] = 0
        elif p > 5 and p <= 11:
            transitionType[ii] = 2
        elif p > 11:
            transitionType[ii] = 1
        elif p >= 0 and p <= 4:
            transitionType[ii] = 1


    elif c == 6:
        if p == 6:
            transitionType[ii] = 0
        elif p > 6 and p <= 12:
            transitionType[ii] = 2
        elif p >= 0 and p <= 5:
            transitionType[ii] = 1


    elif c == 7:
        if p == 7:
            transitionType[ii] = 0
        elif p > 7 and p <= 13:
            transitionType[ii] = 2
        elif p == 0:
            transitionType[ii] = 2
        elif p > 0 and p < 7:
            transitionType[ii] = 1


    elif c == 8:
        if p == 8:
            transitionType[ii] = 0
        elif p > 8 and p <= 14:
            transitionType[ii] = 2
        elif p >= 0 and p <= 1:
            transitionType[ii] = 2
        elif p > 1 and p < 8:
            transitionType[ii] = 1


    elif c == 9:
        if p == 9:
            transitionType[ii] = 0
        elif p > 9 and p <= 15:
            transitionType[ii] = 2
        elif p >= 0 and p <= 2:
            transitionType[ii] = 2
        elif p > 2 and p < 9:
            transitionType[ii] = 1


    elif c == 10:
        if p == 10:
            transitionType[ii] = 0
        elif p > 10 and p <= 15:
            transitionType[ii] = 2
        elif p >= 0 and p <= 3:
            transitionType[ii] = 2
        elif p > 3 and p < 10:
            transitionType[ii] = 1

    elif c == 11:
        if p == 11:
            transitionType[ii] = 0
        elif p > 11 and p <= 15:
            transitionType[ii] = 2
        elif p >= 0 and p <= 4:
            transitionType[ii] = 2
        elif p > 4 and p < 11:
            transitionType[ii] = 1

    elif c == 12:
        if p == 12:
            transitionType[ii] = 0
        elif p > 12 and p <= 15:
            transitionType[ii] = 2
        elif p >= 0 and p <= 5:
            transitionType[ii] = 2
        elif p > 5 and p < 12:
            transitionType[ii] = 1

    elif c == 13:
        if p == 13:
            transitionType[ii] = 0
        elif p > 13 and p <= 15:
            transitionType[ii] = 2
        elif p >= 0 and p <= 6:
            transitionType[ii] = 2
        elif p > 6 and p < 13:
            transitionType[ii] = 1


transitionArray = currentBmuArray-prevBmuArray
posInd = np.where((transitionArray>1))
transitionArray[posInd] = 1
negInd = np.where((transitionArray<0))
transitionArray[negInd] = 1

surveys = np.arange(0,len(transitionType),1)

import pandas as pd
#dataList = [['currentState',currentBmuArray], ['previousState',prevBmuArray],['averageHs',avgHsArray]]
# dataDict = {'currentState':currentBmuArray, 'previousState':prevBmuArray, 'cumulativeWE':avgWEArray,
#             'cumulativeLWP':avgLWPArray, 'averageIr':avgIrArray, 'numberOfDays':avgDayArray,
#             'stormCumuWE':stormCumuWEArray, 'stormCumuLWP':stormCumuLWPArray}
            #'stormCumuWE':stormCumuWEArray, 'stormCumuLWP':stormCumuLWPArray}
# dataDict = {'transition':transitionArray, 'previousState':prevBmuArray, 'averageWE':avgWEArray,
#             'averageLWP':avgLWPArray, 'numberOfDays':avgDayArray}
# dataDict = {'currentState':currentBmuArray, 'previousState':prevBmuArray, '90WE':weTran90,
#             '90HS':hsTran90,'50Ir':irTran50}
dataDict = {'transition':transitionType, 'LWP90':lwpTran90, 'LWP50':lwpTran50, 'LWP25':lwpTran25,
            'WE90':weTran90, 'WE50':weTran50, 'WE25':weTran25,
            'HS90':hsTran90, 'HS50':hsTran50, 'HS25':hsTran25,
            'Ir90':irTran90, 'Ir50':irTran50, 'Ir25':irTran25,
            'FV90': fvTran90, 'FV50': fvTran50, 'FV25': fvTran25,
            'H/L90':hoverlTran90, 'H/L50':hoverlTran50, 'H/L25':hoverlTran25,
            'Residual90':resWaterLevelTran90, 'Residual50':resWaterLevelTran50, 'Residual25':resWaterLevelTran25,
            'windDirection50':windDirectionTran50,
            'windSpeed90':windSpeedTran90, 'windSpeed50':windSpeedTran50, 'windSpeed25':windSpeedTran25}

# dataDict = {'previousState':prevBmuArray,
#             'HS90':hsTran90,
#             'transition':transitionType}#,
            # 'LWP90':lwpTran90, 'LWP50':lwpTran50,
            # 'WE90':weTran90, 'WE50':weTran50,
            # 'HS90':hsTran90, 'HS50':hsTran50,
            # 'Ir90':irTran90, 'Ir50':irTran50,
            # 'FV90': fvTran90, 'FV50': fvTran50,
            # 'H/L90':hoverlTran90, 'H/L50':hoverlTran50,
            # 'Residual90':resWaterLevelTran90, 'Residual50':resWaterLevelTran50,
            # 'windDirection50':windDirectionTran50,
            # 'windSpeed90':windSpeedTran90, 'windSpeed50':windSpeedTran50}
# dataDict = {'transition':transitionType,
#             'prevState':prevBmuArray,
#             'daysBetween':avgDayArray,
#             'LWP90':lwpTran90, 'LWP50':lwpTran50,
#             'WE90':weTran90, 'WE50':weTran50,
#             'HS90':hsTran90, 'HS50':hsTran50,
#             'Ir90':irTran90, 'Ir50':irTran50,
#             'FV90': fvTran90, 'FV50': fvTran50,
#             'H/L90':hoverlTran90, 'H/L50':hoverlTran50,
#             'Residual90':resWaterLevelTran90, 'Residual50':resWaterLevelTran50,
#             'windDirection50':windDirectionTran50,
#             'windSpeed90':windSpeedTran90, 'windSpeed50':windSpeedTran50}
# dataDict = {'currentState':currentBmuArray,
#             'LWP90':lwpTran90, 'LWP50':lwpTran50,
#             'WE90':weTran90, 'WE50':weTran50,
#             'HS90':hsTran90, 'HS50':hsTran50,
#             'Ir90':irTran90, 'Ir50':irTran50,
#             'FV90': fvTran90, 'FV50': fvTran50,
#             'H/L90':hoverlTran90, 'H/L50':hoverlTran50,
#             'Residual90':resWaterLevelTran90, 'Residual50':resWaterLevelTran50,
#             'windDirection50':windDirectionTran50,
#             'windSpeed90':windSpeedTran90, 'windSpeed50':windSpeedTran50}

#             'stormCumuWE':stormCumuWEArray, 'stormCumuLWP':stormCumuLWPArray}
#'averageHs':avgHsArray, 'averageTp':avgTpArray,

df = pd.DataFrame(dataDict)

df = df[df['Residual90'].notna()]
# Need to clean the dataframe...
df = df[df.windSpeed50 != -99]


#

# df['previousState'] = df['previousState'].astype('category')
# df = pd.get_dummies(df,'previousState')
#



# Labels are the values we want to predict
labels = np.array(df['transition'])
# Remove the labels from the features
# axis 1 refers to the columns
features = df.drop('transition', axis = 1)
# Saving feature names for later use
feature_list = list(features.columns)
# Convert to numpy array
features = np.array(features)

# Using Skicit-learn to split data into training and testing sets
from sklearn.model_selection import train_test_split
# Split the data into training and testing sets
train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size = 0.25, random_state = 5)

print('Training Features Shape:', train_features.shape)
print('Training Labels Shape:', train_labels.shape)
print('Testing Features Shape:', test_features.shape)
print('Testing Labels Shape:', test_labels.shape)

# # The baseline predictions are the historical averages
# baseline_preds = test_features[:, feature_list.index('previousState')]
# # Baseline errors, and display average baseline error
# baseline_errors = abs(baseline_preds - test_labels)
# print('Average baseline error: ', round(np.mean(baseline_errors), 2))

# Import the model we are using
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
# Instantiate model with 1000 decision trees
rf = RandomForestClassifier(n_estimators = 1000, random_state = 5)
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix


# Train the model on training data
rf.fit(train_features, train_labels)

# predictions = rf.predict(train_features)
# cm = confusion_matrix(train_labels, predictions)
# scoreTest = rf.score(train_features,train_labels)
# scoreTest2 = accuracy_score(train_labels,predictions)

predictions = rf.predict(test_features)
cm = confusion_matrix(test_labels, predictions)
scoreTest = rf.score(test_features,test_labels)
scoreTest2 = accuracy_score(test_labels,predictions)

print(scoreTest)
print(scoreTest2)
# # Use the forest's predict method on the test data
# Calculate the absolute errors
# errors = abs(predictions - test_labels)# Print out the mean absolute error (mae)
# print('Mean Absolute Error:', round(np.mean(errors), 2), 'degrees.')
plt.figure()
plt.bar(x=feature_list, height=rf.feature_importances_)
plt.xticks(rotation=30, ha='right')

print(pd.crosstab(test_labels, predictions))







import h2o
from h2o.estimators import H2ORandomForestEstimator

h2o.init()
# dataDict = {'previousState':prevBmuArray,
#             'HS90':hsTran90,
#             'transition':transitionType}#,

dataDict = {'FV90': fvTran90,'LWP90':lwpTran90,'WE90':weTran90,'days':avgDayArray,'H/L90':hoverlTran90,'Ir90':irTran90,
            'Residual90':resWaterLevelTran90,'windSpeed90':windSpeedTran90,
            'transition':transitionType}#,
            #  'LWP50':lwpTran50,,'previousState':prevBmuArray
            # 'WE90':weTran90, 'WE50':weTran50,
            # 'HS90':hsTran90, 'HS50':hsTran50,
            #  'Ir50':irTran50,
            # , 'FV50': fvTran50,
            #  'H/L50':hoverlTran50,
            #  'Residual50':resWaterLevelTran50,
            # 'windDirection50':windDirectionTran50,
            #  'windSpeed50':windSpeedTran50}

df = pd.DataFrame(dataDict)

df = df[df['Residual90'].notna()]
# Need to clean the dataframe...
df = df[df.windSpeed90 != -99]
# df['previousState'] = df['previousState'].astype('category')
# df = pd.get_dummies(df,'previousState')
hf = h2o.H2OFrame(df) #, column_types={'previousState':'enum', 'HS90':'numeric', 'transition':'enum'})

#y = 'transition'
#x = list(hf.columns).remove(y)

hf["transition"] = hf["transition"].asfactor()
# hf["previousState"] = hf["previousState"].asfactor()


y = 'transition'
x = hf.col_names
x.remove(y)
print("Response = " + y)
print("Pridictors = " + str(x))

hf['transition'].levels()
train, valid, test = hf.split_frame(ratios=[.7, .2])
print(df.shape)
print(train.shape)
print(valid.shape)
print(test.shape)


model = H2ORandomForestEstimator(ntrees=500, max_depth=20, nfolds=10)

model.train(x=x, y=y, training_frame=train)

performance = model.model_performance(test_data=test)

print(performance)
