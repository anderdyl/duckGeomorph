#   Lets try to average everything in 900-1100 range to see if we add data/smooth alongshore anomalies
#
from getdatatestbed import getDataFRF
import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from glob import glob
import os
from scipy.interpolate import interp1d
import sandBarTool.morphLib as mL

#from sandBarTool import morphLib
#from downloads import pyeemd



geomorphdir = '/media/dylananderson/Elements/filteredFRF_Geomorph/'

files = os.listdir(geomorphdir)

files.sort()

files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]

def getBathy(file, lower, upper):
    bathy = Dataset(file)

    xs_bathy = bathy.variables['xFRF'][:]
    ys_bathy = bathy.variables['yFRF'][:]
    zs_bathy = bathy.variables['elevation'][:]
    ts_bathy = bathy.variables['time'][:]
    pr_bathy = bathy.variables['profileNumber'][:]

    zs_bathy = np.ma.masked_where((pr_bathy > upper), zs_bathy)
    ys_bathy = np.ma.masked_where((pr_bathy > upper), ys_bathy)
    xs_bathy = np.ma.masked_where((pr_bathy > upper), xs_bathy)
    pr_bathy = np.ma.masked_where((pr_bathy > upper), pr_bathy)
    ts_bathy = np.ma.masked_where((pr_bathy > upper), ts_bathy)

    zs_bathy = np.ma.masked_where((pr_bathy < lower), zs_bathy)
    ys_bathy = np.ma.masked_where((pr_bathy < lower), ys_bathy)
    xs_bathy = np.ma.masked_where((pr_bathy < lower), xs_bathy)
    pr_bathy = np.ma.masked_where((pr_bathy < lower), pr_bathy)
    ts_bathy = np.ma.masked_where((pr_bathy < lower), ts_bathy)

    output = dict()
    output['x'] = xs_bathy
    output['y'] = ys_bathy
    output['z'] = zs_bathy
    output['pr'] = pr_bathy
    output['t'] = ts_bathy

    return output


def interpBathy(xSub, zSub, x):

    f = interp1d(xSub, zSub, kind='linear', bounds_error=False)
    newz = f(x)
    newBathy = dict()
    newBathy['x'] = x
    newBathy['z'] = newz
    return newBathy

def interpBathy2(xSub,zSub,x):
    from Loess import Loess

    loess = Loess(xSub, zSub)
    zloess = np.nan * np.ones(len(x),)
    count = 0
    for xxx in x:
        zloess[count] = loess.estimate(xxx, window=5, use_matrix=False, degree=1)
        count = count+1
    newBathy = dict()
    newBathy['x'] = x
    newBathy['z'] = zloess
    return newBathy


subset = files[0:972].copy()
#from palettable.colorbrewer.sequential import Blues_8
#ax.imshow(data, cmap=Blues_8.mpl_colormap

colormap = plt.cm.gist_ncar

import scipy.stats as spstats

labels = []
# xinterp = np.arange(108, 608, 2.5)
xinterp = np.arange(98, 548, 5)
deeperThan = 523

plt.style.use('dark_background')


bathy = dict()
#alllines = np.empty((len(xinterp),))
count = 0
count2 = 0
count3 = 0
count4 = 0
count5 = 0
worstcount = 0
#fig = plt.figure(figsize=(10,10))
for i in range(len(subset)):
    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=1070, upper=1100)
    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=750, upper=950)
    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-10, upper=20)
    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-2, upper=10)
    data = getBathy(os.path.join(geomorphdir, subset[i]), lower=990, upper=1100)

    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-2, upper=10)

    temp = subset[i].split('_')

    surveydate = DT.datetime.strptime(temp[1], '%Y%m%d')
    elevs = data['z']
    cross = data['x']
    crossind = np.argsort(data['x'])
    crossS = cross[crossind]
    elevsS = elevs[crossind]

    xSub = np.ma.MaskedArray.filled(crossS, np.nan)
    zSub = np.ma.MaskedArray.filled(elevsS, np.nan)

    #realValues = ~np.isnan(xSub)

    binnedz, bin_edges, binnumber = spstats.binned_statistic(xSub,zSub, statistic='mean',bins=np.arange(98,566,4))
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width / 2
    xMSub = bin_centers
    zMSub = binnedz
    xSubNewB = xMSub[~np.isnan(zMSub)]
    zSubNewB = zMSub[~np.isnan(zMSub)]
    realValues = ~np.isnan(xSubNewB)
    zSubNew = zSubNewB[~np.isnan(xSubNewB)]
    xSubNew = xSubNewB[~np.isnan(xSubNewB)]

    # realValues = ~np.isnan(xSub)
    # xSubNew = xSub[~np.isnan(xSub)]
    # zSubNew = zSub[~np.isnan(xSub)]



    if any(realValues):

        if np.nanmax(xSubNew) > deeperThan:

            if np.nanmin(xSubNew) < 125:

                if len(xSubNew) > 10:
                    #newdata = interpBathy(xSubNew, zSubNew, xinterp)
                    #newdata = interpBathy2(newdata['x'], newdata['z'], xinterp)
                    newdata = interpBathy2(xSubNew, zSubNew, xinterp)
                    nanValues = np.isnan(newdata['z'])

                    if any(nanValues):
                        #print('Found a transect with nans {}'.format(i))
                        print('Trying to extend the lower half of the profile with a linear fit: {}'.format(surveydate))
                        if count2 == 0:
                            badlines = newdata['z']
                            count2 = count2+1
                            badtime = surveydate
                        else:
                            badlines = np.vstack((badlines, newdata['z']))
                            count2 = count2+1
                            badtime = np.append(badtime, surveydate)

                        #print('Threw out a survey with data due to nans: {}'.format(surveydate))
                        #plt.plot(xSub,zSub)


                        xindex = np.where((xSubNew>=deeperThan))

                        if len(xindex[0]) > 3:

                            xbot = xSubNew[xindex]
                            zbot = zSubNew[xindex]
                            mb = np.polyfit(xbot, zbot, 1)
                            f = np.poly1d(mb)
                            maxXind = np.where((xinterp > np.nanmax(xSubNew)))
                            newX = xinterp[maxXind]
                            newZ = f(newX)
                            xSubNew2 = np.append(xSubNew, newX)
                            zSubNew2 = np.append(zSubNew, newZ)
                            moredata = interpBathy(xSubNew2, zSubNew2, xinterp)
                            #moredata = interpBathy(moredata['x'], moredata['z'], xinterp)
                            #moredata = interpBathy(xSubNew2, moredata['z'], xinterp)

                            del xSubNew2, zSubNew2,newZ,newX,f,mb,zbot,xbot

                            nanValues2 = np.isnan(moredata['z'])

                            if any(nanValues2):
                                print('WHAT HAPPENED AT TRANSECT {}'.format(i))
                                #plt.plot(moredata['x'], moredata['z'], 'r-')
                            else:
                                if count == 0:
                                    alllines = moredata['z']
                                    time = surveydate
                                    count = count + 1
                                else:
                                    alllines = np.vstack((alllines, moredata['z']))
                                    time = np.append(time, surveydate)
                                    count = count + 1


                                if count3 == 0:
                                    extrapolatedlines = newdata['z']
                                    count3 = count3 + 1
                                    extrapolatedtime = surveydate
                                else:
                                    extrapolatedlines = np.vstack((extrapolatedlines, newdata['z']))
                                    count3 = count3 + 1
                                    extrapolatedtime = np.append(extrapolatedtime, surveydate)
                                #plt.plot(moredata['x'], moredata['z'], 'k-')
                                del moredata



                        elif len(xindex[0]) == 3:
                            print('Could not do that: number of points being fit: {}'.format(len(xindex[0])))
                            worstcount = worstcount+1
                        else:
                            print('Could not do that: maximum X value in transect is: {}'.format(np.nanmax(xSub)))
                            worstcount = worstcount+1
                            #plt.plot(xSub, zSub, 'r-')


                    else:
                        if count == 0:
                            alllines = newdata['z']
                            time = surveydate
                            count = count+1
                        else:
                            alllines = np.vstack((alllines, newdata['z']))
                            time = np.append(time,surveydate)
                            count = count+1

                        fig = plt.figure(figsize=(10,10))
                        plt.plot(cross,elevs,'.')
                        plt.plot(xinterp,newdata['z'],'.-')
                        plt.xlim([0, 600])
                        plt.ylim([-7,2])
                        plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
                        plt.savefig('/home/dylananderson/projects/duckGeomorph/interpolatedProfiles/Survey_{}'.format(i))
                        plt.close()
                else:
                    print('Data is less than 10 points long at line {}'.format(i))
            else:
                print('No data onshore of 125 meters at line {}'.format(i))
                print('Most onshore point at {}'.format(np.nanmin(xSubNew)))
                newdata = interpBathy2(xSubNew, zSubNew, xinterp)

                # if count4 == 0:
                #     noneOnshore = newdata['z']
                #     noOnshoreTime = surveydate
                #     count4 = count4 + 1
                # else:
                #     noneOnshore = np.vstack((noneOnshore, newdata['z']))
                #     noOnshoreTime = np.append(noOnshoreTime, surveydate)
                #     count4 = count4 + 1
        else:
            print('No data deeper than 500 meters at line {}'.format(i))
            print('Most offshore point at {}'.format(np.nanmax(xSubNew)))
            #newdata = interpBathy2(xSubNew, zSubNew, xinterp)

            # if count5 == 0:
            #     noneOffshore = newdata['z']
            #     noOffshoreTime = surveydate
            #     count5 = count5 + 1
            # else:
            #     noneOffshore = np.vstack((noneOffshore, newdata['z']))
            #     noOffshoreTime = np.append(noOffshoreTime, surveydate)
            #     count4 = count5 + 1
    else:
        print('Survey with no data at this line {}'.format(i))

        #plt.plot(xSubNew,zSubNew)


#
# bathy = dict()
# #alllines = np.empty((len(xinterp),))
# count = 0
# count2 = 0
# worstcount = 0
# #fig = plt.figure(figsize=(10,10))
# for i in range(len(subset)):
#     #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=1080, upper=1100)
#     data = getBathy(os.path.join(geomorphdir, subset[i]), lower=900, upper=1100)
#
#     temp = subset[i].split('_')
#
#     surveydate = DT.datetime.strptime(temp[1], '%Y%m%d')
#
#     # fig = plt.figure(figsize=(10,10))
#     # plt.plot(data['x'],data['y'],'o')
#     # plt.plot([500,500],[-100,1200],'w--')
#     # plt.xlim([0, 800])
#     # plt.ylim([-100,1200])
#     # plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
#     # plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyExtents/Survey_{}'.format(i))
#     # plt.close()
#
#     elevs = data['z']
#     cross = data['x']
#     xSub = np.ma.MaskedArray.filled(cross, np.nan)
#     zSub = np.ma.MaskedArray.filled(elevs, np.nan)
#     realValues = ~np.isnan(xSub)
#
#     if any(realValues):
#
#         newdata = interpBathy(xSub, zSub, xinterp)
#         nanValues = np.isnan(newdata['z'])
#
#         if any(nanValues):
#             #print('Found a transect with nans {}'.format(i))
#             print('Trying to extend the lower half of the profile with a linear fit: {}'.format(surveydate))
#             if count2 == 0:
#                 badlines = newdata['z']
#                 count2 = count2+1
#                 badtime = surveydate
#             else:
#                 badlines = np.vstack((badlines, newdata['z']))
#                 count2 =count2+1
#                 badtime = np.append(badtime,surveydate)
#
#             #print('Threw out a survey with data due to nans: {}'.format(surveydate))
#             #plt.plot(xSub,zSub)
#
#
#             xindex = np.where((xSub>=650))
#
#             if len(xindex[0]) > 3:
#
#                 xbot = xSub[xindex]
#                 zbot = zSub[xindex]
#                 mb = np.polyfit(xbot, zbot, 1)
#                 f = np.poly1d(mb)
#                 maxXind = np.where((xinterp > np.nanmax(xSub)))
#                 newX = xinterp[maxXind]
#                 newZ = f(newX)
#                 xSubNew = np.append(xSub, newX)
#                 zSubNew = np.append(zSub, newZ)
#                 moredata = interpBathy(xSubNew, zSubNew, xinterp)
#
#                 del xSubNew, zSubNew,newZ,newX,f,mb,zbot,xbot
#
#                 nanValues2 = np.isnan(moredata['z'])
#
#                 if any(nanValues2):
#                     print('WHAT HAPPENED AT TRANSECT {}'.format(i))
#                     #plt.plot(moredata['x'], moredata['z'], 'r-')
#                 else:
#                     if count == 0:
#                         alllines = moredata['z']
#                         time = surveydate
#                         count = count + 1
#                     else:
#                         alllines = np.vstack((alllines, moredata['z']))
#                         time = np.append(time, surveydate)
#                         count = count + 1
#                     #plt.plot(moredata['x'], moredata['z'], 'k-')
#                     del moredata
#
#
#
#             elif len(xindex[0]) == 3:
#                 print('Could not do that: number of points being fit: {}'.format(len(xindex[0])))
#                 worstcount = worstcount+1
#             else:
#                 print('Could not do that: maximum X value in transect is: {}'.format(np.nanmax(xSub)))
#                 worstcount = worstcount+1
#                 #plt.plot(xSub, zSub, 'r-')
#
#
#         else:
#             if count == 0:
#                 alllines = newdata['z']
#                 time = surveydate
#                 count = count+1
#             else:
#                 alllines = np.vstack((alllines, newdata['z']))
#                 time = np.append(time,surveydate)
#                 count = count+1
#     else:
#         print('Survey with no data at this line {}'.format(i))
