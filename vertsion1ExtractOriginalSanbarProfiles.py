#   Lets identify a transect to the north to check on how much nourishment evolution was caught
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



#geomorphdir = '/media/dylananderson/Elements/filteredFRF_Geomorph/'
geomorphdir = '/media/dylananderson/Elements/FRF_Geomorph/'

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


subset = files.copy()
#from palettable.colorbrewer.sequential import Blues_8
#ax.imshow(data, cmap=Blues_8.mpl_colormap

colormap = plt.cm.gist_ncar

import scipy.stats as spstats

labels = []
# xinterp = np.arange(108, 608, 2.5)
#xinterp = np.arange(108, 618, 2.5)
xinterp = np.arange(0, 900, 2)
deeperThan = 200
surveyType = dict()
noneOnshore = np.empty((len(xinterp),))
noneOffshore = np.empty((len(xinterp),))
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
    file_params = subset[i].split('_')

    ### ### Southern Lines
    data1 = getBathy(os.path.join(geomorphdir, subset[i]), lower=-10, upper=20)
    # data2 = getBathy(os.path.join(geomorphdir, subset[i]), lower=20, upper=60)
    # data3 = getBathy(os.path.join(geomorphdir, subset[i]), lower=60, upper=100)
    # ### ### Northern Lines
    # data1 = getBathy(os.path.join(geomorphdir, subset[i]), lower=1070, upper=1100)
    # data2 = getBathy(os.path.join(geomorphdir, subset[i]), lower=1020, upper=1070)
    # data3 = getBathy(os.path.join(geomorphdir, subset[i]), lower=960, upper=1019)

    # ### ### MidNorthern Lines
    # data1 = getBathy(os.path.join(geomorphdir, subset[i]), lower=800, upper=840)
    # data2 = getBathy(os.path.join(geomorphdir, subset[i]), lower=740, upper=799)
    # data3 = getBathy(os.path.join(geomorphdir, subset[i]), lower=680, upper=739)
    ### ### MidSouthern Lines
    # data2 = getBathy(os.path.join(geomorphdir, subset[i]), lower=200, upper=250)
    # data3 = getBathy(os.path.join(geomorphdir, subset[i]), lower=251, upper=310)
    # data1 = getBathy(os.path.join(geomorphdir, subset[i]), lower=311, upper=360)
    # #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=990, upper=1100)

    #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-2, upper=10)

    temp = subset[i].split('_')

    surveydate = DT.datetime.strptime(temp[1], '%Y%m%d')
    elevs = data1['z']
    cross = data1['x']
    crossind = np.argsort(data1['x'])
    crossS = cross[crossind]
    elevsS = elevs[crossind]

    # elevs2 = data2['z']
    # cross2 = data2['x']
    # crossind2 = np.argsort(data2['x'])
    # crossS2 = cross2[crossind2]
    # elevsS2 = elevs2[crossind2]
    #
    # elevs3 = data3['z']
    # cross3 = data3['x']
    # crossind3 = np.argsort(data3['x'])
    # crossS3 = cross3[crossind3]
    # elevsS3 = elevs3[crossind3]

    xSub = np.ma.MaskedArray.filled(crossS, np.nan)
    zSub = np.ma.MaskedArray.filled(elevsS, np.nan)
    # xSub2 = np.ma.MaskedArray.filled(crossS2, np.nan)
    # zSub2 = np.ma.MaskedArray.filled(elevsS2, np.nan)
    # xSub3 = np.ma.MaskedArray.filled(crossS3, np.nan)
    # zSub3 = np.ma.MaskedArray.filled(elevsS3, np.nan)
    #realValues = ~np.isnan(xSub)

    # binnedz, bin_edges, binnumber = spstats.binned_statistic(xSub,zSub, statistic='mean',bins=np.arange(93,578,5))
    # bin_width = (bin_edges[1] - bin_edges[0])
    # bin_centers = bin_edges[1:] - bin_width / 2
    # xMSub = bin_centers
    # zMSub = binnedz
    # xSubNewB = xMSub[~np.isnan(zMSub)]
    # zSubNewB = zMSub[~np.isnan(zMSub)]
    # realValues = ~np.isnan(xSubNewB)
    # zSubNew = zSubNewB[~np.isnan(xSubNewB)]
    # xSubNew = xSubNewB[~np.isnan(xSubNewB)]

    realValues = ~np.isnan(xSub)
    xSubNew = xSub[~np.isnan(xSub)]
    zSubNew = zSub[~np.isnan(xSub)]

    # realValues2 = ~np.isnan(xSub2)
    # xSubNew2 = xSub2[~np.isnan(xSub2)]
    # zSubNew2 = zSub2[~np.isnan(xSub2)]
    #
    # realValues3 = ~np.isnan(xSub3)
    # xSubNew3 = xSub3[~np.isnan(xSub3)]
    # zSubNew3 = zSub3[~np.isnan(xSub3)]


    #newdata = interpBathy(xSubNew, zSubNew, xinterp)

#     fig = plt.figure(figsize=(10, 10))
#     plt.plot(cross, elevs, '.')
# #        plt.plot(xinterp,newdata['z'],'r.-')
#     plt.xlim([0, 600])
#     plt.ylim([-7, 2])
#     plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/' + temp[1][6:8])
#     plt.savefig('/home/dylananderson/projects/duckGeomorph/interpolatedProfiles/Survey_{}'.format(i))
#     plt.close()

    if any(realValues):

        if np.nanmax(xSubNew) > deeperThan:

            if np.nanmin(xSubNew) < 125:



                if len(xSubNew) > 6:


                    newdata = interpBathy(xSubNew, zSubNew, xinterp)
                    #newdata = interpBathy2(newdata['x'],newdata['z'],xinterp)
                    # if any(realValues2):
                    #     newdata2 = interpBathy(xSubNew2,zSubNew2, xinterp)
                    #     #newdata2 = interpBathy2(newdata2['x'], newdata2['z'], xinterp)
                    # else:
                    #     newdata2 = []
                    # if any(realValues3):
                    #     newdata3 = interpBathy(xSubNew3, zSubNew3, xinterp)
                    #     #newdata3 = interpBathy2(newdata3['x'], newdata3['z'], xinterp)
                    # else:
                    #     newdata3 = []
                    #newdata = interpBathy2(newdata['x'], newdata['z'], xinterp)
                    #newdata = interpBathy2(xSubNew, zSubNew, xinterp)

                    # fig = plt.figure(figsize=(10, 10))
                    # plt.plot(cross, elevs, 'b.')
                    # plt.plot(cross2, elevs2, 'b.')
                    # plt.plot(cross3, elevs3, 'b.')
                    # plt.plot(xinterp,newdata['z'],'r.-')
                    # if any(realValues2) & any(realValues3):
                    #
                    #     plt.plot(xinterp,newdata2['z'],'m.-')
                    #     plt.plot(xinterp, newdata3['z'], 'g.-')
                    #     plt.plot(xinterp,np.nanmean((newdata['z'], newdata2['z'],newdata3['z']), axis=0),'ko-')
                    #
                    # else:
                    #     if any(realValues2):
                    #         plt.plot(xinterp, newdata2['z'], 'c.-')
                    #         plt.plot(xinterp, np.nanmean((newdata['z'], newdata2['z']), axis=0), 'k-.')
                    #
                    #     elif any(realValues3):
                    #         plt.plot(xinterp, newdata3['z'], 'y.-')
                    #         plt.plot(xinterp, np.nanmean((newdata['z'], newdata3['z']), axis=0), 'k--')
                    #
                    #     else:
                    #         plt.plot(xinterp, newdata['z'], 'k-')

                    # plt.xlim([0, 600])
                    # plt.ylim([-7, 2])
                    # plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/' + temp[1][6:8] + ' Survey Type = ' + file_params[6])
                    # plt.savefig('/home/dylananderson/projects/duckGeomorph/interpolatedProfiles/Survey_{}'.format(i))
                    # plt.close()

                    # nanValues = np.isnan(newdata['z'])
                    #
                    # if any(nanValues):
                    #     #print('Found a transect with nans {}'.format(i))
                    #     print('Trying to extend the lower half of the profile with a linear fit: {}'.format(surveydate))
                    #     if count2 == 0:
                    #         badlines = newdata['z']
                    #         count2 = count2+1
                    #         badtime = surveydate
                    #     else:
                    #         badlines = np.vstack((badlines, newdata['z']))
                    #         count2 = count2+1
                    #         badtime = np.append(badtime, surveydate)

                        #print('Threw out a survey with data due to nans: {}'.format(surveydate))
                        #plt.plot(xSub,zSub)


                        # xindex = np.where((xSubNew>=deeperThan))
                        #
                        # if len(xindex[0]) > 3:
                        #
                        #     xbot = xSubNew[xindex]
                        #     zbot = zSubNew[xindex]
                        #     mb = np.polyfit(xbot, zbot, 1)
                        #     f = np.poly1d(mb)
                        #     maxXind = np.where((xinterp > np.nanmax(xSubNew)))
                        #     newX = xinterp[maxXind]
                        #     newZ = f(newX)
                        #     xSubNew2 = np.append(xSubNew, newX)
                        #     zSubNew2 = np.append(zSubNew, newZ)
                        #     moredata = interpBathy(xSubNew2, zSubNew2, xinterp)
                        #     #moredata = interpBathy(moredata['x'], moredata['z'], xinterp)
                        #     #moredata = interpBathy(xSubNew2, moredata['z'], xinterp)
                        #
                        #     del xSubNew2, zSubNew2,newZ,newX,f,mb,zbot,xbot
                        #
                        #     nanValues2 = np.isnan(moredata['z'])
                        #
                        #     if any(nanValues2):
                        #         print('WHAT HAPPENED AT TRANSECT {}'.format(i))
                        #         #plt.plot(moredata['x'], moredata['z'], 'r-')
                        #     else:
                        #         if count == 0:
                        #             alllines = moredata['z']
                        #             time = surveydate
                        #             count = count + 1
                        #         else:
                        #             alllines = np.vstack((alllines, moredata['z']))
                        #             time = np.append(time, surveydate)
                        #             count = count + 1
                        #
                        #
                        #         if count3 == 0:
                        #             extrapolatedlines = newdata['z']
                        #             count3 = count3 + 1
                        #             extrapolatedtime = surveydate
                        #         else:
                        #             extrapolatedlines = np.vstack((extrapolatedlines, newdata['z']))
                        #             count3 = count3 + 1
                        #             extrapolatedtime = np.append(extrapolatedtime, surveydate)
                        #         #plt.plot(moredata['x'], moredata['z'], 'k-')
                        #         del moredata
                        #
                        #
                        #
                        # elif len(xindex[0]) == 3:
                        #     print('Could not do that: number of points being fit: {}'.format(len(xindex[0])))
                        #     worstcount = worstcount+1
                        # else:
                        #     print('Could not do that: maximum X value in transect is: {}'.format(np.nanmax(xSub)))
                        #     worstcount = worstcount+1
                        #     #plt.plot(xSub, zSub, 'r-')

                    if count == 0:


                        alllines = newdata['z']
                        surveyType = file_params[6]
                        time = surveydate
                        count = count+1

                    else:
                        alllines = np.vstack((alllines, newdata['z']))
                        surveyType = np.append(surveyType,file_params[6])
                        time = np.append(time,surveydate)
                        count = count+1
                    # else:
                    #     if count == 0:
                    #         if any(realValues2) & any(realValues3):
                    #             avgProf = np.nanmean((newdata['z'],newdata2['z'],newdata3['z']),axis=0)
                    #         else:
                    #             if any(realValues2):
                    #                 avgProf = np.nanmean((newdata['z'], newdata2['z']), axis=0)
                    #             elif any(realValues3):
                    #                 avgProf = np.nanmean((newdata['z'], newdata3['z']), axis=0)
                    #             else:
                    #                 avgProf = newdata['z']
                    #
                    #         alllines = avgProf
                    #         surveyType = file_params[6]
                    #         time = surveydate
                    #         count = count+1
                    #     else:
                    #         if any(realValues2) & any(realValues3):
                    #             avgProf = np.nanmean((newdata['z'],newdata2['z'],newdata3['z']),axis=0)
                    #         else:
                    #             if any(realValues2):
                    #                 avgProf = np.nanmean((newdata['z'], newdata2['z']), axis=0)
                    #             elif any(realValues3):
                    #                 avgProf = np.nanmean((newdata['z'], newdata3['z']), axis=0)
                    #             else:
                    #                 avgProf = newdata['z']
                    #         alllines = np.vstack((alllines, avgProf))
                    #         surveyType = np.append(surveyType,file_params[6])
                    #         time = np.append(time,surveydate)
                    #         count = count+1

                        # fig = plt.figure(figsize=(10,10))
                        # plt.plot(cross,elevs,'.')
                        # plt.plot(xinterp,newdata['z'],'.-')
                        # plt.xlim([0, 600])
                        # plt.ylim([-7,2])
                        # plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
                        # plt.savefig('/home/dylananderson/projects/duckGeomorph/interpolatedProfiles/Survey_{}'.format(i))
                        # plt.close()
                else:
                    print('Data is less than 10 points long at line {}'.format(i))
            else:
                print('No data onshore of 125 meters at line {}'.format(i))
                print('Most onshore point at {}'.format(np.nanmin(xSubNew)))
                newdata = interpBathy(xSubNew, zSubNew, xinterp)

                # if count4 == 0:
                #     noneOnshore = newdata['z']
                #     noOnshoreTime = surveydate
                #     count4 = count4 + 1
                # else:
                #     noneOnshore = np.vstack((noneOnshore, newdata['z']))
                #     noOnshoreTime = np.append(noOnshoreTime, surveydate)
                #     count4 = count4 + 1
        else:
            print('No data deeper than 200 meters at line {}'.format(i))
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

        # plt.plot(xSubNew,zSubNew)
plt.style.use('dark_background')


#
#
# bathy = dict()
# #alllines = np.empty((len(xinterp),))
# count = 0
# count2 = 0
# worstcount = 0
# #fig = plt.figure(figsize=(10,10))
# for i in range(len(subset)):
#     #data = getBathy(os.path.join(geomorphdir, subset[i]), lower=1080, upper=1100)
#     data = getBathy(os.path.join(geomorphdir, subset[i]), lower=-10, upper=20)
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
#             xindex = np.where((xSub>=deeperThan))
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
volumes = np.empty(len(alllines,))
from numpy import trapz
for i in range(len(alllines)):
    volumes[i] = trapz(alllines[i,:]+10,dx=5)


# # Plotting all of the lines
t1 = 0
t2 = -1
subsetOfAllines = alllines[t1:t2]
fig, ax = plt.subplots(2,1)
ax[0].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subsetOfAllines)))))
for i in range(len(subsetOfAllines)):
    ax[0].plot(xinterp, subsetOfAllines[i,:], label=time[i])

#ax[0].legend(loc='upper right')
ax[0].set_xlim([50, 850])
ax[0].set_ylim([-8, 4])
ax[0].set_title('Cross-shore profile variability')

plt.set_cmap('RdBu_r')

tg, xg = np.meshgrid(time, xinterp)
plt0 = ax[1].pcolor(xg,tg,(alllines-np.nanmean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
fig.colorbar(plt0, ax=ax[1], orientation='horizontal')
ax[1].set_ylim([time[t1], time[t2]])
ax[1].set_title('Surveys (dev.)')

plt.figure()
#plt.plot(time,volumes)
groups = np.unique(surveyType)
for name in groups:
    matched_indexes = []
    i = 0
    length = len(surveyType)

    while i < length:
        if name == surveyType[i]:
            matched_indexes.append(i)
        i += 1

    plt.scatter(time[matched_indexes],volumes[matched_indexes],marker='o',label=name)
plt.legend()
plt.ylabel('Volumes (m^3/m alongshore)')

tg, xg = np.meshgrid(time, xinterp)

morphoPickle = 'sandbarsOriginalSouthTransect.pickle'
output = {}
output['time'] = time
output['alllines'] = alllines
output['xinterp'] = xinterp
# output['avgFY'] = -fygrid
# output['Mgrid'] = Mgrid
# output['lw'] = lw
# output['norm'] = norm
import pickle
with open(morphoPickle,'wb') as f:
    pickle.dump(output, f)


plt.style.use('default')
plt.figure()
ax = plt.subplot2grid((3,3,),(0,0),rowspan=2,colspan=3)
ax.plot(xinterp, alllines[0, :], color=[0.7, 0.7, 0.7, 1], linewidth=0.5,label='surveys')
for xx in range(len(alllines)):
    ax.plot(xinterp, alllines[xx,:], color=[0.7,0.7,0.7,1],linewidth=0.5)
ax.plot(xinterp, np.nanmean(alllines,axis=0),'k-',label='mean')
ax.plot(xinterp, np.nanmean(alllines,axis=0) + np.nanstd(alllines,axis=0),'k--',label='std')
ax.plot(xinterp, np.nanmean(alllines,axis=0) - np.nanstd(alllines,axis=0),'k--')
ax.legend()
ax2 = plt.subplot2grid((3,3,),(2,0),rowspan=1,colspan=3)
ax2.plot(xinterp, np.nanstd(alllines,axis=0),'k--')
ax2.set_ylabel('std (m)')
ax2.set_xlabel('xFRF coords. (m)')
ax.set_ylabel('Elevation relative to NAVD88 (m)')

# import pickle
# #dbfile = open('oct22wam1571770800FlowFields_unfiltered.pickle', 'rb')
# dbfile = open('oct22wam1571770800FlowFields_filtered.pickle', 'rb')
#
# data = pickle.load(dbfile)
# dbfile.close()