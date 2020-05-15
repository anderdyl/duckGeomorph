#   Lets load a couple DEMs and get a work flow running for 2D complex EOFs
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



geomorphdir = '/media/dylananderson/Elements/FRF_DEMs/'
files = os.listdir(geomorphdir)
files.sort()
files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]

def getBathy(file):
    bathy = Dataset(file)

    xs_bathy = bathy.variables['xFRF'][:]
    ys_bathy = bathy.variables['yFRF'][:]
    zs_bathy = bathy.variables['elevation'][:]
    ts_bathy = bathy.variables['time'][:]
    #pr_bathy = bathy.variables['profileNumber'][:]

    # zs_bathy = np.ma.masked_where((pr_bathy > upper), zs_bathy)
    # ys_bathy = np.ma.masked_where((pr_bathy > upper), ys_bathy)
    # xs_bathy = np.ma.masked_where((pr_bathy > upper), xs_bathy)
    # pr_bathy = np.ma.masked_where((pr_bathy > upper), pr_bathy)
    # ts_bathy = np.ma.masked_where((pr_bathy > upper), ts_bathy)
    #
    # zs_bathy = np.ma.masked_where((pr_bathy < lower), zs_bathy)
    # ys_bathy = np.ma.masked_where((pr_bathy < lower), ys_bathy)
    # xs_bathy = np.ma.masked_where((pr_bathy < lower), xs_bathy)
    # pr_bathy = np.ma.masked_where((pr_bathy < lower), pr_bathy)
    # ts_bathy = np.ma.masked_where((pr_bathy < lower), ts_bathy)

    output = dict()
    output['x'] = xs_bathy
    output['y'] = ys_bathy
    output['z'] = zs_bathy
    # output['pr'] = pr_bathy
    output['t'] = ts_bathy

    return output



subset = files[0:425]
plt.style.use('dark_background')


bathy = dict()
#alllines = np.empty((len(xinterp),))
count = 0
count2 = 0
worstcount = 0
#fig = plt.figure(figsize=(10,10))
allBathy = []
time = []
for i in range(len(subset)):
    data = getBathy(os.path.join(geomorphdir, subset[i]))

    temp = subset[i].split('_')

    surveydate = DT.datetime.strptime(temp[1], '%Y%m%d')
    xind = np.where((data['x'] > 100) & (data['x'] < 550))
    yind = np.where((data['y'] > 600) & (data['y'] < 1000))
    portion = data['z'][0,yind[0][0]:yind[0][-1]+1,xind[0][0]:xind[0][-1]+1]
    nanValues = np.isnan(portion)

    if not np.ma.array(portion).mask.any():
        allBathy.append(np.ma.filled(portion))
        time.append(surveydate)

        # plt.figure(figsize=(10,10))
        # plt.pcolor(data['x'][xind[0]],data['y'][yind[0]],portion,vmin=-7,vmax=2)
        # plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
        # if i < 10:
        #     plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDEMextents/Survey_00{}'.format(i))
        # elif i < 100:
        #     plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDEMextents/Survey_0{}'.format(i))
        # else:
        #     plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDEMextents/Survey_{}'.format(i))
        # plt.close()
#
#     plt.figure(figsize=(10,10))
#     plt.pcolor(data['x'],data['y'],data['z'][0,:,:],vmin=-7,vmax=2)
#     plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
#     if i < 10:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDEMextents/Survey_00{}'.format(i))
#     elif i < 100:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDEMextents/Survey_0{}'.format(i))
#     else:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDEMextents/Survey_{}'.format(i))
#     plt.close()
#
#
# geomorphdir = '/home/dylananderson/projects/duckGeomorph/bathyDEMextents/'
#
# files = os.listdir(geomorphdir)
#
# files.sort()
#
# files_path = [os.path.join(geomorphdir,x) for x in os.listdir(geomorphdir)]
#
# files_path.sort()
# import cv2
#
# frame = cv2.imread(files_path[0])
# height, width, layers = frame.shape
# forcc = cv2.VideoWriter_fourcc(*'XVID')
# video = cv2.VideoWriter('croppedNoNANsDEMsIntBathy.avi', forcc, 8, (width, height))
# for image in files_path:
#     video.write(cv2.imread(image))
# cv2.destroyAllWindows()
# video.release()



portionx = data['x'][xind[0]]
portiony = data['y'][yind[0]]

allBathyMean = np.mean(allBathy, axis=0)
plt.figure(figsize=(10,10))
plt.pcolor(portionx,portiony,allBathyMean,vmin=-7,vmax=2)

plt.figure()
plt.plot(portionx,allBathyMean[10,:])




for i in range(len(allBathy)):
    data = allBathy[i]-allBathyMean
    t = time[i]

    plt.figure(figsize=(10,10))
    plt.pcolor(portionx,portiony,data,vmin=-2,vmax=2,cmap='RdBu_r')
    plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
    plt.title('{}'.format(t))
    plt.colorbar()
    if i < 10:
        plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDiff/Survey_00{}'.format(i))
    elif i < 100:
        plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDiff/Survey_0{}'.format(i))
    else:
        plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDiff/Survey_{}'.format(i))
    plt.close()

#
#     plt.figure(figsize=(10,10))
#     plt.pcolor(data['x'],data['y'],data['z'][0,:,:],vmin=-7,vmax=2)
#     plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
#     if i < 10:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDEMextents/Survey_00{}'.format(i))
#     elif i < 100:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDEMextents/Survey_0{}'.format(i))
#     else:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDEMextents/Survey_{}'.format(i))
#     plt.close()
#
#
geomorphdir = '/home/dylananderson/projects/duckGeomorph/bathyDiff/'

files = os.listdir(geomorphdir)

files.sort()

files_path = [os.path.join(geomorphdir,x) for x in os.listdir(geomorphdir)]

files_path.sort()
import cv2

frame = cv2.imread(files_path[0])
height, width, layers = frame.shape
forcc = cv2.VideoWriter_fourcc(*'XVID')
video = cv2.VideoWriter('croppedBathyDiff.avi', forcc, 8, (width, height))
for image in files_path:
    video.write(cv2.imread(image))
cv2.destroyAllWindows()
video.release()

