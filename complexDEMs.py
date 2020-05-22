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

def P2R(radii, angles):
    return radii * np.exp(1j*angles)

def R2P(x):
    return np.abs(x), np.angle(x)


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
    yind = np.where((data['y'] > 0) & (data['y'] < 1000))
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
plt.pcolor(portionx,portiony,allBathyMean,vmin=-7,vmax=1)

plt.figure()
plt.plot(portionx,allBathyMean[15,:])



demean = np.zeros((np.shape(allBathy)[1]*np.shape(allBathy)[2],np.shape(allBathy)[0]))
for i in range(len(allBathy)):
    data = allBathy[i]-allBathyMean
    demean[:,i] = data.flatten()
    t = time[i]

# demean needs to have dimensions (time,num of points in DEM)
demean = demean

#demean = alllines - np.mean(alllines,axis=0)

from scipy.signal import hilbert
#from scipy.fftpack import hilbert
data = (hilbert(demean.T))

# data = (hilbert(demean))

#datareal = np.imag(hilbert(alllines))
#dataimag = np.real(hilbert(alllines))

# data should have dims (time, num of points)
# data = data.T
# c should have dims (time,time)
c = np.matmul(np.conj(data).T,data)/np.shape(data)[0]


import scipy.linalg as la
import numpy.linalg as npla

lamda, loadings = la.eigh(c)

lamda2, loadings2 = npla.eig(c)

ind = np.argsort(lamda[::-1])

lamda[::-1].sort()

loadings = loadings[:,ind]

# pcs should have dims (time,num of points)
pcs = np.dot(data, loadings)# / np.sqrt(lamda)
loadings = loadings# * np.sqrt(lamda)
# pcsreal = np.real(pcs)
# pcsimag = np.imag(pcs)
# eofreal = np.real(loadings)
# eofimag = np.imag(loadings)
pcsreal = np.real(pcs[:,0:200])
pcsimag = np.imag(pcs[:,0:200])
eofreal = np.real(loadings[:,0:200])
eofimag = np.imag(loadings[:,0:200])

S = np.power(loadings*np.conj(loadings),0.5) * np.sqrt(lamda)

theta = np.arctan2(eofimag,eofreal)
theta2 = theta*180/np.pi

Rt = np.power(pcs*np.conj(pcs),0.5) / np.sqrt(lamda)

phit = np.arctan2(pcsimag,pcsreal)
phit2 = phit*180/np.pi

# mode = 0

# fig, ax = plt.subplots(2,2)
#
# ax[0,0].plot(xinterp, S[:,mode],'o')
# ax[0,1].plot(xinterp, theta2[:,mode],'o')
# ax[1,0].plot(time,Rt[:,mode],'o')
# ax[1,1].plot(time,phit2[:,mode],'o')
totalV = np.sum(lamda)
percentV = lamda / totalV

mode = 1

fig, ax = plt.subplots(2,2)
axS = ax[0,0].pcolor(portionx,portiony,np.real(np.reshape(S[:,mode],(np.shape(allBathy)[1],np.shape(allBathy)[2]))))
fig.colorbar(axS,ax=ax[0,0])
axT = ax[0,1].pcolor(portionx,portiony,np.real(np.reshape(theta2[:,mode],(np.shape(allBathy)[1],np.shape(allBathy)[2]))))
fig.colorbar(axT,ax=ax[0,1])

ax[1,0].plot(time,Rt[:,mode],'o')
ax[1,1].plot(time,phit2[:,mode],'o')





mode = 1

mag = np.real(S[:,mode])
phase = np.real(theta2[:,mode])*np.pi/180

def P2R(radii, angles):
    return radii * np.exp(1j*angles)

u = np.zeros((np.shape(mag)))
v = np.zeros((np.shape(mag)))
for i in range(len(mag)):
    utemp = P2R(np.real(mag[i]), np.real(phase[i]))
    u[i] = np.real(utemp)
    v[i] = np.imag(utemp)

U = np.real(np.reshape(u,(np.shape(allBathy)[1],np.shape(allBathy)[2])))
V = np.real(np.reshape(v,(np.shape(allBathy)[1],np.shape(allBathy)[2])))


fig = plt.figure(figsize=(14,6))

ax1 = plt.subplot2grid((2,4),(0,0),rowspan=2,colspan=2)
ax1.quiver(portionx, portiony, U, V, mag)
ax2 = plt.subplot2grid((2,4),(0,2),rowspan=1,colspan=2)
ax2.plot(time,Rt[:,mode],'o')
ax3 = plt.subplot2grid((2,4),(1,2),rowspan=1,colspan=2)
ax3.plot(time,phit2[:,mode],'o')
ax1.set_title('Mode 4: {}%'.format(np.round(100*percentV[3])))
ax1.set_xlabel('xFRF (m)')
ax1.set_ylabel('yFRF (m)')
ax2.set_ylabel('magnitude')
ax3.set_ylabel('phase')
ax3.set_xlabel('time')
eofPred = np.zeros(np.shape(S[:,0]))
# for pp in range(len(time)):
#
#     mode = 1
#
#     # eofPred[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
#     #         timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
#
#     mag = np.real(S[:, mode])
#     phase = np.real(theta[:, mode]) # * np.pi / 180
#     mult = np.real(Rt[pp,mode])
#     offset = phit[pp,mode]
#
#     eofPred[:] = Rt[pp, mode] * np.sin(phit[pp, mode]) * S[:, mode] * np.sin(theta[:, mode]) + Rt[
#         pp, mode] * np.cos(phit[pp, mode]) * S[:, mode] * np.cos(theta[:, mode])
#
#     fig = plt.figure(figsize=(14, 6))
#     ax1 = plt.subplot2grid((2, 4), (0, 0), rowspan=2, colspan=2)
#     axPC = ax1.pcolor(portionx, portiony, np.real(np.reshape(eofPred, (np.shape(allBathy)[1], np.shape(allBathy)[2]))),vmin=-1,vmax=1,cmap='RdBu')
#     fig.colorbar(axPC,ax=ax1)
#     #
#     # u = np.zeros((np.shape(mag)))
#     # v = np.zeros((np.shape(mag)))
#     # for i in range(len(mag)):
#     #     utemp = P2R(np.multiply(np.real(mag[i]),mult), np.add(np.real(phase[i]),offset))
#     #     u[i] = np.real(utemp)
#     #     v[i] = np.imag(utemp)
#     #
#     # U = np.real(np.reshape(u, (np.shape(allBathy)[1], np.shape(allBathy)[2])))
#     # V = np.real(np.reshape(v, (np.shape(allBathy)[1], np.shape(allBathy)[2])))
#     #
#     # fig = plt.figure(figsize=(14, 6))
#     #
#     # ax1 = plt.subplot2grid((2, 4), (0, 0), rowspan=2, colspan=2)
#     # ax1.quiver(portionx, portiony, U, V, mag)
#     ax2 = plt.subplot2grid((2, 4), (0, 2), rowspan=1, colspan=2)
#     ax2.plot(time, Rt[:, mode], 'o')
#     ax2.plot(time[pp],Rt[pp,mode],'mo')
#     ax3 = plt.subplot2grid((2, 4), (1, 2), rowspan=1, colspan=2)
#     ax3.plot(time, phit2[:, mode], 'o')
#     ax3.plot(time[pp], phit2[pp,mode],'mo')
#     ax1.set_title('Mode 2: {}%'.format(time[pp]))
#     ax1.set_xlabel('xFRF (m)')
#     ax1.set_ylabel('yFRF (m)')
#     ax2.set_ylabel('magnitude')
#     ax3.set_ylabel('phase')
#     ax3.set_xlabel('time')
#     if pp < 10:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/ceofRealMovie2/Survey_00{}'.format(pp))
#     elif pp < 100:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/ceofRealMovie2/Survey_0{}'.format(pp))
#     else:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/ceofRealMovie2/Survey_{}'.format(pp))
#     plt.close()



cpc1 = P2R(Rt[:, 0], phit[:, 0])
cpc2 = P2R(Rt[:, 1], phit[:, 1])
cpc3 = P2R(Rt[:, 2], phit[:, 2])
cpc4 = P2R(Rt[:, 3], phit[:, 3])
cpc5 = P2R(Rt[:, 4], phit[:, 4])
cpc6 = P2R(Rt[:, 5], phit[:, 5])
cpc7 = P2R(Rt[:, 6], phit[:, 6])
cpc8 = P2R(Rt[:, 7], phit[:, 7])

modes = 8

CPCs = np.vstack((np.real(cpc1), np.imag(cpc1), np.real(cpc2), np.imag(cpc2), np.real(cpc3), np.imag(cpc3),
                  np.real(cpc4), np.imag(cpc4), np.real(cpc5), np.imag(cpc5), np.real(cpc6), np.imag(cpc6),
                  np.real(cpc7), np.imag(cpc7), np.real(cpc8), np.imag(cpc8),))
CPCs = CPCs.T

var_explained = np.array((percentV[0], percentV[0], percentV[1], percentV[1], percentV[2], percentV[2],
                          percentV[3], percentV[3], percentV[4], percentV[4], percentV[5], percentV[5],
                          percentV[6], percentV[6], percentV[7], percentV[7]))


fig = plt.figure()
plt.plot(np.real(cpc1),np.imag(cpc1))
plt.plot(np.real(cpc2),np.imag(cpc2))
plt.plot(np.real(cpc3),np.imag(cpc3))
plt.plot(np.real(cpc4),np.imag(cpc4))
plt.plot(np.real(cpc5),np.imag(cpc5))
plt.plot(np.real(cpc6),np.imag(cpc6))
plt.plot(np.real(cpc7),np.imag(cpc7))

fig = plt.figure()
plt.plot(np.real(cpc1)*var_explained[0],np.imag(cpc1)*var_explained[0])
plt.plot(np.real(cpc2)*var_explained[2],np.imag(cpc2)*var_explained[2])
plt.plot(np.real(cpc3)*var_explained[4],np.imag(cpc3)*var_explained[4])
plt.plot(np.real(cpc4)*var_explained[6],np.imag(cpc4)*var_explained[6])
plt.plot(np.real(cpc5)*var_explained[8],np.imag(cpc5)*var_explained[8])

temp = CPCs*var_explained
fig = plt.figure()
plt.plot(temp[:,0], temp[:,1])
plt.plot(temp[:,2], temp[:,3])
plt.plot(temp[:,4], temp[:,5])
plt.plot(temp[:,6], temp[:,7])
plt.plot(temp[:,8], temp[:,9])


from sklearn.cluster import KMeans
numClusters = 10

PCsub = CPCs*var_explained

#  KMEANS
kma = KMeans(n_clusters=numClusters, n_init=2000).fit(PCsub)

# groupsize
_, group_size = np.unique(kma.labels_, return_counts=True)

# groups
d_groups = {}
for k in range(numClusters):
    d_groups['{0}'.format(k)] = np.where(kma.labels_ == k)

# # centroids
# centroids = np.dot(kma.cluster_centers_, EOFsub)
# centroids = kma.cluster_centers_
centroids = np.zeros((numClusters, PCsub.shape[1]))
for k in range(numClusters):
    print(PCsub[d_groups['{0}'.format(k)],0:1])
    centroids[k,:] = np.mean(PCsub[d_groups['{0}'.format(k)],:], axis=1)/var_explained

bmus = kma.labels_



from matplotlib import gridspec
fig2 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(3, 5, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
dems = np.zeros((numClusters,len(S[:,0])))
magPhaseCentroid = np.zeros((np.shape(centroids)))
colormap = 'RdBu'
for i in range(numClusters):
    #getind = sorted[i]
    # cpc1RT, cpc1Phi = R2P(complex(centroids[i, 0], centroids[i, 1]))
    # cpc2RT, cpc2Phi = R2P(complex(centroids[i, 2], centroids[i, 3]))
    # cpc3RT, cpc3Phi = R2P(complex(centroids[i, 4], centroids[i, 5]))
    # cpc4RT, cpc4Phi = R2P(complex(centroids[i, 6], centroids[i, 7]))
    # cpc5RT, cpc5Phi = R2P(complex(centroids[i, 8], centroids[i, 9]))
    # cpc6RT, cpc6Phi = R2P(complex(centroids[i, 10], centroids[i, 11]))

    cpcRt = np.zeros((modes,))
    cpcPhi = np.zeros((modes,))

    cpcRt[0], cpcPhi[0] = R2P(complex(centroids[i, 0], centroids[i, 1]))
    cpcRt[1], cpcPhi[1] = R2P(complex(centroids[i, 2], centroids[i, 3]))
    cpcRt[2], cpcPhi[2] = R2P(complex(centroids[i, 4], centroids[i, 5]))
    cpcRt[3], cpcPhi[3] = R2P(complex(centroids[i, 6], centroids[i, 7]))
    cpcRt[4], cpcPhi[4] = R2P(complex(centroids[i, 8], centroids[i, 9]))
    cpcRt[5], cpcPhi[5] = R2P(complex(centroids[i, 10], centroids[i, 11]))

    cpcRt[6], cpcPhi[6] = R2P(complex(centroids[i, 8], centroids[i, 9]))
    cpcRt[7], cpcPhi[7] = R2P(complex(centroids[i, 10], centroids[i, 11]))

    magPhaseCentroid[i,0:16] = [cpcRt[0], cpcPhi[0], cpcRt[1], cpcPhi[1], cpcRt[2], cpcPhi[2],
                                cpcRt[3], cpcPhi[3], cpcRt[4], cpcPhi[4], cpcRt[5], cpcPhi[5],
                                cpcRt[6], cpcPhi[6], cpcRt[7], cpcPhi[7]]


    profile = 0 * np.ones(len(S[:,0]), )
    for mode in range(6):
        profile = profile + cpcRt[mode] * np.sin(cpcPhi[mode]) * S[:, mode] * np.sin(theta[:, mode]) + cpcRt[mode]* np.cos(cpcPhi[mode]) * S[:, mode] * np.cos(theta[:, mode])

    dems[i,:] = profile
    #profile = kma.centroids[i,:] * pred_std + pred_mean

    ax = plt.subplot(gs[gr, gc])

    # ax.plot(xinterp,profile+np.mean(alllines,axis=0))
    ax.pcolor(portionx,portiony, np.real(np.reshape(profile, (np.shape(allBathy)[1], np.shape(allBathy)[2]))),vmin=-1.25,vmax=1.25,cmap=colormap)

    # ax.set_xlim([80, 820])
    # ax.set_ylim([-8, 2.25])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    ax.text(480,900, group_size[i], color='k',fontweight='bold')

    if gc > 0:
        ax.set_yticks([])

    if gr < (5-1):
        ax.set_xticks([])
    #  counter
    gc += 1
    if gc >= 5:
        gc = 0
        gr += 1

import matplotlib.cm as cm
import matplotlib.colors as mcolors
normalize = mcolors.Normalize(vmin=-1.25, vmax=1.25)
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
# s_map.set_array(colorparam)
fig2.subplots_adjust(right=0.92)
cbar_ax = fig2.add_axes([0.94, 0.10, 0.02, 0.7])
cbar = fig2.colorbar(s_map, cax=cbar_ax)
cbar.set_label('Deviation from Mean (m)')


magEOF1cen = magPhaseCentroid[:,0]
phaseEOF1cen = magPhaseCentroid[:,1]*180/np.pi
sortedPeaks = np.sort(phaseEOF1cen)
sortedPeakInd = np.argsort(phaseEOF1cen)
sortedPhases = phaseEOF1cen[sortedPeakInd]



from matplotlib import gridspec
fig2 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(numClusters, 1, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
colormap = 'RdBu_r'
for i in range(numClusters):

    profile = dems[sortedPeakInd[i],:]
    ax = plt.subplot(gs[gr, gc])

    # ax.plot(xinterp,profile+np.mean(alllines,axis=0))
    ax.pcolor(portionx,portiony, np.real(np.reshape(profile, (np.shape(allBathy)[1], np.shape(allBathy)[2]))),vmin=-1.25,vmax=1.25,cmap=colormap)

    # ax.set_xlim([80, 820])
    # ax.set_ylim([-8, 2.25])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    ax.text(480,900, group_size[sortedPeakInd[i]], color='k',fontweight='bold')
    # ax.text(480,900, magEOF1cen[i], color='k',fontweight='bold')
    if i == 0:
        plt.title('Clusters of DEMs')
    if gc > 0:
        ax.set_yticks([])

    if gr < (numClusters-1):
        ax.set_xticks([])
    #  counter
    gr += 1
    if gr >= numClusters:
        gc += 1
        gr = 0


import matplotlib.cm as cm
import matplotlib.colors as mcolors
normalize = mcolors.Normalize(vmin=-1.25, vmax=1.25)
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
# s_map.set_array(colorparam)
fig2.subplots_adjust(right=0.86)
cbar_ax = fig2.add_axes([0.9, 0.10, 0.02, 0.7])
cbar = fig2.colorbar(s_map, cax=cbar_ax)
cbar.set_label('Deviation from Mean (m)')




from matplotlib import gridspec
fig2 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(5, 2, wspace=0.0, hspace=0.0)
gr, gc = 0, 0
colormap = 'RdBu_r'
for i in range(numClusters):

    profile = dems[sortedPeakInd[i],:]
    ax = plt.subplot(gs[gr, gc])

    # ax.plot(xinterp,profile+np.mean(alllines,axis=0))
    ax.pcolor(portionx,portiony, np.real(np.reshape(profile, (np.shape(allBathy)[1], np.shape(allBathy)[2]))),vmin=-1.25,vmax=1.25,cmap=colormap)

    # ax.set_xlim([80, 820])
    # ax.set_ylim([-8, 2.25])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    ax.text(480,900, group_size[sortedPeakInd[i]], color='k',fontweight='bold')
    # ax.text(480,900, magEOF1cen[i], color='k',fontweight='bold')
    if i == 0:
        plt.title('Clusters of DEMs')
    if gc > 0:
        ax.set_yticks([])

    if gr < 4:
        ax.set_xticks([])
    #  counter
    gr += 1
    if gr >= 5:
        gc += 1
        gr = 0


import matplotlib.cm as cm
import matplotlib.colors as mcolors
normalize = mcolors.Normalize(vmin=-1.25, vmax=1.25)
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
# s_map.set_array(colorparam)
fig2.subplots_adjust(right=0.86)
cbar_ax = fig2.add_axes([0.9, 0.10, 0.02, 0.7])
cbar = fig2.colorbar(s_map, cax=cbar_ax)
cbar.set_label('Deviation from Mean (m)')



#
#     plt.figure(figsize=(10,10))
#     plt.pcolor(portionx,portiony,data,vmin=-2,vmax=2,cmap='RdBu_r')
#     plt.title('Survey Date = ' + temp[1][0:4] + '/' + temp[1][4:6] + '/'+ temp[1][6:8])
#     plt.title('{}'.format(t))
#     plt.colorbar()
#     if i < 10:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDiff/Survey_00{}'.format(i))
#     elif i < 100:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDiff/Survey_0{}'.format(i))
#     else:
#         plt.savefig('/home/dylananderson/projects/duckGeomorph/bathyDiff/Survey_{}'.format(i))
#     plt.close()

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
# geomorphdir = '/home/dylananderson/projects/duckGeomorph/bathyDiff/'
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
# video = cv2.VideoWriter('croppedBathyDiff.avi', forcc, 8, (width, height))
# for image in files_path:
#     video.write(cv2.imread(image))
# cv2.destroyAllWindows()
# video.release()

