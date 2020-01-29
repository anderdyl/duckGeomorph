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
import xarray as xr
#from sandBarTool import morphLib
#from downloads import pyeemd

#file = 'merged1572548400.nc'
#file = 'merged1572550200.nc'
#file = 'test1572550200.nc'
file = '/home/dylananderson/projects/drifters/testMerged_highres2.nc'
#data = xr.open_dataset('merged1572548400.nc')
data = xr.open_dataset(file)

#data = xr.open_dataset('merged1572539400.nc')
merged = data['merged'].values
x = data['xFRF'].values
y = data['yFRF'].values
imTime = data['t'].values

stack = np.fliplr(merged[200:250,400,0:240].T)

ystack = y[200:250]
xstack = x[400]*np.ones((np.shape(ystack)))
#time = np.arange(1960, 2020, 1/12)

#x = np.arange(100,600,10)

# z = np.nan * np.zeros((len(time),len(x)))
# for t in range(len(time)):
#     temp = time[t]-np.floor(time[t])
#     if temp > 0.66:
#         z[t,:] = 2.*np.cos(-2*np.pi*x/250 + 2*np.pi*time[t]/10) + 0.5*np.cos(2*np.pi*x/200 + 2*np.pi*time[t]/5)
#
#     elif temp > 0.33:
#         z[t,:] = 2.*np.cos(-2*np.pi*x/250 + 2*np.pi*time[t]/10) + 0.5*np.cos(2*np.pi*x/200 + 2*np.pi*time[t]/5)
#
#     else:
#         z[t,:] = 2.*np.cos(-2*np.pi*x/250 + 2*np.pi*time[t]/10) + 0.5*np.cos(2*np.pi*x/200 + 2*np.pi*time[t]/5)
#

# z = np.nan * np.zeros((len(time),len(x)))
# for t in range(len(time)):
#     if time[t] > 2015:
#         z[t,:] = 2.*np.cos(-2*np.pi*x/250 + 2*np.pi*time[t]/10)
#     elif time[t] > 2010:
#         z[t,:] = 2.*np.cos(2*np.pi*x/250 + 2*np.pi*time[t]/10)
#     elif time[t] > 2005:
#         z[t,:] = 2.*np.cos(-2*np.pi*x/250 + 2*np.pi*time[t]/10)
#     elif time[t] > 1995:
#         z[t,:] = 2.*np.cos(2*np.pi*x/250 + 2*np.pi*time[t]/10)
#     else:
#         z[t,:] = 2.*np.cos(-2*np.pi*x/250 + 2*np.pi*time[t]/10)

#z = np.nan * np.zeros((len(time),len(x)))
#X, T = np.meshgrid(x,time)
#z = 2.*np.cos(-2*np.pi*X/250 + 2*np.pi*T/10)

#for t in range(len(time)):
#    z[t,:] = 2.*np.cos(-2*np.pi*x/1000 + 2*np.pi*time[t]/1) + 2*np.cos(-2*np.pi*x/250 + 2*np.pi*time[t]/10)

time = imTime
z = stack
x = ystack


fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(np.flipud((z)), cmap='bwr', extent=[x[0], x[-1], time[0], time[-1]])
ax.set_aspect(15)

alllines = z
xinterp = x
demean = alllines - np.mean(alllines,axis=0)
from scipy.signal import hilbert
#from scipy.fftpack import hilbert
data = (hilbert(demean.T))
#data = (hilbert(demean))


#datareal = np.imag(hilbert(alllines))
#dataimag = np.real(hilbert(alllines))
data = data.T
c = np.matmul(np.conj(data).T,data)/np.shape(data)[0]


import scipy.linalg as la
import numpy.linalg as npla

lamda, loadings = la.eigh(c)

lamda2, loadings2 = npla.eig(c)

ind = np.argsort(lamda[::-1])

lamda[::-1].sort()

loadings = loadings[:,ind]

pcs = np.dot(data, loadings)# / np.sqrt(lamda)
loadings = loadings# * np.sqrt(lamda)
pcsreal = np.real(pcs[:,:])
pcsimag = np.imag(pcs[:,:])
eofreal = np.real(loadings[:,:])
eofimag = np.imag(loadings[:,:])
S = np.power(loadings*np.conj(loadings),0.5) * np.sqrt(lamda)

theta = np.arctan2(eofimag,eofreal)
theta2 = theta*180/np.pi

Rt = np.power(pcs*np.conj(pcs),0.5) / np.sqrt(lamda)

phit = np.arctan2(pcsimag,pcsreal)
phit2 = phit*180/np.pi

mode = 0

fig, ax = plt.subplots(2,2)

ax[0,0].plot(xinterp, S[:,mode],'o')
ax[0,1].plot(xinterp, theta2[:,mode],'o')
ax[1,0].plot(time,Rt[:,mode],'o')
#ax[1,0].plot(time,Rt[mode,:],'o')
ax[1,1].plot(time,phit2[:,mode],'o')
#ax[1,1].plot(time,phit2[mode,:],'o')

mode = 1

fig, ax = plt.subplots(2,2)

ax[0,0].plot(xinterp, S[:,mode],'o')
ax[0,1].plot(xinterp, theta2[:,mode],'o')
ax[1,0].plot(time,Rt[:,mode],'o')
ax[1,1].plot(time,phit2[:,mode],'o')

mode = 2

fig, ax = plt.subplots(2,2)

ax[0,0].plot(xinterp, S[:,mode],'o')
ax[0,1].plot(xinterp, theta2[:,mode],'o')
ax[1,0].plot(time,Rt[:,mode],'o')
ax[1,1].plot(time,phit2[:,mode],'o')

mode = 3

fig, ax = plt.subplots(2,2)

ax[0,0].plot(xinterp, S[:,mode],'o')
ax[0,1].plot(xinterp, theta2[:,mode],'o')
ax[1,0].plot(time,Rt[:,mode],'o')
ax[1,1].plot(time,phit2[:,mode],'o')


totalV = np.sum(lamda)
percentV = lamda / totalV


timeind = np.arange(0,len(time))
#timeind = np.arange(3, 50)
RtSubset = Rt[timeind, :]
phitSubset = phit[timeind, :]
phit2Subset = phit2[timeind, :]
timeSubset = time[timeind]
alllinesSubset = alllines[timeind, :]

eofPred = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred2 = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred3 = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred4 = np.nan * np.ones((np.shape(alllinesSubset)))


for timestep in range(len(timeind)):
    mode = 0
    eofPred[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
            timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    mode = 1
    eofPred2[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
            timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    mode = 2
    eofPred3[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
            timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    mode = 3
    eofPred4[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
            timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])

t1 = 0
t2 = -1

fig, ax = plt.subplots(1,5)
#plt.set_cmap('RdBu')#bwr')
#plt.set_cmap('bwr')

tg, xg = np.meshgrid(time, xinterp)
#xg, tg = np.meshgrid(time, xinterp)
#alllines.T
#eofPred.T
#eofPred2.T
#eofPred3.T
#eofPred4.T


plt0 = ax[0].pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T)#, vmin=-2.5, vmax=2.5)
fig.colorbar(plt0, ax=ax[0], orientation='horizontal')
#ax[0].set_ylim([time[t1], time[t2]])
ax[0].set_title('Hypothetical Morphology')

plt1 = ax[1].pcolor(xg,tg,eofPred.T)#, vmin=-2.5, vmax=2.5)
#ax[1].set_ylim([time[t1], time[t2]])
fig.colorbar(plt1, ax=ax[1], orientation='horizontal')
ax[1].set_title('CEOF1 {:.2f}%'.format(percentV[0]))
ax[1].get_yaxis().set_ticks([])

plt2 = ax[2].pcolor(xg,tg,eofPred2.T)#, vmin=-2.5, vmax=2.5)
#ax[2].set_ylim([time[t1], time[t2]])
fig.colorbar(plt2, ax=ax[2], orientation='horizontal')
ax[2].set_title('CEOF2 {:.2f}%'.format(percentV[1]))
ax[2].get_yaxis().set_ticks([])

plt3 = ax[3].pcolor(xg,tg,eofPred3.T)#, vmin=-2.5, vmax=2.5)
#ax[3].set_ylim([time[t1], time[t2]])
fig.colorbar(plt3, ax=ax[3], orientation='horizontal')
ax[3].set_title('CEOF3 {:.2f}%'.format(percentV[1]))
ax[3].get_yaxis().set_ticks([])

plt4 = ax[4].pcolor(xg,tg,eofPred4.T)#, vmin=-2.5, vmax=2.5)
#ax[4].set_ylim([time[t1], time[t2]])
fig.colorbar(plt4, ax=ax[4], orientation='horizontal')
ax[4].set_title('CEOF4 {:.2f}%'.format(percentV[1]))
ax[4].get_yaxis().set_ticks([])