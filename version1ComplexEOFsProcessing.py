
from getdatatestbed import getDataFRF
import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from glob import glob
import os
from scipy.interpolate import interp1d
import sandBarTool.morphLib as mL

import pickle
# dbfile = open('sandbarsSouthernTransect_referencedMHHW_reallylongWithNegatives.pickle', 'rb')
dbfile = open('sandbarsSouthernTransect_referencedMHHW_5lineAvgLonger2.pickle', 'rb')

#dbfile = open('sandbarsNorthernTransect_referencedMHHW_reallylong.pickle', 'rb')

data = pickle.load(dbfile)
dbfile.close()
alllinesOG = data['alllines']
xinterpOG = data['xinterp']
time = data['time']

def interpBathy(xSub, zSub, x):
    f = interp1d(xSub, zSub, kind='linear', bounds_error=False)
    newz = f(x)
    newBathy = dict()
    newBathy['x'] = x
    newBathy['z'] = newz
    return newBathy

xinterp = np.arange(0,510,5)
alllines = np.nan*np.ones((len(time),len(xinterp)))

for hh in range(len(time)):
    lowerRes = interpBathy(xinterpOG, alllinesOG[hh,:], xinterp)
    alllines[hh,:] = lowerRes['z']


nanmean = np.nanmean(alllines,axis=0)
finder = np.where(np.isnan(alllines))
alllines[finder] = nanmean[finder[1]]
demean = alllines - np.mean(alllines,axis=0)


plt.style.use('default')
plt.figure()
for xx in range(len(alllines)):
    plt.plot(xinterp, alllines[xx,:], color=[0.7,0.7,0.7,1],linewidth=0.5)
plt.plot(xinterp, np.nanmean(alllines,axis=0),'k-')
plt.plot(xinterp, np.nanmean(alllines,axis=0) + np.nanstd(alllines,axis=0),'k--')
plt.plot(xinterp, np.nanmean(alllines,axis=0) - np.nanstd(alllines,axis=0),'k--')



from scipy.signal import hilbert
data = (hilbert(demean.T))

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



# plt.style.use('default')
# mode = 0
# fig = plt.figure(figsize=(11,6))
# ax1 = plt.subplot2grid((2,2),(0,0),rowspan=1,colspan=1)
# ax1.plot(xinterp, S[:,mode],'.')
# ax1.set_ylabel('amplitude')
# ax1.set_xlabel('cross-shore (m)')
# text = ax1.text(-0.15,1.05, "A)", transform=ax1.transAxes)
#
# ax2 = plt.subplot2grid((2,2),(1,0),rowspan=1,colspan=1)
# ax2.plot(xinterp, theta2[:,mode],'.')
# ax2.set_ylabel('phase')
# ax2.set_xlabel('cross-shore (m)')
# text2 = ax2.text(-0.15,1.05, "B)", transform=ax2.transAxes)
#
# ax3 = plt.subplot2grid((2,2),(0,1),rowspan=1,colspan=1)
# ax3.plot(time,Rt[:,mode],'.')
# ax3.set_ylabel('magnitude')
# ax3.set_xlabel('time')
# text3 = ax3.text(-0.15,1.05, "C)", transform=ax3.transAxes)
#
# ax4 = plt.subplot2grid((2,2),(1,1),rowspan=1,colspan=1)
# ax4.plot(time,phit2[:,mode],'.')
# ax4.set_ylabel('phase')
# ax4.set_xlabel('time')
# text4 = ax4.text(-0.15,1.05, "D)", transform=ax4.transAxes)
#

#
# ax1.set_xlim([120,600])
# ax1.set_title('Cross-shore Surveys')
# ax1.set_ylabel('Years')
# ax1.set_xlabel('cross-shore distance (m)')
# fig.subplots_adjust(right=0.84)
# cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
# cbar = fig.colorbar(plt0, cax=cbar_ax)
# # cbar.set_label('WE (Hs$^2$Tp)')
# cbar.set_label('Deviation from mean profile (m)')




mode = 0

fig, ax = plt.subplots(2,2)

ax[0,0].plot(xinterp, S[:,mode],'o')
ax[0,0].set_ylabel('Spatial Magnitude (m)')
ax[0,0].set_xlabel('Cross-shore (m)')
ax[1,0].plot(xinterp, theta2[:,mode],'o')
ax[1,0].set_ylabel('Spatial Phase (deg)')
ax[1,0].set_xlabel('Cross-shore (m)')
ax[0,1].plot(time,Rt[:,mode],'o')
ax[0,1].set_ylabel('Temporal Magnitude (m)')
ax[0,1].set_xlabel('Time')
ax[1,1].plot(time,phit2[:,mode],'o')
ax[1,1].set_ylabel('Temporal Phase (deg)')
ax[1,1].set_xlabel('Time')



PC1 = Rt[:, mode]*np.sin(phit[:, mode]) + Rt[:, mode]*np.cos(phit[:, mode])

#PC1a = Rt[:, mode]*np.sin(phit[:, mode])
#PC1b = Rt[:, mode]*np.cos(phit[:, mode])
#
# plt.figure()
# plt.plot(time,PC1,'.')
# plt.plot(time,PC1a,'.')
# plt.plot(time,PC1b,'.')


totalV = np.sum(lamda)
percentV = lamda / totalV

ztemp = 0*np.ones(len(xinterp),)
timestep = 200
for mode in range(2):
    # ztemp = ztemp + Rt[timestep,mode]*np.sin(phit[timestep,mode]) * S[:,mode]*np.sin(theta[:,mode]) + Rt[timestep,mode]*np.cos(phit[timestep,mode]) * S[:,mode]*np.cos(theta[:,mode])
    ztemp = ztemp + Rt[timestep,mode]*S[:,mode]*np.cos(phit[timestep,mode] - theta[:,mode])

#Plotting a comparison
fig2 = plt.figure()
#plt.plot(xinterp,np.mean(alllines,axis=0)+ztemp+ztemp1+ztemp2+ztemp3+ztemp4+ztemp5+ztemp6+ztemp7+ztemp8)
plt.plot(xinterp,np.mean(alllines,axis=0)+ztemp)
#plt.plot(xinterp,np.mean(alllines,axis=0)+ztemp2)
plt.plot(xinterp,np.mean(alllines,axis=0))
plt.plot(xinterp,alllines[timestep,:])

#plt.plot(xinterp,ztemp)
#plt.plot(xinterp,ztemp1)
#plt.plot(xinterp,ztemp2)
#plt.plot(xinterp,ztemp+ztemp1+ztemp2)


def P2R(radii, angles):
    return radii * np.exp(1j*angles)

def R2P(x):
    return np.abs(x), np.angle(x)



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
    # eofPred[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
    #         timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    eofPred[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
    mode = 1
    # eofPred2[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
    #         timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    eofPred2[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])

    mode = 2
    # eofPred3[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
    #         timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    eofPred3[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])

    mode = 3
    # eofPred4[timestep, :] = RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
    #         timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])
    eofPred4[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])

    # ztemp = 0 * np.ones(len(xinterp), )
    # for mode in range(8):
    #     ztemp = ztemp + RtSubset[timestep, mode] * np.sin(phitSubset[timestep, mode]) * S[:, mode] * np.sin(theta[:, mode]) + RtSubset[
    #         timestep, mode] * np.cos(phitSubset[timestep, mode]) * S[:, mode] * np.cos(theta[:, mode])



t1 = 0
t2 = 300
#
# plt.style.use('default')
#
# fig = plt.figure(figsize=(7,11))
# ax1 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
# plt.set_cmap('RdBu_r')
#
# tg, xg = np.meshgrid(time, xinterp)
# plt0 = ax1.pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1, vmax=1)
# # cb = fig.colorbar(plt0, ax=ax1, fraction=0.045, pad=0.1)
# ax1.set_ylim([time[t1], time[t2]])
# ax1.set_xlim([120,600])
# ax1.set_title('Cross-shore Surveys')
# ax1.set_ylabel('Years')
# ax1.set_xlabel('cross-shore distance (m)')
# fig.subplots_adjust(right=0.84)
# cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
# cbar = fig.colorbar(plt0, cax=cbar_ax)
# # cbar.set_label('WE (Hs$^2$Tp)')
# cbar.set_label('Deviation from mean profile (m)')



t1 =495
t2 = -24
#
# plt.style.use('default')
#
# fig = plt.figure()
# ax1 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
# #plt.set_cmap('RdBu')#bwr')
# plt.set_cmap('bwr')
#
# tg, xg = np.meshgrid(time, xinterp)
# plt0 = ax1.pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.6, vmax=1.6)
# fig.colorbar(plt0, ax=ax1, orientation='horizontal')
# ax1.set_ylim([time[t1], time[t2]])
# ax1.set_title('Cross-shore Surveys (deviation from mean profile)')
#
#
# fig, ax = plt.subplots(1,5)
# #plt.set_cmap('RdBu')#bwr')
# plt.set_cmap('RdBu_r')
#
# tg, xg = np.meshgrid(time, xinterp)
# plt0 = ax[0].pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
# fig.colorbar(plt0, ax=ax[0], orientation='horizontal')
# ax[0].set_ylim([time[t1], time[t2]])
# ax[0].set_title('Surveys (dev.)')
#
# plt1 = ax[1].pcolor(xg,tg,eofPred.T, vmin=-1.05, vmax=1.05)
# ax[1].set_ylim([time[t1], time[t2]])
# fig.colorbar(plt1, ax=ax[1], orientation='horizontal')
# ax[1].set_title('CEOF1 {:.2f}'.format(percentV[0]))
# ax[1].get_yaxis().set_ticks([])
#
# plt2 = ax[2].pcolor(xg,tg,eofPred2.T, vmin=-.85, vmax=.85)
# ax[2].set_ylim([time[t1], time[t2]])
# fig.colorbar(plt2, ax=ax[2], orientation='horizontal')
# ax[2].set_title('CEOF2 {:.2f}'.format(percentV[1]))
# ax[2].get_yaxis().set_ticks([])
#
# plt3 = ax[3].pcolor(xg,tg,eofPred3.T, vmin=-.65, vmax=.65)
# ax[3].set_ylim([time[t1], time[t2]])
# ax[3].set_title('CEOF3 {:.2f}'.format(percentV[2]))
# ax[3].get_yaxis().set_ticks([])
#
# fig.colorbar(plt3, ax=ax[3], orientation='horizontal')
#
# plt4 = ax[4].pcolor(xg,tg,eofPred4.T, vmin=-.65, vmax=.65)
# ax[4].set_ylim([time[t1], time[t2]])
# ax[4].set_title('CEOF4 {:.2f}'.format(percentV[3]))
# ax[4].get_yaxis().set_ticks([])
# fig.colorbar(plt4, ax=ax[4], orientation='horizontal')
#


t1 = 0
t2 = -1
fig, ax = plt.subplots(1,5)
#plt.set_cmap('RdBu')#bwr')
plt.set_cmap('RdBu_r')

tg, xg = np.meshgrid(time, xinterp)
plt0 = ax[0].pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
fig.colorbar(plt0, ax=ax[0], orientation='horizontal')
ax[0].set_ylim([time[t1], time[t2]])
ax[0].set_title('Surveys (dev.)')

plt1 = ax[1].pcolor(xg,tg,eofPred.T, vmin=-.75, vmax=0.75)
ax[1].set_ylim([time[t1], time[t2]])
fig.colorbar(plt1, ax=ax[1], orientation='horizontal')
ax[1].set_title('CEOF1 {:.2f}'.format(percentV[0]))
ax[1].get_yaxis().set_ticks([])

plt2 = ax[2].pcolor(xg,tg,eofPred2.T, vmin=-.85, vmax=.85)
ax[2].set_ylim([time[t1], time[t2]])
fig.colorbar(plt2, ax=ax[2], orientation='horizontal')
ax[2].set_title('CEOF2 {:.2f}'.format(percentV[1]))
ax[2].get_yaxis().set_ticks([])

plt3 = ax[3].pcolor(xg,tg,eofPred3.T, vmin=-.65, vmax=.65)
ax[3].set_ylim([time[t1], time[t2]])
ax[3].set_title('CEOF3 {:.2f}'.format(percentV[2]))
ax[3].get_yaxis().set_ticks([])

fig.colorbar(plt3, ax=ax[3], orientation='horizontal')

plt4 = ax[4].pcolor(xg,tg,eofPred4.T, vmin=-.65, vmax=.65)
ax[4].set_ylim([time[t1], time[t2]])
ax[4].set_title('CEOF4 {:.2f}'.format(percentV[3]))
ax[4].get_yaxis().set_ticks([])
fig.colorbar(plt4, ax=ax[4], orientation='horizontal')

plt.tight_layout(pad=0.5)

plt.show()


totalV = np.sum(lamda)
percentV = lamda / totalV
import datetime
def datetime2matlabdn(dt):
   ord = dt.toordinal()
   mdn = dt + datetime.timedelta(days = 366)
   frac = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   return mdn.toordinal() + frac

matlabTime = np.zeros((len(time),))
for qq in range(len(time)):
    matlabTime[qq] = datetime2matlabdn(time[qq])
#
# morphoPickle = 'ceofsSouthernTransectLatestLonger.pickle'
# output = {}
# output['time'] = time
# output['alllines'] = alllines
# output['xinterp'] = xinterp
# output['S'] = S
# output['Rt'] = Rt
# output['thetaRadians'] = theta
# output['thetaDegrees'] = theta2
# output['phiRadian'] = phit
# output['phiDegrees'] = phit2
# output['totalV'] = totalV
# output['lamda'] = lamda
# output['percentV'] = percentV
#
#
#
# import pickle
# with open(morphoPickle,'wb') as f:
#     pickle.dump(output, f)


import scipy.io
output = dict()
output['alllines'] = alllines
output['xinterp'] = xinterp
output['time'] = matlabTime
output['S'] = S
output['Rt'] = Rt
output['thetaRadians'] = theta
output['thetaDegrees'] = theta2
output['phiRadian'] = phit
output['phiDegrees'] = phit2
output['totalV'] = totalV
output['lamda'] = lamda
output['percentV'] = percentV
scipy.io.savemat('ceofsSouthernTransectLatestLonger2.mat',output)


from scipy.stats import linregress

recon1 = eofPred.T
recon2 = eofPred2.T
recon3 = eofPred3.T
recon4 = eofPred4.T

cumrecon = recon1+recon2
surveys = (alllines-np.mean(alllines, axis=0)).T# - recon4

plt.figure()
plt.plot(surveys[0,:],recon1[0,:],'o')

rsquared = np.nan * np.ones((len(xinterp),))
rvalues = np.nan * np.ones((len(xinterp),))
slopes = np.nan * np.ones((len(xinterp),))
intercepts = np.nan * np.ones((len(xinterp),))
rsquared2 = np.nan * np.ones((len(xinterp),))
rvalues2 = np.nan * np.ones((len(xinterp),))
slopes2 = np.nan * np.ones((len(xinterp),))
intercepts2 = np.nan * np.ones((len(xinterp),))

for hh in range(len(xinterp)):
    slope, intercept, r, p, se = linregress(surveys[hh,:], cumrecon[hh,:])
    rsquared[hh] = np.square(r)
    rvalues[hh] = r
    slopes[hh] = slope
    intercepts[hh] = intercept
    slope2, intercept2, r2, p2, se2 = linregress(surveys[hh,:], recon1[hh,:])
    rsquared2[hh] = np.square(r2)
    rvalues2[hh] = r2
    slopes2[hh] = slope2
    intercepts2[hh] = intercept2

plt.figure()
plt.plot(xinterp,rsquared)
plt.plot(xinterp,rsquared2)










