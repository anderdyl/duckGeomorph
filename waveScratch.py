


from getdatatestbed import getDataFRF
import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from glob import glob
import os
from scipy.interpolate import interp1d
import sandBarTool.morphLib as mL
#from downloads import pyeemd


wavedir = '/home/dylananderson/projects/dataFRF/'

# Need to sort the files to ensure correct temporal order...
files = os.listdir(wavedir)
files.sort()
files_path = [os.path.join(os.path.abspath(wavedir), x) for x in files]


def get8marray(file):
    waves = Dataset(file)

    waveHs = waves.variables['waveHs'][:]
    waveTp = waves.variables['waveTp'][:]
    waveTm = waves.variables['waveTm'][:]
    waveTm1 = waves.variables['waveTm1'][:]
    waveTm2 = waves.variables['waveTm2'][:]
    waveTpfb = waves.variables['waveTpfb'][:]
    waveMeanDirection = waves.variables['waveMeanDirection'][:]
    waveMeanDirectionPeakFrequency = waves.variables['waveMeanDirectionPeakFrequency'][:]
    wavePeakDirectionPeakFrequency = waves.variables['wavePeakDirectionPeakFrequency'][:]
    wavePrincipleDirection = waves.variables['wavePrincipleDirection'][:]
    time = waves.variables['time'][:]


    output = dict()
    output['waveHs'] = waveHs
    output['waveTp'] = waveTp
    output['waveTm'] = waveTm
    output['waveTm1'] = waveTm1
    output['waveTm2'] = waveTm2
    output['waveTpfb'] = waveTpfb
    output['waveMeanDirection'] = waveMeanDirection
    output['wavePeakDirectionPeakFrequency'] = wavePeakDirectionPeakFrequency
    output['waveMeanDirectionPeakFrequency'] = waveMeanDirectionPeakFrequency
    output['wavePrincipleDirection'] = wavePrincipleDirection
    output['t'] = time

    return output



sept = get8marray(files_path[0])
septTime = [DT.datetime.fromtimestamp(x) for x in sept['t']]
oct = get8marray(files_path[1])
octTime = [DT.datetime.fromtimestamp(x) for x in oct['t']]
nov = get8marray(files_path[2])
novTime = [DT.datetime.fromtimestamp(x) for x in nov['t']]

alltime = septTime + octTime + novTime
waveHs = np.append(np.append(sept['waveHs'],oct['waveHs']),nov['waveHs'])
waveTp = np.append(np.append(sept['waveTp'],oct['waveTp']),nov['waveTp'])
waveTm = np.append(np.append(sept['waveTm'],oct['waveTm']),nov['waveTm'])


waveDirpdpf = np.append(np.append(sept['wavePeakDirectionPeakFrequency'],oct['wavePeakDirectionPeakFrequency']),nov['wavePeakDirectionPeakFrequency'])
waveDirmdpf = np.append(np.append(sept['waveMeanDirectionPeakFrequency'],oct['waveMeanDirectionPeakFrequency']),nov['waveMeanDirectionPeakFrequency'])
waveDirmd =np.append(np.append(sept['waveMeanDirection'],oct['waveMeanDirection']),nov['waveMeanDirection'])

waveNormpdpf = waveDirpdpf-72
waveNormmdpf = waveDirmdpf-72
waveNormmd = waveDirmd-72

#neg = np.where((waveNormmd>180))
#waveNormmd[neg[0]] = waveNormmd-360

LWP_pdpf = 1025*np.square(waveHs)*waveTp*(9.81/(64*np.pi))*np.cos(waveNormpdpf*(np.pi/180))*np.sin(waveNormpdpf*(np.pi/180))
LWP_mdpf = 1025*np.square(waveHs)*waveTp*(9.81/(64*np.pi))*np.cos(waveNormmdpf*(np.pi/180))*np.sin(waveNormmdpf*(np.pi/180))
LWP_md = 1025*np.square(waveHs)*waveTp*(9.81/(64*np.pi))*np.cos(waveNormmd*(np.pi/180))*np.sin(waveNormmd*(np.pi/180))
LWP_means = 1025*np.square(waveHs)*waveTm*(9.81/(64*np.pi))*np.cos(waveNormmd*(np.pi/180))*np.sin(waveNormmd*(np.pi/180))

#cumsum((1025).*((hs.^2).*tp.*(9.81.^2)./(64.*pi)).*cosd(phi).*sind(phi))

fig, ax = plt.subplots(3,1)
ax[0].plot(alltime, waveHs)
ax[1].plot(alltime, waveTp)
ax[1].plot(alltime, waveTm)

ax[2].plot(alltime, waveNormpdpf)
ax[2].plot(alltime, waveNormmdpf)
ax[2].plot(alltime, waveNormmd)

fig = plt.figure(figsize=(10,10))
plt.plot(alltime,(LWP_pdpf))
plt.plot(alltime,(LWP_mdpf))
plt.plot(alltime,(LWP_md))
plt.plot(alltime,(LWP_means))

#plt.plot(alltime,np.cumsum(LWP_pdpf))
#plt.plot(alltime,np.cumsum(LWP_mdpf))
#plt.plot(alltime,np.cumsum(LWP_md))



import glob

argus_files = glob.glob('/mnt/gaia/peeler/argus/argus02bFullFrame/2019/c1/**/*.raw')

argus_times = [x.split('/')[9].split('.')[0] for x in argus_files]

argusDT = [DT.datetime.fromtimestamp(int(x)) for x in argus_times]

fig2 = plt.figure(figsize=(10,10))
plt.plot(alltime,waveHs)

smallWaves = np.where((waveHs < 0.66))
smallTime = [alltime[x] for x in smallWaves[0]]
plt.plot(smallTime,waveHs[smallWaves[0]],'ro')

for x in argusDT:
    plt.plot([x, x],[0, 2],'k-')