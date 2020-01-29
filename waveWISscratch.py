

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
    waveTm = waves.variables['waveTm'][:]
    waveTm1 = waves.variables['waveTm1'][:]
    waveTm2 = waves.variables['waveTm2'][:]
    waveMeanDirection = waves.variables['waveMeanDirection'][:]
    waveMeanDirectionSwell = waves.variables['waveMeanDirectionSwell'][:]
    time = waves.variables['time'][:]


    output = dict()
    output['waveHs'] = waveHs
    output['waveTp'] = waveTp
    output['waveTm'] = waveTm
    output['waveTm1'] = waveTm1
    output['waveTm2'] = waveTm2
    output['waveMeanDirection'] = waveMeanDirection
    output['waveMeanDirectionSwell'] = waveMeanDirectionSwell
    output['t'] = time

    return output


Hs = []
Tm = []
Dm = []
time = []
for i in files_path:
    waves = getWIS(i)
    Hs = np.append(Hs,waves['waveHs'])
    Tm = np.append(Tm,waves['waveTm'])
    Dm = np.append(Dm,waves['waveMeanDirection'])
    time = np.append(time,waves['t'])



waveNorm = Dm - 72
neg = np.where((waveNorm > 180))
waveNorm[neg[0]] = waveNorm[neg[0]]-360
offpos = np.where((waveNorm>90))
offneg = np.where((waveNorm<-90))
waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
waveNorm[offneg[0]] = waveNorm[offneg[0]]*0



LWP_means = 1025*np.square(Hs)*Tm*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
t = [DT.datetime.fromtimestamp(x) for x in time]

fig = plt.figure(figsize=(10,10))
plt.plot(t,np.cumsum(LWP_means))