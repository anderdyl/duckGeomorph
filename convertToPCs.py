import sys
import pickle
from sklearn.decomposition import PCA
import numpy as np
from scipy.interpolate import interp1d
from shapely.geometry import LineString
import scipy.io



def updatePCs(stormNumber):

    with open(r'nourishmentProfileEOFs2.pickle', "rb") as input_file:
        outputEOFs = pickle.load(input_file)

    zInt = outputEOFs['zInt']
    x = outputEOFs['x']
    dist = x
    originalIPCA = PCA(n_components=9)
    origPCs = originalIPCA.fit_transform(zInt)

    dataPredS = scipy.io.loadmat('latestScarpProfile.mat')
    distm = dataPredS['dist']
    predScarpX = dataPredS['predScarpX']
    predScarpZ = dataPredS['predScarpZ']
    zInterp = dataPredS['zInterp']

    newPCs = originalIPCA.transform(zInterp)#surrogateProf.reshape(1,-1))

    #newPCs =  [pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9]
    newPCsProfile = originalIPCA.inverse_transform(newPCs)

    latestPCs = dict()
    latestPCs['newPCs'] = newPCs
    latestPCs['newPCsProfile'] = newPCsProfile
    latestPCs['dist'] = dist
    scipy.io.savemat('latestProfileFit.mat', latestPCs)

    return newPCs



if __name__ == '__main__':

    stormNumber = float(sys.argv[1])

    newPCs = updatePCs(stormNumber)
