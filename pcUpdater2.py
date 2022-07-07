import sys
import pickle
from sklearn.decomposition import PCA
import numpy as np
from scipy.interpolate import interp1d
from shapely.geometry import LineString
import scipy.io


def updatePCs(pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9):

    with open(r'nourishmentProfileEOFs2.pickle', "rb") as input_file:
        outputEOFs = pickle.load(input_file)

    zInt = outputEOFs['zInt']
    x = outputEOFs['x']
    dist = x
    originalIPCA = PCA(n_components=9)
    origPCs = originalIPCA.fit_transform(zInt)
    newPCs =  [pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9]
    newPCsProfile = originalIPCA.inverse_transform(newPCs)

    latestPCs = dict()
    latestPCs['newPCs'] = newPCs
    latestPCs['newPCsProfile'] = newPCsProfile
    latestPCs['dist'] = dist
    scipy.io.savemat('latestProfileFit2.mat', latestPCs)

    return newPCs



if __name__ == '__main__':

    pc1 = float(sys.argv[1])
    pc2 = float(sys.argv[2])
    pc3 = float(sys.argv[3])
    pc4 = float(sys.argv[4])
    pc5 = float(sys.argv[5])
    pc6 = float(sys.argv[6])
    pc7 = float(sys.argv[7])
    pc8 = float(sys.argv[8])
    pc9 = float(sys.argv[9])

    newPCs = updatePCs(pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9)
