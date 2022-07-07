import sys
import pickle
from sklearn.decomposition import PCA
import numpy as np
from scipy.interpolate import interp1d
from shapely.geometry import LineString
import scipy.io


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n



def scarpOnProfile(pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,a1,b1,d1):

    # d1 = 34.5
    # a1 = -0.54
    # b1 = 0.447
    # pc1 = -7.1007
    # pc2 = -4.6730
    # pc3 = 3.0054
    # pc4 = 1.3754
    # pc5 = 3.6793
    # pc6 = 2.7119
    # pc7 = -0.6917
    # pc8 = 0.9720
    # pc9 = 0.4590
    # pc1 = -18.0095
    # pc2 = 5.7950
    # pc3 = 3.4899
    # pc4 =  - 1.0298
    # pc5 = 6.6321
    # pc6 =  - 1.3476
    # pc7 =  - 1.1262
    # pc8 =  - 0.3094
    # pc9 =  - 0.4824
    # d1 = 36.02
    # a1 = -2.065
    # b1 = 0.22
    distx = d1
    a = a1
    b = b1

    with open(r'nourishmentProfileEOFs2.pickle', "rb") as input_file:
        outputEOFs = pickle.load(input_file)

    zInt = outputEOFs['zInt']
    x = outputEOFs['x']
    dist = x
    originalIPCA = PCA(n_components=9)
    zeroLine = np.zeros((np.shape(dist)))

    origPCs = originalIPCA.fit_transform(zInt-np.mean(zInt,axis=0))

    startProfilePCs =  [pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9]
    pcProfile = originalIPCA.inverse_transform(startProfilePCs)+np.mean(zInt,axis=0)

    f = zeroLine
    g = pcProfile
    first_line = LineString(np.column_stack((dist, f)))
    second_line = LineString(np.column_stack((dist, g)))
    intersection = first_line.intersection(second_line)

    if intersection.geom_type == 'Multi3Point':
        print('shit, multipoint')
    elif intersection.geom_type == 'Point':
        x, y = intersection.xy
        zeroCrossingX = x[0]

    # Lets make a new line
    ogDuneInd = np.where((dist < (zeroCrossingX - distx)))

    scarpXDune = dist[ogDuneInd[0]]
    scarpZDune = pcProfile[ogDuneInd]

    # finding the points between the scarp and the zero point...
    # xIndex = np.where((dist <= zeroCrossingX) & (dist > (zeroCrossingX - dist1[0][profIndex])))

    # find the points between the scarp and the zero point + 100 m
    xIndex = np.where((dist <= zeroCrossingX + 100) & (dist > (zeroCrossingX - distx)))

    z2 = a * np.power(dist[xIndex[0]] - dist[xIndex[0][0]], b)
    newx2 = dist[xIndex[0]] - dist[xIndex[0][0]]
    scarpXBeach = newx2 + zeroCrossingX - distx
    scarpZBeach = z2 + pcProfile[xIndex[0][0]]

    # # finding negatives and making them zzero
    # findNegatives = np.where((scarpZBeach < 0))
    # scarpZBeach[findNegatives] = scarpZBeach[findNegatives]*0

    diff = scarpZBeach - pcProfile[xIndex]
    keepers = np.where(diff <= 0)

    scarpXBeach = scarpXBeach[keepers]
    scarpZBeach = scarpZBeach[keepers]

    # if scarpZBeach[-1] < -2:
    #     print('the bottom of the scarp is in -2 territory')
    #
    #     if diff[-1] < 0:
    #
    #         # shortenInd = np.where(np.diff(np.where(diff < 0))[0] > 1)
    #         # scarpXBeach = scarpXBeach[0:(np.where(diff<0)[0][shortenInd[0][0]]+2)]
    #         # scarpZBeach = scarpZBeach[0:(np.where(diff<0)[0][shortenInd[0][0]]+2)]
    #         scarpXBeach = scarpXBeach[0:-30]
    #         scarpZBeach = scarpZBeach[0:-30]
    #
    #     alter = np.where((scarpZBeach < 0))
    #     below = np.where((pcProfile < -2))
    #     xlinear = [scarpXBeach[alter[0][0]], dist[below[0][0]]]
    #     zlinear = [0, -2]
    #     flinear = interp1d(xlinear, zlinear)
    #     offshore = np.where((dist > scarpXBeach[alter[0][0]]))
    #     if offshore[0][0] == below[0][0]:
    #         xflat = dist[offshore[0][0]:int(offshore[0][0]+1)]
    #         zflat = flinear(xflat)
    #     else:
    #         xflat = dist[offshore[0][0]:below[0][0]]
    #         zflat = flinear(xflat)
    #
    #     keepers2 = np.where((scarpZBeach > 0))
    #     scarpXBeach = scarpXBeach[keepers2[0]]
    #     scarpZBeach = scarpZBeach[keepers2[0]]
    #
    #     ogSubaqueous = np.where((dist > xflat[-1]))
    #     scarpXSubaqueous = dist[ogSubaqueous[0]]
    #     scarpZSubaqueous = pcProfile[ogSubaqueous[0]]
    #
    #     predScarp1X = np.hstack((scarpXDune, scarpXBeach))
    #     predScarp15X = np.hstack((predScarp1X, xflat))
    #     predScarp2X = np.hstack((predScarp15X, scarpXSubaqueous))
    #     predScarp1Z = np.hstack((scarpZDune, scarpZBeach))
    #     predScarp15Z = np.hstack((predScarp1Z, zflat))
    #     predScarp2Z = np.hstack((predScarp15Z, scarpZSubaqueous))
    #
    # else:
    #     ogSubaqueous = np.where((dist > scarpXBeach[-1]))
    #     scarpXSubaqueous = dist[ogSubaqueous[0]]
    #     scarpZSubaqueous = pcProfile[ogSubaqueous[0]]
    #     predScarp1X = np.hstack((scarpXDune, scarpXBeach))
    #     predScarp2X = np.hstack((predScarp1X, scarpXSubaqueous))
    #     predScarp1Z = np.hstack((scarpZDune, scarpZBeach))
    #     predScarp2Z = np.hstack((predScarp1Z, scarpZSubaqueous))

    ogSubaqueous = np.where((dist > scarpXBeach[-1]))
    scarpXSubaqueous = dist[ogSubaqueous[0]]
    scarpZSubaqueous = pcProfile[ogSubaqueous[0]]
    predScarp1X = np.hstack((scarpXDune, scarpXBeach))
    predScarp2X = np.hstack((predScarp1X, scarpXSubaqueous))
    predScarp1Z = np.hstack((scarpZDune, scarpZBeach))
    predScarp2Z = np.hstack((predScarp1Z, scarpZSubaqueous))

    allNegs = np.where(predScarp2Z < 0)
    predScarp2Z[allNegs] = 0 * np.ones((len(allNegs), 1))

    finterp = interp1d(predScarp2X, predScarp2Z, 'cubic')

    surrogateProf = finterp(dist)

    smoothed = moving_average(surrogateProf, 7)
    extended = np.hstack((surrogateProf[0:3], smoothed))
    extended2 = np.hstack((extended, surrogateProf[-4:-1]))
    # surrogates2[hh,:] = extended2 #
    # surrogates2[hh, :] = surrogateProf
    newPCs = originalIPCA.transform(surrogateProf.reshape(1,-1))
    newScarpProfile = surrogateProf
    newPCsProfile = originalIPCA.inverse_transform(newPCs)+np.mean(zInt,axis=0)

    latestPCs = dict()
    latestPCs['newPCs'] = newPCs
    latestPCs['pcProfile'] = newPCsProfile
    latestPCs['oldPcProfile'] = pcProfile
    latestPCs['scarpProfile'] = surrogateProf
    latestPCs['distance'] = dist
    scipy.io.savemat('latestStormFit.mat', latestPCs)

    # distanceX = dict()
    # distanceX['distance'] = dist
    # scipy.io.savemat('distanceX.mat', distanceX)

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

    a1 = float(sys.argv[10])
    b1 = float(sys.argv[11])
    d1 = float(sys.argv[12])

    newPCs = scarpOnProfile(pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,a1,b1,d1)


