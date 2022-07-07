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



def updateProfile(stormNumber):

    dataPredS = scipy.io.loadmat('latestScarpFit.mat')
    distx = dataPredS['distx']
    a = dataPredS['a']
    b = dataPredS['b']

    with open(r'nourishmentProfileEOFs2.pickle', "rb") as input_file:
        outputEOFs = pickle.load(input_file)

    zInt = outputEOFs['zInt']
    x = outputEOFs['x']
    dist = x
    originalIPCA = PCA(n_components=9)
    zeroLine = np.zeros((np.shape(dist)))

    origPCs = originalIPCA.fit_transform(zInt)
    if stormNumber == 1:
        dataPred = scipy.io.loadmat('startProfile.mat')
        startProfilePCs = dataPred['startingPCs']
        pcProfile = originalIPCA.inverse_transform(startProfilePCs)
        scarpProfile = originalIPCA.inverse_transform(startProfilePCs)

    else:
        dataPred = scipy.io.loadmat('latestStormFit.mat')
        # newPCs = dataPred['newPCs']
        pcProfile = dataPred['pcProfile']
        scarpProfile = dataPred['scarpProfile']



    f = zeroLine
    g = pcProfile[0]
    first_line = LineString(np.column_stack((dist, f)))
    second_line = LineString(np.column_stack((dist, g)))
    intersection = first_line.intersection(second_line)

    if intersection.geom_type == 'Multi3Point':
        print('shit, multipoint')
    elif intersection.geom_type == 'Point':
        x, y = intersection.xy
        zeroCrossingX = x[0]

    f2 = zeroLine
    g2 = scarpProfile[0]
    first_line2 = LineString(np.column_stack((dist, f2)))
    second_line2 = LineString(np.column_stack((dist, g2)))
    intersection2 = first_line2.intersection(second_line2)

    if intersection.geom_type == 'Point':
        x2, y2 = intersection2.xy
        print('x2 = {}, and y2 = {}'.format(x2,y2))
        if x2 > x:
            zeroCrossingX = x2[0]
            print('choosing EOF shoreline: {}'.format(stormNumber))
        else:
            zeroCrossingX = x[0]
            print('zeroCrossingX = {}'.format(zeroCrossingX))

    # Lets make a new line
    ogDuneInd = np.where((dist < (zeroCrossingX - distx)))

    scarpXDune = dist[ogDuneInd[1]]
    scarpZDune = pcProfile[ogDuneInd]

    # finding the points between the scarp and the zero point...
    # xIndex = np.where((dist <= zeroCrossingX) & (dist > (zeroCrossingX - dist1[0][profIndex])))

    # find the points between the scarp and the zero point + 100 m
    xIndex = np.where((dist <= zeroCrossingX + 100) & (dist > (zeroCrossingX - distx)))

    z2 = a * np.power(dist[xIndex[1]] - dist[xIndex[1][0]], b)
    newx2 = dist[xIndex[1]] - dist[xIndex[1][0]]
    scarpXBeach = newx2 + zeroCrossingX - distx
    scarpZBeach = z2 + pcProfile[xIndex[1][0]]

    # # finding negatives and making them zzero
    # findNegatives = np.where((scarpZBeach < 0))
    # scarpZBeach[findNegatives] = scarpZBeach[findNegatives]*0

    diff = scarpZBeach - pcProfile[xIndex]
    keepers = np.where(diff[0] <= 0)

    scarpXBeach = scarpXBeach[0][keepers[0]]
    scarpZBeach = scarpZBeach[0][keepers[0]]

    if scarpZBeach[-1] < -2:
        print('the bottom of the scarp is in -2 territory')
        alter = np.where((scarpZBeach < 0))
        below = np.where((scarpProfile < -2))
        xlinear = [scarpXBeach[alter[0][0]], dist[below[1][0]]]
        zlinear = [0, -2]
        flinear = interp1d(xlinear, zlinear)
        offshore = np.where((dist > scarpXBeach[alter[0][0]]))
        if offshore[0][0] == below[1][0]:
            xflat = dist[offshore[0][0]:int(offshore[0][0]+1)]
            zflat = flinear(xflat)
        else:
            xflat = dist[offshore[0][0]:below[1][0]]
            zflat = flinear(xflat)

        keepers2 = np.where((scarpZBeach > 0))
        scarpXBeach = scarpXBeach[keepers2[0]]
        scarpZBeach = scarpZBeach[keepers2[0]]

        ogSubaqueous = np.where((dist > xflat[-1]))
        scarpXSubaqueous = dist[ogSubaqueous[0]]
        scarpZSubaqueous = scarpProfile[0][ogSubaqueous[0]]

        predScarp1X = np.hstack((scarpXDune, scarpXBeach))
        predScarp15X = np.hstack((predScarp1X, xflat))
        predScarp2X = np.hstack((predScarp15X, scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune, scarpZBeach))
        predScarp15Z = np.hstack((predScarp1Z, zflat))
        predScarp2Z = np.hstack((predScarp15Z, scarpZSubaqueous))

    else:
        ogSubaqueous = np.where((dist > scarpXBeach[-1]))
        scarpXSubaqueous = dist[ogSubaqueous[0]]
        scarpZSubaqueous = scarpProfile[0][ogSubaqueous[0]]
        predScarp1X = np.hstack((scarpXDune, scarpXBeach))
        predScarp2X = np.hstack((predScarp1X, scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune, scarpZBeach))
        predScarp2Z = np.hstack((predScarp1Z, scarpZSubaqueous))

    finterp = interp1d(predScarp2X, predScarp2Z, 'cubic')

    surrogateProf = finterp(dist)

    smoothed = moving_average(surrogateProf, 7)
    extended = np.hstack((surrogateProf[0:3], smoothed))
    extended2 = np.hstack((extended, surrogateProf[-4:-1]))
    # surrogates2[hh,:] = extended2 #
    # surrogates2[hh, :] = surrogateProf
    newPCs = originalIPCA.transform(surrogateProf.reshape(1,-1))
    newScarpProfile = surrogateProf
    newPCsProfile = originalIPCA.inverse_transform(newPCs)

    latestPCs = dict()
    latestPCs['newPCs'] = newPCs
    latestPCs['pcProfile'] = newPCsProfile
    latestPCs['scarpProfile'] = surrogateProf
    latestPCs['distance'] = dist
    scipy.io.savemat('latestStormFit.mat', latestPCs)

    # distanceX = dict()
    # distanceX['distance'] = dist
    # scipy.io.savemat('distanceX.mat', distanceX)

    return newPCs



if __name__ == '__main__':

    stormNumber = float(sys.argv[1])
    newPCs = updateProfile(stormNumber)
