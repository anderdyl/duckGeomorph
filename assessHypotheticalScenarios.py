import sys
from sklearn.decomposition import PCA
from scipy.interpolate import interp1d
from shapely.geometry import LineString
from shapely.geometry import Polygon
import scipy.io
import numpy as np
import pickle
from scipy.optimize import curve_fit
from scipy.stats.kde import gaussian_kde
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.ticker import FuncFormatter
from numpy.random import seed
from numpy.random import rand
# seed random number generator
seed(1)

import time

startTime = time.time()


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def objective(t,a,b,c):
    return a*t*t + b*t + c

def objective2(t,a,b):
    return a*t + b

#
# with open(r'hypotheticalTrainingConditions2.pickle', "rb") as input_file:
#     outputData = pickle.load(input_file)
#
# dataN = outputData['dataN']
# dataT = outputData['dataT']
# dataS = outputData['dataS']
#
# M2amp = outputData['M2amp']
# S2amp = outputData['S2amp']
# K1amp = outputData['K1amp']
# N2amp = outputData['N2amp']
#
# M2period = outputData['M2period']
# N2period = outputData['N2period']
# K1period = outputData['K1period']
# S2period = outputData['S2period']
#
# M2speed = outputData['M2speed']
# N2speed = outputData['N2speed']
# K1speed = outputData['K1speed']
# S2speed = outputData['S2speed']
#
# PCs1 = outputData['PCs1']
# PCs2 = outputData['PCs2']
# PCs3 = outputData['PCs3']
# PCs4 = outputData['PCs4']
# PCs5 = outputData['PCs5']
# PCs6 = outputData['PCs6']
# PCs7 = outputData['PCs7']
# PCs8 = outputData['PCs8']
# PCs9 = outputData['PCs9']
#
#
# ht1 = np.arange(2000,2045,1/8766)
# ht2 = np.arange(2045,2070,2/8766)
# ht = np.hstack((ht1,ht2))
# time = np.array([2020,2025,2030,2035,2040,2045,2050,2055,2060,2065,2070,2075,2080,2085,2090,2095,2100])
# USACE = np.array([0.025,0.057,0.09,0.125,0.161,0.199,0.238,0.278,0.32,0.363,0.407,0.453,0.5,0.548,0.598,0.649,0.701])
# USACElow = np.array([0.004,0.027,0.051,0.075,0.099,0.123,0.147,0.171,0.194,0.218,0.242,0.266,0.29,0.314,0.337,0.361,0.385])
# popt, _ = curve_fit(objective,time,USACE)
# popt2, _ = curve_fit(objective,time,USACElow)
# a,b,c = popt
# a2,b2,c2 = popt2
# x_new = np.arange(2000,2070,1/8766)
# usaceLow = objective(ht,a,b,c)
# usaceMed = objective(ht,a2,b2,c2)
#
# slrNumbers = np.floor(rand(len(dataN))*len(usaceMed))
#
# dataSLR = np.empty((len(dataN),1))
# for ii in range(len(dataN)):
#     top = usaceMed[int(slrNumbers[ii])]
#     bot = usaceLow[int(slrNumbers[ii])]
#
#     slrTemp = np.random.uniform(bot,top)
#     dataSLR[ii] = slrTemp
#
# dataNplusSLR = np.hstack((dataN,dataSLR))
# dataNplusSLRplusS = np.hstack((dataNplusSLR,dataS))
# dataNplusSLRplusSplusT = np.hstack((dataNplusSLRplusS,dataT))
#
# allHypos = dataNplusSLRplusSplusT
#
# outMat = dict()
# outMat['allHypos'] = allHypos
#
# scipy.io.savemat('allHypos.mat',outMat)


dataPred = scipy.io.loadmat('initialHyposTest.mat')
distx = dataPred['distx']
a = dataPred['a']
b = dataPred['b']
num = dataPred['num']
designParams2 = dataPred['designParams2']


numberOfSamples = 1000000

a = a[0:numberOfSamples,0]
b = b[0:numberOfSamples,0]
distx = distx[0:numberOfSamples,0]

dataS = designParams2

with open(r'nourishmentProfileEOFs2.pickle', "rb") as input_file:
    outputEOFs = pickle.load(input_file)

zInt = outputEOFs['zInt']
x = outputEOFs['x']
dist = x
originalIPCA = PCA(n_components=9)
origPCs = originalIPCA.fit_transform(zInt)

zeroLine = np.zeros((np.shape(dist)))
actualArea = []
initialShoreline = []
for ii in range(len(a)):

    if np.remainder(ii,20) == 0:
        print('done with {} profile areas'.format(ii))

    distxLine = distx[ii]
    aLine = a[ii]
    bLine = b[ii]

    startProfilePCs = designParams2[ii,0:9]
    pcProfile = originalIPCA.inverse_transform(startProfilePCs)
    scarpProfile = originalIPCA.inverse_transform(startProfilePCs)

    f = zeroLine
    g = pcProfile#[0]
    first_line = LineString(np.column_stack((dist, f)))
    second_line = LineString(np.column_stack((dist, g)))
    intersection = first_line.intersection(second_line)

    if intersection.geom_type == 'Multi3Point':
        print('shit, multipoint')
    elif intersection.geom_type == 'Point':
        x, y = intersection.xy
        zeroCrossingX = x[0]

    initialShoreline.append(zeroCrossingX)

    f2 = zeroLine
    g2 = scarpProfile#[0]
    first_line2 = LineString(np.column_stack((dist, f2)))
    second_line2 = LineString(np.column_stack((dist, g2)))
    intersection2 = first_line2.intersection(second_line2)

    if intersection.geom_type == 'Point':
        x2, y2 = intersection2.xy
        #print('x2 = {}, and y2 = {}'.format(x2,y2))
        if x2 > x:
            zeroCrossingX = x2[0]
            #print('choosing EOF shoreline: {}'.format(ii))
        else:
            zeroCrossingX = x[0]
            #print('zeroCrossingX = {}'.format(zeroCrossingX))

    # Lets make a new line
    ogDuneInd = np.where((dist < (zeroCrossingX - distxLine)))

    scarpXDune = dist[ogDuneInd[0]]
    scarpZDune = pcProfile[ogDuneInd[0]]

    # finding the points between the scarp and the zero point...
    # xIndex = np.where((dist <= zeroCrossingX) & (dist > (zeroCrossingX - dist1[0][profIndex])))

    # find the points between the scarp and the zero point + 100 m
    xIndex = np.where((dist <= zeroCrossingX + 70) & (dist > (zeroCrossingX - distxLine)))

    z2 = aLine * np.power(dist[xIndex[0]] - dist[xIndex[0][0]], bLine)
    newx2 = dist[xIndex[0]] - dist[xIndex[0][0]]
    scarpXBeach = newx2 + zeroCrossingX - distxLine
    scarpZBeach = z2 + pcProfile[xIndex[0][0]]

    # # finding negatives and making them zzero
    # findNegatives = np.where((scarpZBeach < 0))
    # scarpZBeach[findNegatives] = scarpZBeach[findNegatives]*0

    diff = scarpZBeach - pcProfile[xIndex]
    keepers = np.where(diff <= 0)

    scarpXBeach10 = scarpXBeach[keepers[0]]
    scarpZBeach10 = scarpZBeach[keepers[0]]

    bottom = scarpZBeach10[-1]
    if bottom < -2.25:
        #print('the bottom of the scarp is in -2 territory')
        alter = np.where((scarpZBeach10 < 0))
        below = np.where((scarpProfile < -2))
        # below = np.where((scarpProfile < -2))

        xlinear = [scarpXBeach10[alter[0][0]], dist[below[0][0]]]
        zlinear = [0, -2]
        flinear = interp1d(xlinear, zlinear)
        offshore = np.where((dist > scarpXBeach10[alter[0][0]]))
        if offshore[0][0] == below[0][0]:
            xflat = dist[offshore[0][0]:int(offshore[0][0]+1)]
            zflat = flinear(xflat)
        else:
            xflat = dist[offshore[0][0]:below[0][0]]
            zflat = flinear(xflat)

        keepers2 = np.where((scarpZBeach10 > 0))
        scarpXBeach20 = scarpXBeach10[keepers2[0]]
        scarpZBeach20 = scarpZBeach10[keepers2[0]]

        ogSubaqueous = np.where((dist > xflat[-1]))
        scarpXSubaqueous = dist[ogSubaqueous[0]]
        scarpZSubaqueous = scarpProfile[ogSubaqueous[0]]

        predScarp1X = np.hstack((scarpXDune, scarpXBeach20))
        predScarp15X = np.hstack((predScarp1X, xflat))
        predScarp2X = np.hstack((predScarp15X, scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune, scarpZBeach20))
        predScarp15Z = np.hstack((predScarp1Z, zflat))
        predScarp2Z = np.hstack((predScarp15Z, scarpZSubaqueous))

    else:
        ogSubaqueous = np.where((dist > scarpXBeach10[-1]))
        scarpXSubaqueous = dist[ogSubaqueous[0]]
        scarpZSubaqueous = scarpProfile[ogSubaqueous[0]]
        predScarp1X = np.hstack((scarpXDune, scarpXBeach10))
        predScarp2X = np.hstack((predScarp1X, scarpXSubaqueous))
        predScarp1Z = np.hstack((scarpZDune, scarpZBeach10))
        predScarp2Z = np.hstack((predScarp1Z, scarpZSubaqueous))

    finterp = interp1d(predScarp2X, predScarp2Z, 'cubic')

    surrogateProf = finterp(dist)



    # Difference between the XBeach volumes
    x_y_curveTest1 = list()
    x_y_curveTest2 = list()
    z1 = g
    z2 = surrogateProf
    findNegatives2 = np.where((z2<0))
    z2[findNegatives2] = z2[findNegatives2]*0
    finder = np.where((z1<(z2-.15)))

    for n in range(finder[0][0]):
        x_y_curveTest1.append(list((dist[n], z1[n])))
        x_y_curveTest2.append(list((dist[n], z2[n])))

    polygon_points = []  # creates a empty list where we will append the points to create the polygon

    for xyvalue in x_y_curveTest1:
        polygon_points.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 1

    for xyvalue in x_y_curveTest2[::-1]:
        polygon_points.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 2 in the reverse order (from last point to first point)

    for xyvalue in x_y_curveTest1[0:1]:
        polygon_points.append([xyvalue[0], xyvalue[1]])  # append the first point in curve 1 again, to it "closes" the polygon

    polygon = Polygon(polygon_points)
    area = polygon.area
    actualArea.append(area)
    #(area)
    #p5.text(zeroCrossing[0][profIndex]-25,3.5,'{} m^2'.format(np.round(area*10)/10))
    #fdu = integrate.simps()

    #
    # # Difference between the XBeach volumes
    # x_y_curveTest12 = list()
    # x_y_curveTest22 = list()
    # z12 = preProfiles[profIndex, :]
    # f = interp1d(dist,z12,kind='cubic')
    # z12new = f(np.ma.filled(predScarp2X))
    # z22 = predScarp2Z
    # finder = np.where((z12new<(z22-.1)))
    # if finder[0][0]<(finder[0][1]-5):
    #     for n in range(finder[0][1]):
    #         x_y_curveTest12.append(list((predScarp2X[n], z12new[n])))
    #         x_y_curveTest22.append(list((predScarp2X[n], z22[n])))
    # else:
    #     for n in range(finder[0][0]):
    #         x_y_curveTest12.append(list((predScarp2X[n], z12new[n])))
    #         x_y_curveTest22.append(list((predScarp2X[n], z22[n])))
    #
    # polygon_points2 = []  # creates a empty list where we will append the points to create the polygon
    #
    # for xyvalue in x_y_curveTest12:
    #     polygon_points2.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 1
    #
    # for xyvalue in x_y_curveTest22[::-1]:
    #     polygon_points2.append([xyvalue[0], xyvalue[1]])  # append all xy points for curve 2 in the reverse order (from last point to first point)
    #
    # for xyvalue in x_y_curveTest12[0:1]:
    #     polygon_points2.append([xyvalue[0], xyvalue[1]])  # append the first point in curve 1 again, to it "closes" the polygon
    #
    # polygon2 = Polygon(polygon_points2)
    # area2 = polygon2.area
    # predictedArea.append(area2)


# plt.figure()
# plt.plot(dist,surrogateProf)
# plt.plot(dist,g)
    # smoothed = moving_average(surrogateProf, 7)
    # extended = np.hstack((surrogateProf[0:3], smoothed))
    # extended2 = np.hstack((extended, surrogateProf[-4:-1]))
    #
    #
    # # surrogates2[hh,:] = extended2 #
    # # surrogates2[hh, :] = surrogateProf
    # newPCs = originalIPCA.transform(surrogateProf.reshape(1,-1))
    # newScarpProfile = surrogateProf
    # newPCsProfile = originalIPCA.inverse_transform(newPCs)



area = np.asarray(actualArea)


endTime = time.time()

print(endTime-startTime)

asdfg


import pandas as pd
import seaborn as sns
from sklearn.tree import DecisionTreeRegressor
data = np.vstack((designParams2.T,area))
data = np.vstack((data,a))
data = np.vstack((data,b))
data = np.vstack((data,distx))
data = np.vstack((data,initialShoreline))



df = pd.DataFrame(data.T,columns=['EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','MSL','Hs','Tp','NTR','Dur','Dir','M2','N2','K1','S2','tide','area','a','b','distx','shore'])

#X = df[['Hs','Tp','Dir','NTR','Dur','MSL','EOF1','tide']]
#y = df['area']

asdfg

df['WL'] = df['MSL'] + df['NTR']# + df['tide']
df['WP'] = 1000*9.81*9.81*df['Hs']*df['Hs']*df['Tp']/64/np.pi
df['shore'] = df['shore']-80
# df['Xbins'] = np.digitize(df.area,np.arange(0,100,10))

# d1 = df.assign(A_cut=pd.cut(df.area, np.arange(0,110,10), labels=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
#                B_cut=pd.cut(df.EOF1, np.arange(-20,50,10), labels=[1, 2, 3, 4, 5, 6]))
#
# d2 = d1.assign(cartesian=pd.Categorical(d1.filter(regex='_cut').apply(tuple,1)))
#
# tempMean = d2.groupby('cartesian').Hs.mean().unstack()
# tempHsMean = df.groupby([pd.cut(df.area, np.arange(0,110,10)),pd.cut(df.EOF1, np.arange(-20,50,10))]).Hs.mean().unstack()
from matplotlib import cm

areaBins = np.arange(0,102.5,10)
# eof1Bins = np.arange(-20,45.5,5)
# wpMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.EOF1, eof1Bins)]).WP.mean().unstack()
# wlMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.EOF1, eof1Bins)]).WL.mean().unstack()
# aMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.EOF1, eof1Bins)]).a.mean().unstack()
# bMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.EOF1, eof1Bins)]).b.mean().unstack()

shoreBins = np.arange(10,151,15)
wpMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.shore, shoreBins)]).WP.mean().unstack()
wlMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.shore, shoreBins)]).WL.mean().unstack()
aMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.shore, shoreBins)]).a.mean().unstack()
bMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.shore, shoreBins)]).b.mean().unstack()
distxMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.shore, shoreBins)]).distx.mean().unstack()

durMean = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.shore, shoreBins)]).Dur.median().unstack()-36
# durCat = pd.qcut(durMean,q=3)

upper_quantile = df.WP.quantile([.9])
wpMeanUpper = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.shore, shoreBins)]).WP.mean([df.WP > upper_quantile]).unstack()




xs = areaBins[0:-1]+3.5 #np.arange(3.5,100,7.5)
ys = shoreBins[0:-1]+5#np.arange(15,150,10)
xax,yax = np.meshgrid(xs,ys)
# xax,yax = np.meshgrid(np.arange(2.5,100,5),np.arange(-17.5,43,5))
# plt.scatter(np.arange(0,110,10),np.arange(-20,50,10),wpMean.T)
#plt.scatter(xax,yax,s=wlMean.T*500,c=wpMean.T,vmin=10000,vmax=400000,cmap=cm.OrRd)

fig = plt.figure()
ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
colorparam = np.zeros((len(xs)*len(ys),))
colormap = cm.plasma
normalize = mcolors.Normalize(vmin=.01, vmax=4.5)
counter = 0
pc1 = 0
pc2 = 0
pc3 = 0
for m in range(1,len(xs)):
    for n in range(1,len(ys)):
        xTemp = np.arange(0,1.4,0.025)
        # xTemp = np.arange(0,distxMean.iat[m,n]/15,0.025)
        zTemp = aMean.iat[m,n] * np.power(xTemp, bMean.iat[m,n])
        factor = distxMean.iat[m,n]/6
        xStretch = (xTemp*factor)+xs[m]-0.7
        # xStretch = (xTemp*4)+xs[m]-2
        yStretch = (zTemp*factor)+ys[n] - np.nanmin(zTemp*factor) - np.nanmax(np.abs(zTemp*factor))/2
        colorparam[counter] = wpMean.iat[m,n]/100000
        color = colormap(normalize(colorparam[counter]))
        # ax.plot(xStretch,yStretch,linewidth=(wlMean.iat[m,n]*6+1)/2,color=color)
        # ax.plot(xStretch,yStretch,linewidth=3,color=color)
        if durMean.iat[m,n] > 120:
            if pc1 == 0:
                p1 = ax.plot(xStretch,yStretch,linewidth=4.5,color=color,label='>5 days')
                pc1 = pc1+1
            else:
                p1 = ax.plot(xStretch,yStretch,linewidth=4.5,color=color,label='_nolegend_')

        elif durMean.iat[m,n] > 48 and durMean.iat[m,n] < 120:
            if pc2 == 0:
                p2 = ax.plot(xStretch,yStretch,linewidth=3,color=color,label='2-5 days')
                pc2 = pc2+1
            else:
                p2 = ax.plot(xStretch,yStretch,linewidth=3,color=color,label='_nolegend_')
        else:
            if pc3 == 0:
                p3 = ax.plot(xStretch,yStretch,linewidth=1.5,color=color,label='<2 days')
                pc3 = pc3+1
            else:
                p3 = ax.plot(xStretch,yStretch,linewidth=1.5,color=color,label='_nolegend_')
        counter = counter+1

# plt.tight_layout()
# ax.set_aspect('equal', 'box')
ax.set_xlabel('Eroded Volume ($m^3/m$)')
ax.set_ylabel('Initial Beach Width (m)')
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
s_map.set_array(colorparam)
fig.subplots_adjust(right=0.82)
cbar_ax = fig.add_axes([0.84, 0.15, 0.02, 0.7])
#fmt = lambda x, pos: '{:.1%}'.format(x)
cbar = fig.colorbar(s_map, cax=cbar_ax)#, format=FuncFormatter(fmt))
cbar.set_label('$\overline{WP}_{max}$ $(W/m x 10^6)$')
# cbar.ax.set_yticklabels(np.arange(cbar_min, cbar_max+cbar_step, cbar_step), fontsize=16, weight='bold')
#leg = ax.legend([p1,p2,p3],['>100 hrs','50-100 hrs','<50 hrs'])
leg = ax.legend()#[p1,p2,p3])
leg.get_lines()[0].set_linewidth(1.5)
leg.get_lines()[1].set_linewidth(3)
leg.get_lines()[2].set_linewidth(4.5)
ax.set_title('Sensitivity of Post-storm Sub-aerial Profile Shape')





#
#
# def q3(x):
#     return x.quantile(0.75)
# dfTest = pd.DataFrame(np.random.randn(100, 4), columns=['A', 'B', 'C', 'D'])#, index=mindex)
#
# #cMeanTest = dfTest.groupby([pd.qcut(dfTest.A, 2),pd.cut(dfTest.B, 2)]).C.mean().unstack()
#
# wpUpper = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.shore, shoreBins)]).WP.quantile([0.9]).unstack(level=1)
#
# wpCounts = df.groupby([pd.cut(df.area, areaBins),pd.cut(df.shore, shoreBins)]).WP.count().unstack(level=1)
#
# aUpper = np.nan*np.ones((np.shape(wpUpper)))
# for m in range(1,len(xs)):
#     for n in range(1,len(ys)):
#         if ~np.isnan(wpUpper.iat[m,n]):
#             aUpper[m, n] =
#             #df.groupby([pd.cut(df.area, areaBins), pd.cut(df.shore, shoreBins), pd.cut(df.WP, np.arange(wpUpper.iloc[m,n],10000000+wpUpper.iloc[m,n],9999999))]).a.mean().unstack(level=1).iat[m,n]
#         # dfSubset = dfTest.groupby([pd.cut(dfTest.A, 4),pd.cut(dfTest.B, 4),pd.cut(dfTest.C,np.arange(1.2,1000,900))])
#
# #
# for temp in dfSubset.groups:
#     #upperTest = group.C.quantile([.9])
#     #print(upperTest)
#     print(temp)
# print()
# print(dfSubset.groups.keys())
# labs = [f'bin_{x+1}' for x in range(3)]
#
# g = pd.cut(dfTest.groupby([pd.qcut(dfTest.A, 2),pd.cut(dfTest.B, 2)])['C'].mean(), bins=3, labels = labs)
#
# pd.cut(squares.groupby('square')['brightness'].transform('mean'), bins=4, labels=labs)
# dfs = dict(tuple(dfSubset))
#
#
# # dfTest['Xbins'] = np.digitize(dfTest.A,np.arange(-1.5,1.5,.25))
# d1Test = dfTest.assign(A_cut=pd.cut(dfTest.A, np.arange(-2.5,2.6,.5), labels=[1, 2, 3, 4, 5,6,7,8,9,10]),
#                B_cut=pd.cut(dfTest.B, np.arange(-2.5,2.6,.5), labels=[1, 2, 3, 4, 5,6,7,8,9,10]))
# d2 = d1Test.assign(cartesian=pd.Categorical(d1Test.filter(regex='_cut').apply(tuple,1)))
#
# d3 = d1Test.groupby('cartesian')
#
# # d1 = df.assign(A_cut=pd.cut(df.area, np.arange(0,110,10), labels=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
# #                B_cut=pd.cut(df.EOF1, np.arange(-20,50,10), labels=[1, 2, 3, 4, 5, 6]))
# #
# # d2 = d1.assign(cartesian=pd.Categorical(d1.filter(regex='_cut').apply(tuple,1)))
# #
# d2.sort_values(by='cartesian')
#
# for M in d2.cartesian:
#     print(M)
#     #print(d2.loc([d2['cartesian']==M]))


#
# cols = ['foo', 'bar', 'new']
# df['combined'] = df[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)


#
#
# import pandas as pd
# import seaborn as sns
# from sklearn.tree import DecisionTreeRegressor
# data = np.vstack((designParams2.T,area))
#
# df = pd.DataFrame(data.T,columns=['Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','M2','N2','K1','S2','tide','area'])
#
# X = df[['Hs','Tp','Dir','NTR','Dur','MSL','EOF1','tide']]
# y = df['area']
#
# # instantiating the model
# dtree = DecisionTreeRegressor()
#
# # fitting the model
# model = dtree.fit(X,y)



# # import export_graphviz
# from sklearn.tree import export_graphviz
#
# # export the decision tree to a tree.dot file
# # for visualizing the plot easily anywhere
# export_graphviz(dtree, out_file='tree.dot',
#                 feature_names=['Hs','Tp','Dir','NTR','Dur','MSL','EOF1','tide'])

# from sklearn import tree
# plt.figure(figsize=(40,20))
# _ = tree.plot_tree(dtree, feature_names = X.columns,
#              filled=True, fontsize=6, rounded = True)
# plt.show()



#
# from graphviz import Source
# from sklearn import tree
# #Source( tree.export_graphviz(dtree, out_file=None, feature_names=X.columns))
#
# graph = Source( tree.export_graphviz(dtree, out_file=None, feature_names=X.columns))
# graph.format = 'png'
# graph.render('dtree_render',view=True)
#

#
# from sklearn.externals.six import StringIO
# from sklearn.tree import export_graphviz
# import pydotplus
# from IPython.display import Image
#
#
# dot_data = StringIO()
#
# export_graphviz(
#     model,
#     out_file = dot_data,
#     filled=True, rounded=True, proportion=False,
#     special_characters=True,
#     feature_names=X.columns,
#     class_names=["cat", "dog"]
# )
#
# graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
#
# Image(graph.create_png())







numberOfSamples = 200000
minimumLevel = .00025
levelJump = 0.002

plt.figure()

# STILL WATER LEVEL
col = 10-1
ax131 = plt.subplot2grid((6,4),(0,0),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(a):np.max(a):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([a[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax131.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax131.set_xlim([-1.5,0.2])
ax131.set_ylabel('SWL')

ax132 = plt.subplot2grid((6,4),(0,1),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(b):np.max(b):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([b[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax132.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)

ax133 = plt.subplot2grid((6,4),(0,2),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(distx):np.max(distx):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([distx[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax133.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax133.set_xlim([5,40])

ax134 = plt.subplot2grid((6,4),(0,3),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(area):np.max(area):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([area[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax134.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax134.set_xlim([0,40])


# TIDE
col = 20-1
ax231 = plt.subplot2grid((6,4),(1,0),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(a):np.max(a):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([a[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax231.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax231.set_xlim([-1.5,0.2])
ax231.set_ylabel('tide')

ax232 = plt.subplot2grid((6,4),(1,1),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(b):np.max(b):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([b[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax232.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)

ax233 = plt.subplot2grid((6,4),(1,2),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(distx):np.max(distx):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([distx[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax233.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax233.set_xlim([5,40])

ax234 = plt.subplot2grid((6,4),(1,3),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(area):np.max(area):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([area[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax234.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax234.set_xlim([0,40])



# HS
col = 11-1
ax331 = plt.subplot2grid((6,4),(2,0),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(a):np.max(a):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([a[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax331.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax331.set_xlim([-1.5,0.2])
ax331.set_ylim([1.9,5])
ax331.set_ylabel('Hs')

ax332 = plt.subplot2grid((6,4),(2,1),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(b):np.max(b):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([b[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax332.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax332.set_ylim([1.9,5])

ax333 = plt.subplot2grid((6,4),(2,2),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(distx):np.max(distx):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([distx[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax333.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax333.set_xlim([5,40])
ax333.set_ylim([1.9,5])

ax334 = plt.subplot2grid((6,4),(2,3),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(area):np.max(area):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([area[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax334.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax334.set_xlim([0,40])
ax334.set_ylim([1.9,5])


# TP
col = 12-1
ax431 = plt.subplot2grid((6,4),(3,0),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(a):np.max(a):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([a[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax431.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax431.set_xlim([-1.5,0.2])
ax431.set_ylabel('Tp')

ax432 = plt.subplot2grid((6,4),(3,1),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(b):np.max(b):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([b[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax432.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)

ax433 = plt.subplot2grid((6,4),(3,2),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(distx):np.max(distx):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([distx[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax433.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax433.set_xlim([5,40])

ax434 = plt.subplot2grid((6,4),(3,3),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(area):np.max(area):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([area[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax434.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax434.set_xlim([0,40])

# Dur
col = 14-1
ax531 = plt.subplot2grid((6,4),(4,0),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(a):np.max(a):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([a[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax531.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax531.set_xlim([-1.5,0.2])
ax531.set_ylabel('Dur')

ax532 = plt.subplot2grid((6,4),(4,1),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(b):np.max(b):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([b[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax532.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)

ax533 = plt.subplot2grid((6,4),(4,2),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(distx):np.max(distx):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([distx[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax533.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax533.set_xlim([5,40])

ax534 = plt.subplot2grid((6,4),(4,3),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(area):np.max(area):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([area[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax534.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax534.set_xlim([0,40])

# EOF1
col = 1-1
ax631 = plt.subplot2grid((6,4),(5,0),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(a):np.max(a):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([a[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax631.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax631.set_xlim([-1.5,0.2])
ax631.set_ylabel('EOF1')
ax631.set_xlabel('a')

ax632 = plt.subplot2grid((6,4),(5,1),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(b):np.max(b):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([b[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax632.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax632.set_xlabel('b')

ax633 = plt.subplot2grid((6,4),(5,2),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(distx):np.max(distx):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([distx[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax633.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax633.set_xlim([5,40])
ax633.set_xlabel('distx')

ax634 = plt.subplot2grid((6,4),(5,3),rowspan=1,colspan=1)
xx, yy = np.mgrid[np.min(area):np.max(area):100j, np.min(dataS[0:numberOfSamples,col]):np.max(dataS[0:numberOfSamples,col]):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([area[0:numberOfSamples], dataS[0:numberOfSamples,col]])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
levels = np.arange(np.min(f),np.max(f),np.max(f)/100)
cset = ax634.contour(xx, yy, f, levels[0:50:2], cmap=cm.plasma)
ax634.set_xlim([0,40])
ax634.set_xlabel('volume lost')

# ax1 = plt.subplot2grid((5,5),(0,0),rowspan=1,colspan=1)
# dist_space = np.arange(1,8,0.1)
# kde1 = gaussian_kde(sampleHs)
# ax1.plot(dist_space, kde1(dist_space),linewidth=2,label='Historical')
# kde2 = gaussian_kde(hypo[0:1250,10])
# ax1.plot(dist_space, kde2(dist_space),'--',linewidth=2,label='MDA Selections')
# ax1.yaxis.set_ticks([])
# ax1.yaxis.set_ticklabels([])
# ax1.legend()
