
from getdatatestbed import getDataFRF
import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from glob import glob
import os
from scipy.interpolate import interp1d
import sandBarTool.morphLib as mL
import scipy.io
from matplotlib import gridspec
import pickle
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.metrics import silhouette_score
import scipy.cluster.hierarchy as shc



dbfile = open('sandbarsSouthernTransect_referencedMHHW_5lineAvg.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()
alllines = data['alllines']
xinterp = data['xinterp']
time = data['time']


dbfile2 = open('ceofsSouthernTransectLatest.pickle', 'rb')
data2 = pickle.load(dbfile2)
dbfile2.close()

timeS = data2['time']
alllinesS = data2['alllines']
xinterpS = data2['xinterp']
SS = data2['S']
RS = data2['Rt']
thetaS = data2['thetaRadians']
theta2S = data2['thetaDegrees']
phitS = data2['phiRadian']
phit2S = data2['phiDegrees']
totalVS = data2['totalV']
lambaS = data2['lamda']
percentVS = data2['percentV']


def P2R(radii, angles):
    return radii * np.exp(1j*angles)

def R2P(x):
    return np.abs(x), np.angle(x)

def pol2cart(theta, rho):
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return(x, y)

x1, y1 = pol2cart(phitS[:,0],RS[:,0])
# x1_2, y1_2 = P2R(phitS[:,1],RS[:,1])
x2, y2 = pol2cart(phitS[:,1],RS[:,1])
x3, y3 = pol2cart(phitS[:,2],RS[:,2])
x4, y4 = pol2cart(phitS[:,3],RS[:,3])

modes = 2
data = np.empty((len(x1),int(2*modes)))
data[:,0] = np.real(x1)*(percentVS[0]/percentVS[0])
data[:,1] = np.real(y1)*(percentVS[0]/percentVS[0])
data[:,2] = np.real(x2)*(percentVS[1]/percentVS[0])
data[:,3] = np.real(y2)*(percentVS[1]/percentVS[0])
# data[:,4] = np.real(x3)*(percentVS[2]/percentVS[0])
# data[:,5] = np.real(y3)*(percentVS[2]/percentVS[0])
# data[:,6] = np.real(x4)*(percentVS[3]/percentVS[0])
# data[:,7] = np.real(y4)*(percentVS[3]/percentVS[0])

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler

###### Regular PCA
# X = alllines
# # Scaling the data so that all the features become comparable
# scaler = StandardScaler()
# X_scaled = scaler.fit_transform(X)
#
# # Normalizing the data so that the data approximately
# # follows a Gaussian distribution
# X_normalized = normalize(X_scaled)
#
# # Converting the numpy array into a pandas DataFrame
# X_normalized = pd.DataFrame(X_normalized)
#
#
# pca = PCA(n_components = 2)
# X_principal = pca.fit_transform(X_normalized)
# X_principal = pd.DataFrame(X_principal)
# X_principal.columns = ['P1', 'P2']
#
#
# plt.figure(figsize =(8, 8))
# plt.title('Visualising the data')
# Dendrogram = shc.dendrogram((shc.linkage(X_principal, method ='ward')))
#
# ac9 = AgglomerativeClustering(n_clusters=9)
#
# # Visualizing the clustering
# plt.figure(figsize=(6, 6))
# plt.scatter(X_principal['P1'], X_principal['P2'],
#             c=ac9.fit_predict(X_principal), cmap='rainbow')
# plt.show()


###### DBSCAN
#X = StandardScaler().fit_transform(X)
#
# # #############################################################################
# # Compute DBSCAN
# db = DBSCAN(eps=0.5, min_samples=5).fit(X)
# core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
# core_samples_mask[db.core_sample_indices_] = True
# labels = db.labels_
#
# # Number of clusters in labels, ignoring noise if present.
# n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
# n_noise_ = list(labels).count(-1)
#
# print('Estimated number of clusters: %d' % n_clusters_)
# print('Estimated number of noise points: %d' % n_noise_)
# # print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
# # print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
# # print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
# # print("Adjusted Rand Index: %0.3f"
# #       % metrics.adjusted_rand_score(labels_true, labels))
# # print("Adjusted Mutual Information: %0.3f"
# #       % metrics.adjusted_mutual_info_score(labels_true, labels))
# print("Silhouette Coefficient: %0.3f"
#       % metrics.silhouette_score(X, labels))
#
# # #############################################################################
# # Plot result
# import matplotlib.pyplot as plt
#
# # Black removed and is used for noise instead.
# unique_labels = set(labels)
# colors = [plt.cm.Spectral(each)
#           for each in np.linspace(0, 1, len(unique_labels))]
# for k, col in zip(unique_labels, colors):
#     if k == -1:
#         # Black used for noise.
#         col = [0, 0, 0, 1]
#
#     class_member_mask = (labels == k)
#
#     xy = X[class_member_mask & core_samples_mask]
#     plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
#              markeredgecolor='k', markersize=14)
#
#     xy = X[class_member_mask & ~core_samples_mask]
#     plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
#              markeredgecolor='k', markersize=6)
#
# plt.title('Estimated number of clusters: %d' % n_clusters_)
# plt.show()
#





def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )




def transition_matrix(transitions):
    n = 1+ max(transitions) #number of states
    M = [[0]*n for _ in range(n)]
    M2 = [[0]*n for _ in range(n)]
    for (i,j) in zip(transitions,transitions[1:]):
        M[i][j] += 1
        M2[i][j] += 1
    #now convert to probabilities:
    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    return M, M2


def conditional_transition_matrix(transitions,num_of_states):
    n = num_of_states #1+ max(transitions) #number of states
    M = [[0]*n for _ in range(n)]
    M2 = [[0]*n for _ in range(n)]
    for (i,j) in zip(transitions,transitions[1:]):
        M[i][j] += 1
        M2[i][j] += 1
    #now convert to probabilities:
    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    return M, M2



#



# import datetime
# def datetime2matlabdn(dt):
#    ord = dt.toordinal()
#    mdn = dt + datetime.timedelta(days = 366)
#    frac = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
#    return mdn.toordinal() + frac

# matlabTime = np.zeros((len(time),))
# for qq in range(len(time)):
#     matlabTime[qq] = datetime2matlabdn(time[qq])
#



from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering


def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)




X = data

# setting distance_threshold=0 ensures we compute the full tree.
model = AgglomerativeClustering(distance_threshold=0, n_clusters=None)
model1 = AgglomerativeClustering(n_clusters=1)
model2 = AgglomerativeClustering(n_clusters=2)
model3 = AgglomerativeClustering(n_clusters=3)
model4 = AgglomerativeClustering(n_clusters=4)
model5 = AgglomerativeClustering(n_clusters=5)
model6 = AgglomerativeClustering(n_clusters=6)
model7 = AgglomerativeClustering(n_clusters=7)
model8 = AgglomerativeClustering(n_clusters=8)
model9 = AgglomerativeClustering(n_clusters=9)
model10 = AgglomerativeClustering(n_clusters=10)
model11 = AgglomerativeClustering(n_clusters=11)
model12 = AgglomerativeClustering(n_clusters=12)
model13 = AgglomerativeClustering(n_clusters=13)
model14 = AgglomerativeClustering(n_clusters=14)
model15 = AgglomerativeClustering(n_clusters=15)
model16 = AgglomerativeClustering(n_clusters=16)
model17 = AgglomerativeClustering(n_clusters=17)
model18 = AgglomerativeClustering(n_clusters=18)
model19 = AgglomerativeClustering(n_clusters=19)
model20 = AgglomerativeClustering(n_clusters=20)

k = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

from sklearn.metrics import silhouette_score
# Appending the silhouette scores of the different models to the list
silhouette_scores = []
silhouette_scores.append(silhouette_score(data, model2.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model3.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model4.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model5.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model6.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model7.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model8.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model9.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model10.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model11.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model12.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model13.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model14.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model15.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model16.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model17.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model18.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model19.fit_predict(data)))
silhouette_scores.append(silhouette_score(data, model20.fit_predict(data)))



plt.figure()
# Plotting a bar graph to compare the results
plt.plot(k, silhouette_scores)
plt.xlabel('Number of clusters', fontsize=20)
plt.ylabel('S(i)', fontsize=20)
plt.show()


from sklearn.metrics import davies_bouldin_score
# Appending the silhouette scores of the different models to the list
davies_bouldin = []
davies_bouldin.append(davies_bouldin_score(data, model2.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model3.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model4.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model5.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model6.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model7.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model8.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model9.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model10.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model11.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model12.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model13.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model14.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model15.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model16.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model17.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model18.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model19.fit_predict(data)))
davies_bouldin.append(davies_bouldin_score(data, model20.fit_predict(data)))

plt.figure()
# Plotting a bar graph to compare the results
plt.plot(k, davies_bouldin)
plt.xlabel('Number of clusters', fontsize=20)
plt.ylabel('davies_bouldin', fontsize=20)
plt.show()


from sklearn.metrics import calinski_harabasz_score
# Appending the silhouette scores of the different models to the list
calinski_harabasz = []
calinski_harabasz.append(calinski_harabasz_score(data, model2.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model3.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model4.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model5.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model6.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model7.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model8.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model9.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model10.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model11.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model12.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model13.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model14.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model15.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model16.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model17.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model18.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model19.fit_predict(data)))
calinski_harabasz.append(calinski_harabasz_score(data, model20.fit_predict(data)))

plt.figure()
# Plotting a bar graph to compare the results
plt.plot(k, calinski_harabasz)
plt.xlabel('Number of clusters', fontsize=20)
plt.ylabel('calinski_harabasz', fontsize=20)
plt.show()

import scipy.cluster.hierarchy as shc
Dendrogram = shc.dendrogram((shc.linkage(data, method ='ward')))
plt.figure(figsize =(6, 6))
plt.scatter(data[:,0], data[:,1],
            c = model10.fit_predict(data), cmap ='rainbow')
plt.show()


plt.figure()
model = AgglomerativeClustering(distance_threshold=0, n_clusters=None)
model = model.fit(data)
plt.title('Hierarchical Clustering Dendrogram')
# plot the top three levels of the dendrogram
plot_dendrogram(model, truncate_mode='level', p=3)
#plt.plot(np.arange(0,200),4.5*np.ones(200,))
plt.xlabel("Number of points in node")
plt.show()



yc = model11.fit_predict(data)
bmu2 = model2.fit_predict(data)
bmu3 = model3.fit_predict(data)
bmu4 = model4.fit_predict(data)
bmu5 = model5.fit_predict(data)
bmu6 = model6.fit_predict(data)
bmu7 = model7.fit_predict(data)
bmu8 = model8.fit_predict(data)

import itertools
ii = itertools.count(data.shape[0])
clusters = [{'node_id': next(ii), 'left': x[0], 'right':x[1]} for x in model.children_]

working = dict(enumerate(model.children_, model.n_leaves_))

import copy
n_points = data.shape[0]
members = {i:[i] for i in range(n_points)}
for cluster in clusters:
    node_id = cluster["node_id"]
    members[node_id] = copy.deepcopy(members[cluster["left"]])
    members[node_id].extend(copy.deepcopy(members[cluster["right"]]))

on_split = {c["node_id"]: [c["left"], c["right"]] for c in clusters}
up_merge = {c["left"]: {"into": c["node_id"], "with": c["right"]} for c in clusters}
up_merge.update({c["right"]: {"into": c["node_id"], "with": c["left"]} for c in clusters})

links = scipy.cluster.hierarchy.linkage(data)

X = data
plt.figure()
plt.scatter(X[yc==0, 0], X[yc==0, 1], s=20, label ='Cluster 1')
plt.scatter(X[yc==1, 0], X[yc==1, 1], s=20, label ='Cluster 2')
plt.scatter(X[yc==2, 0], X[yc==2, 1], s=20, label ='Cluster 3')
plt.scatter(X[yc==3, 0], X[yc==3, 1], s=20, label ='Cluster 4')
plt.scatter(X[yc==4, 0], X[yc==4, 1], s=20, label ='Cluster 5')
plt.scatter(X[yc==5, 0], X[yc==5, 1], s=20, label ='Cluster 6')
plt.scatter(X[yc==6, 0], X[yc==6, 1], s=20, label ='Cluster 7')
plt.scatter(X[yc==7, 0], X[yc==7, 1], s=20, label ='Cluster 8')
plt.scatter(X[yc==8, 0], X[yc==8, 1], s=20, label ='Cluster 9')
plt.scatter(X[yc==9, 0], X[yc==9, 1], s=20, label ='Cluster 10')
plt.scatter(X[yc==10, 0], X[yc==10, 1], s=20, label ='Cluster 11')
plt.scatter(X[yc==11, 0], X[yc==11, 1], s=20, label ='Cluster 12')
plt.scatter(X[yc==12, 0], X[yc==12, 1], s=20, label ='Cluster 13')

plt.title('Clusters of Customers (Hierarchical Clustering Model)')
# plt.xlabel('Annual Income(k$)')
# plt.ylabel('Spending Score(1-100')
plt.show()

numClusters = 11
plt.figure()
centroids = np.empty((numClusters,modes*2))
indices = list()
groupSize = np.empty((numClusters,))
for ii in range(numClusters):
    index = np.where(yc == ii)
    indices.append(index)
    groupSize[ii] = len(index[0])

    tempData = data[index[0],:]

    centroids[ii,:] = np.nanmean(tempData,axis=0)
    avgProf = np.nanmean(alllines[index[0],:],axis=0)
    plt.plot(xinterpS, avgProf)
plt.show()


profiles = np.zeros((numClusters,len(xinterp)))
magPhaseCentroid = np.zeros((np.shape(centroids)))
orderAngle = np.empty((numClusters,))
orderMag = np.empty((numClusters,))
orderAngle2 = np.empty((numClusters,))
orderMag2 = np.empty((numClusters,))

for i in range(numClusters):

    cpcRt = np.zeros((modes,))
    cpcPhi = np.zeros((modes,))

    cpcRt[0], cpcPhi[0] = R2P(complex(centroids[i, 0], centroids[i, 1]))
    cpcRt[1], cpcPhi[1] = R2P(complex(centroids[i, 2]/(percentVS[1]/percentVS[0]), centroids[i, 3]/(percentVS[1]/percentVS[0])))
    # cpcRt[2], cpcPhi[2] = R2P(complex(centroids[i, 4]/(percentVS[2]/percentVS[0]), centroids[i, 5]/(percentVS[2]/percentVS[0])))
    orderMag[i] = cpcRt[0]
    orderAngle[i] = cpcPhi[0]
    orderMag2[i] = cpcRt[1]
    orderAngle2[i] = cpcPhi[1]

    magPhaseCentroid[i,0:4] = [cpcRt[0], cpcPhi[0]*180/np.pi, cpcRt[1], cpcPhi[1]*180/np.pi]#, cpcRt[2],cpcPhi[2]*180/np.pi]
                                # cpcRt[3], cpcPhi[3], cpcRt[4], cpcPhi[4], cpcRt[5], cpcPhi[5],
                                # cpcRt[6], cpcPhi[6], cpcRt[7], cpcPhi[7]]

    profile = 0 * np.ones(len(xinterp), )
    for mode in range(modes):
        # profile = profile + cpcRt[mode] * np.sin(cpcPhi[mode]) * S[:, mode] * np.sin(theta[:, mode]) + cpcRt[mode]* np.cos(cpcPhi[mode]) * S[:, mode] * np.cos(theta[:, mode])
        profile = profile + cpcRt[mode] * SS[:, mode] * np.cos(cpcPhi[mode] - thetaS[:, mode])


    profiles[i,:] = profile+np.mean(alllines,axis=0)


negAngles = np.where((orderAngle<0))
orderAngle = orderAngle*180/np.pi
orderAngle[negAngles] = orderAngle[negAngles]+360
val = np.sort(orderAngle)
order = np.argsort(orderAngle)

negAngles2 = np.where((orderAngle2<0))
orderAngle2 = orderAngle2*180/np.pi
orderAngle2[negAngles2] = orderAngle2[negAngles2]+360

#madeUpOrder = [6,]
# progression = centroids[:,1]
# negAngles2 = np.where(progression<0)
# progression[negAngles2[0]] = orderAngle[negAngles2[0]]+2*np.pi
# val = np.sort(progression)
# order = np.argsort(progression)

plt.style.use('default')
import matplotlib.cm as cm
#colors = cm.rainbow(np.linspace(0, 1, numClusters))
#colors2 = cm.gray(np.linspace(0, 1, numClusters))
colors = cm.hsv(np.linspace(0, 1, (numClusters+1)))

plt.figure()
for ii in range(numClusters):
    finder = order[ii]
    plt.plot(centroids[finder,0],centroids[finder,1],'o',color=colors[ii,:],label=orderAngle[finder])
    #plt.plot(orderMag[finder[0][0],0],centroids[finder[0][0],1],'o',color=colors[ii,:],label=ii)

plt.legend()
plt.show()

plt.figure()
for ii in range(numClusters):
    finder = order[ii]

    index = np.where(yc == finder)

    tempData = data[index[0],:]
    plt.plot(tempData[:,0],tempData[:,1],'o',color=colors[ii,:],label=orderAngle[finder])
plt.legend()
plt.show()

import cmocean
import matplotlib.colors
cmap = cm.rainbow
norm = matplotlib.colors.Normalize(vmin=0,vmax=360)
plt.figure()
ax1 = plt.subplot2grid((6,6), (0,0), rowspan=3, colspan=6)
for ii in range(numClusters):
    finder = order[ii]
    # ax1.plot(xinterpS, profiles[finder,:], color=cmap(norm(orderAngle[finder])), label = r'{} @ {} and {} @ {}'.format(np.round(100*orderMag[finder])/100, np.round(orderAngle[finder]), np.round(100*orderMag2[finder])/100, np.round(orderAngle2[finder])))
    # ax1.plot(xinterpS, profiles[finder,:],10,orderAngle[finder]*np.ones((np.shape(xinterpS))),vmin=0,vmax=360,cmap=cmocean.cm.phase,label = r'{} @ {} and {} @ {}'.format(np.round(100*orderMag[finder])/100,np.round(orderAngle[finder]),np.round(100*orderMag2[finder])/100,np.round(orderAngle2[finder])))
    ax1.plot(xinterpS, profiles[finder,:],color=colors[ii,:],label = r'{} @ {} and {} @ {}'.format(np.round(100*orderMag[finder])/100,np.round(orderAngle[finder]),np.round(100*orderMag2[finder])/100,np.round(orderAngle2[finder])))

ax1.legend()
ax1.set_title('Centroids')

ax2 = plt.subplot2grid((6,6),(3,0),rowspan=3,colspan=3)
ids = np.arange(6,11,1)
for ii in ids:
    finder = order[ii]
    ax2.plot(xinterpS, profiles[finder,:],color=colors[ii,:],label = r'{} @ {} and {} @ {}'.format(np.round(100*orderMag[finder])/100,np.round(orderAngle[finder]),np.round(100*orderMag2[finder])/100,np.round(orderAngle2[finder])))
ax2.legend()
ax2.set_title('Emergence of New Bar')

ax3 = plt.subplot2grid((6,6),(3,3),rowspan=3,colspan=3)
ids = np.arange(0,6,1)
for ii in ids:
    finder = order[ii]
    ax3.plot(xinterpS, profiles[finder,:],color=colors[ii,:],label = r'{} @ {} and {} @ {}'.format(np.round(100*orderMag[finder])/100,np.round(orderAngle[finder]),np.round(100*orderMag2[finder])/100,np.round(orderAngle2[finder])))
ax3.legend()
ax3.set_title('Rapid Migration')

plt.show()

colorsOrig = copy.deepcopy(colors)


def recta(x1, y1, x2, y2):
    a = (y1 - y2) / (x1 - x2)
    b = y1 - a * x1
    return (a, b)

def curva_b(xa, ya, xb, yb, xc, yc):
    (x1, y1, x2, y2) = (xa, ya, xb, yb)
    (a1, b1) = recta(xa, ya, xb, yb)
    (a2, b2) = recta(xb, yb, xc, yc)
    puntos = []

    for i in range(0, 1000):
        if x1 == x2:
            continue
        else:
            (a, b) = recta(x1, y1, x2, y2)
        x = i*(x2 - x1)/1000 + x1
        y = a*x + b
        puntos.append((x,y))
        x1 += (xb - xa)/1000
        y1 = a1*x1 + b1
        x2 += (xc - xb)/1000
        y2 = a2*x2 + b2
    return puntos

def hanging_line(point1, point2):
    import numpy as np

    a = (point2[1] - point1[1])/(np.cosh(point2[0]) - np.cosh(point1[0]))
    b = point1[1] - a*np.cosh(point1[0])
    x = np.linspace(point1[0], point2[0], 100)
    y = a*np.cosh(x) + b

    return (x,y)


sorted_bmus = np.tile(0,(len(yc),), )
for i in range(numClusters):
    posc = np.where(yc == order[i])
    sorted_bmus[posc] = int(i)
    #sorted_time[posc] = time[posc]


#time = data['time']
dt = timeS[1:]-timeS[0:-1]
days = [t.days for t in dt]#[item for sublist in m for item in sublist]
days2 = np.array(days)

bins = dict()
date = dict()
nextbin = dict()
nextdate = dict()
prevbin = dict()
prevdate = dict()
daysBetween = dict()
for xx in range(numClusters):
    binind = np.where((sorted_bmus==xx))
    bins[xx] = binind
    date[xx] = time[binind]
    nextbinind = np.where((sorted_bmus[0:-1]==xx))
    K = 1
    res = [x + K for x in nextbinind]
    nextbin[xx] = sorted_bmus[res]
    nextdate[xx] = time[res]
    prevbinind = np.where((sorted_bmus==xx))
    res2 = [x - K for x in prevbinind]
    prevbin[xx] = sorted_bmus[res2]
    prevdate[xx] = time[res2]
    daysBetween[xx] = days2[res2]


m,m2 = transition_matrix(sorted_bmus)
for row in m: print(' '.join('{0:.2f}'.format(x) for x in row))


flat_list = [item for sublist in m for item in sublist]
flatarray = np.asarray(flat_list)
flatarray.resize(numClusters, numClusters)


tooLong = np.where(days2>65)
transits = []
for hh in range(len(tooLong[0])-1):
    if hh == 0:
        ind = np.where((time < timeS[tooLong[0][hh+1]]))
        tran1, tran2 = conditional_transition_matrix(sorted_bmus[ind], numClusters)
    elif hh == len(tooLong):
        ind = np.where((time > timeS[tooLong[0][hh]]))
        tran1, tran2 = conditional_transition_matrix(sorted_bmus[ind], numClusters)
    else:
        ind = np.where((time > timeS[tooLong[0][hh]]) & (time < timeS[tooLong[0][hh+1]]))
        tran1, tran2 = conditional_transition_matrix(sorted_bmus[ind], numClusters)

    if len(ind[0] > 1):
        flat_listCond = [item for sublist in tran2 for item in sublist]
        flatarrayCond = np.asarray(flat_listCond)
        flatarrayCond.resize(numClusters, numClusters)


        transits.append(flatarrayCond)

allTransits = np.sum(transits,axis=0)

def convertToProbabilities(totalTransitions):
    mm, nn = np.shape(totalTransitions)
    output = np.zeros((mm,nn))
    for qq in range(mm):
        s = sum(totalTransitions[qq])
        if s > 0:
            output[qq,0:len(totalTransitions[qq])] = [f / s for f in totalTransitions[qq]]
    return output


probs = convertToProbabilities(allTransits)

for row in probs: print(' '.join('{0:.2f}'.format(x) for x in row))

flat_list = [item for sublist in probs for item in sublist]
flatarray = np.asarray(flat_list)
flatarray.resize(numClusters, numClusters)







tooHigh = np.where((orderAngle > 180))
orderAngle[tooHigh[0]] = orderAngle[tooHigh[0]] - 360


tooHigh2 = np.where((orderAngle2 > 180))
orderAngle2[tooHigh2[0]] = orderAngle2[tooHigh2[0]] - 360


colors = np.vstack((colorsOrig[6:-1,:],colorsOrig[0:6,:]))


colors = cm.gist_rainbow(np.linspace(0, 1, (numClusters)))

#colors = cm.gist_ncar(np.linspace(0, 1, (numClusters)))

#colors = cm.rainbow(np.linspace(0, 1, (numClusters+1)))


alfa = 0.29

fig = plt.figure(figsize=(14,14))
ax10 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
for ii in range(11): #range(numClusters):
    finder = ii#order[ii]
    index = np.where(sorted_bmus == finder)

    if ii == 0:

        mode2phis = phit2S[index[0], 1]
        index2 = np.where((mode2phis< 0))
        mode2phis[index2[0]] = mode2phis[index2[0]]+360
        ax10.scatter(phit2S[index[0], 0], mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')
    elif ii == 1:
        ax10.scatter(phit2S[index[0], 0], phit2S[index[0], 1], 15, color=colors[ii, :],alpha=alfa,edgecolor='none')
        ax10.scatter(phit2S[index[0], 0], phit2S[index[0],1]+360, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')

    elif ii == 2:
        mode1phis = phit2S[index[0], 0]
        index1 = np.where((mode1phis< -50))
        mode1phis[index1[0]] = mode1phis[index1[0]]+360
        mode2phis = phit2S[index[0], 1]
        index2 = np.where((mode2phis< 0))
        mode2phis[index2[0]] = mode2phis[index2[0]]+360
        ax10.scatter(mode1phis, mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')
        # ax10.scatter(mode1phis, mode2phis+360, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')

    elif ii == 3:
        mode1phis = phit2S[index[0], 0]
        index1 = np.where((mode1phis< 0))
        mode1phis[index1[0]] = mode1phis[index1[0]]+360
        mode2phis = phit2S[index[0], 1]
        ax10.scatter(mode1phis, mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')
        ax10.scatter(mode1phis-360, mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')

    elif ii == 4:
        mode1phis = phit2S[index[0], 0]
        mode2phis = phit2S[index[0], 1]
        index2 = np.where((mode2phis< 0))
        mode2phis[index2[0]] = mode2phis[index2[0]]+360
        ax10.scatter(mode1phis, mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')
    elif ii == 5:
        mode1phis = phit2S[index[0], 0]
        index1 = np.where((mode1phis< -50))
        mode1phis[index1[0]] = mode1phis[index1[0]]+360
        mode2phis = phit2S[index[0], 1]
        index2 = np.where((mode2phis< 0))
        mode2phis[index2[0]] = mode2phis[index2[0]]+360
        ax10.scatter(mode1phis, mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')
    elif ii == 6:
        mode1phis = phit2S[index[0], 0]
        index1 = np.where((mode1phis>-50))
        mode1phis[index1[0]] = mode1phis[index1[0]]-360
        mode2phis = phit2S[index[0], 1]
        index2 = np.where((mode2phis> 0))
        mode2phis[index2[0]] = mode2phis[index2[0]]-360
        ax10.scatter(mode1phis, mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')
        ax10.scatter(mode1phis+360, mode2phis+360, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')

    elif ii == 7:
        mode1phis = phit2S[index[0], 0]
        mode2phis = phit2S[index[0], 1]
        index2 = np.where((mode2phis> 75))
        mode2phis[index2[0]] = mode2phis[index2[0]]-360
        ax10.scatter(mode1phis, mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')
        ax10.scatter(mode1phis+360, mode2phis+360, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')

    elif ii == 8:
        mode1phis = phit2S[index[0], 0]
        index1 = np.where((mode1phis<-130))
        mode1phis[index1[0]] = mode1phis[index1[0]]+360
        mode2phis = phit2S[index[0], 1]
        ax10.scatter(mode1phis, mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')

    elif ii == 10:
        mode1phis = phit2S[index[0], 0]
        mode2phis = phit2S[index[0], 1]
        index2 = np.where((mode2phis< 0))
        mode2phis[index2[0]] = mode2phis[index2[0]]+360
        ax10.scatter(mode1phis, mode2phis, 15, color=colors[ii, :],alpha=alfa,edgecolor='none')
    else:
        ax10.scatter(phit2S[index[0], 0], phit2S[index[0], 1],15,color=colors[ii,:],alpha=alfa,edgecolor='none')





#
# probMultiplier = 200
# for ii in range(11): #range(numClusters):
#     finder = order[ii]
#     #index = np.where(yc == finder)
#     if ii == 5:
#         ax10.scatter(orderAngle[finder],orderAngle2[finder]+360,probMultiplier*flatarray[ii,ii],color=colors[ii,:],edgecolor='k')
#     elif ii == 6:
#         ax10.scatter(orderAngle[finder],orderAngle2[finder],probMultiplier*flatarray[ii,ii],color=colors[ii,:],edgecolor='k')
#         ax10.scatter(orderAngle[finder]+360,orderAngle2[finder]+360,probMultiplier*flatarray[ii,ii],color=colors[ii,:],edgecolor='k')
#     elif ii == 7:
#         ax10.scatter(orderAngle[finder],orderAngle2[finder],probMultiplier*flatarray[ii,ii],color=colors[ii,:],edgecolor='k')
#         ax10.scatter(orderAngle[finder]+360,orderAngle2[finder]+360,probMultiplier*flatarray[ii,ii],color=colors[ii,:],edgecolor='k')
#
#     else:
#         ax10.scatter(orderAngle[finder],orderAngle2[finder],probMultiplier*flatarray[ii,ii],color=colors[ii,:],edgecolor='k')
#

smallSize = 1
midSize = 3
largeSize = 6

arrowSize1=15
arrowSize2=25
arrowSize3=35
lineSizeMultiplier = 15
for ii in range(11):
    if ii == 0:
        goto = np.array((1,2))
        # goto = np.array((1,2,8))
        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder],orderAngle2[finder]]
            point2 = [orderAngle[finder2],orderAngle2[finder2]]
            if qq == 1:
                point2 = [orderAngle[finder2], orderAngle2[finder2]+360]

                point3x = orderAngle[finder]+10#2*np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),3)
                point3y = orderAngle2[finder]+30#-np.divide(np.abs(orderAngle2[finder]+orderAngle2[finder2]),2)
            elif qq == 4:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),3)
                point3y = orderAngle2[finder]+5*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),2)

            elif qq == 8:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]+orderAngle[finder2]),3)
                point3y = orderAngle2[finder]-2*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),2)
            else:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),2)
                point3y = orderAngle2[finder]+np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),2)
            print('working on {}'.format(qq))
        #x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0],point1[1], point3x, point3y, point2[0],point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    elif ii == 1:
        # goto = np.array((0,2,3,4,5))
        goto = np.array((2,3))
        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder],orderAngle2[finder]]
            point2 = [orderAngle[finder2],orderAngle2[finder2]]
            if qq == 1:
                point3x = orderAngle[finder]+2*np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),3)
                point3y = orderAngle2[finder]-np.divide(np.abs(orderAngle2[finder]+orderAngle2[finder2]),2)
            elif qq == 4:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),3)
                point3y = orderAngle2[finder]+2*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),3)
            elif qq == 5:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]+orderAngle[finder2]),3)
                point3y = orderAngle2[finder]-2*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),2)
                point2 = [orderAngle[finder2],orderAngle2[finder2]+360]
            elif qq == 2:
                # point2 = [orderAngle[finder2], orderAngle2[finder2]+360]

                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),2)
                point3y = orderAngle2[finder]+2*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),3)
            else:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),2)
                point3y = orderAngle2[finder]+1*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),3)
            print('working on {}'.format(qq))
        #x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0],point1[1], point3x, point3y, point2[0],point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    elif ii == 2:
        # goto = np.array((0,1,4,5,7))
        goto = np.array((0,4,5))
        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder],orderAngle2[finder]]
            point2 = [orderAngle[finder2],orderAngle2[finder2]]
            if qq == 1:
                point3x = orderAngle[finder]-1*np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),10)
                point3y = orderAngle2[finder]-np.divide(np.abs(orderAngle2[finder]+orderAngle2[finder2]),2)
            elif qq == 4:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),3)
                point3y = orderAngle2[finder]+2*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),3)
            elif qq == 5:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),3)
                point3y = orderAngle2[finder]+2*np.divide(np.abs(orderAngle2[finder]-(orderAngle2[finder2]+360)),2)
                point2 = [orderAngle[finder2],orderAngle2[finder2]+360]
            elif qq == 0:
                point3x = orderAngle[finder]-np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),2)
                point3y = orderAngle2[finder]+1*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),4)
            elif qq == 7:
                point2 = [orderAngle[finder2]+360, orderAngle2[finder2]+360]
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-point2[0]),2)
                point3y = orderAngle2[finder]+1*np.divide(np.abs(orderAngle2[finder]-point2[1]),3)
            else:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),2)
                point3y = orderAngle2[finder]+1*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),3)
            print('working on {}'.format(qq))
        #x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0],point1[1], point3x, point3y, point2[0],point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    elif ii == 3:
        # goto = np.array((1,2,4,5,8))
        goto = np.array((2,8))
        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder],orderAngle2[finder]]
            point2 = [orderAngle[finder2],orderAngle2[finder2]]
            if qq == 1:
                point3x = orderAngle[finder]-1*np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),10)
                point3y = orderAngle2[finder]-np.divide(np.abs(orderAngle2[finder]+orderAngle2[finder2]),2)
            elif qq == 4:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),3)
                point3y = orderAngle2[finder]+2*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),3)
            elif qq == 8:
                point1 = [orderAngle[finder]-360, orderAngle2[finder]]
                point3x = -100#orderAngle[finder]-np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),3)
                point3y = 80#orderAngle2[finder]+2*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),3)
            elif qq == 5:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),3)
                point3y = orderAngle2[finder]+2*np.divide(np.abs(orderAngle2[finder]-(orderAngle2[finder2]+360)),2)
                point2 = [orderAngle[finder2],orderAngle2[finder2]+360]
            elif qq == 0:
                point3x = orderAngle[finder]-np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),2)
                point3y = orderAngle2[finder]+1*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),4)
            else:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),2)
                point3y = orderAngle2[finder]+1*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),3)
            print('working on {}'.format(qq))
        #x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0],point1[1], point3x, point3y, point2[0],point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    elif ii == 4:
        # goto = np.array((2,5,6))
        goto = np.array((5,6))

        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder],orderAngle2[finder]]
            point2 = [orderAngle[finder2],orderAngle2[finder2]]
            if qq == 2:
                point3x = orderAngle[finder]-1*np.divide(np.abs(orderAngle[finder]-orderAngle[finder2]),10)
                point3y = orderAngle2[finder]-4*np.divide(np.abs(orderAngle2[finder]-orderAngle2[finder2]),3)
            elif qq == 5:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-(orderAngle[finder2])),3)
                point3y = orderAngle2[finder]+2*np.divide(np.abs(orderAngle2[finder]-(orderAngle2[finder2]+360)),2)
                point2 = [orderAngle[finder2],orderAngle2[finder2]+360]

            elif qq == 6:
                point3x = orderAngle[finder]+np.divide(np.abs(orderAngle[finder]-(orderAngle[finder2]+360)),3)
                point3y = orderAngle2[finder]+2*np.divide(np.abs(orderAngle2[finder]-(orderAngle2[finder2]+360)),3)
                point2 = [orderAngle[finder2]+360,orderAngle2[finder2]+360]

            print('working on {}'.format(qq))
        #x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0],point1[1], point3x, point3y, point2[0],point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    elif ii == 5:
        # goto = np.array((2,3,6,7))
        goto = np.array((6,7))

        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder],orderAngle2[finder]+360]
            point2 = [orderAngle[finder2],orderAngle2[finder2]]
            if qq == 2:
                point1 = [orderAngle[finder], orderAngle2[finder]+360]
                point3x = point1[0]-2*np.divide(np.abs(point1[0]-orderAngle[finder2]),4)
                point3y = (orderAngle2[finder]+360)+2*np.divide(np.abs((orderAngle2[finder]+360)-orderAngle2[finder2]),3)
            if qq == 3:
                point1 = [orderAngle[finder], orderAngle2[finder]+360]
                point3x = point1[0]-1*np.divide(np.abs(point1[0]-orderAngle[finder2]),10)
                point3y = (orderAngle2[finder]+360)-2*np.divide(np.abs((orderAngle2[finder]+360)-orderAngle2[finder2]),3)
            elif qq == 6:
                point2 = [orderAngle[finder2]+360,orderAngle2[finder2]+360]
                point3x = point1[0]+np.divide(np.abs(point1[0]-(orderAngle[finder2]+360)),3)
                point3y = (orderAngle[finder])+80#*np.divide(np.abs((orderAngle[finder]+360)-(orderAngle2[finder2]+360)),2)

            elif qq == 7:
                point2 = [orderAngle[finder2]+360,orderAngle2[finder2]+360]
                point3x = point1[0]+np.divide(np.abs(point1[0]-(orderAngle[finder2]+360)),3)
                point3y = (orderAngle[finder])+2*np.divide(np.abs((orderAngle[finder])-(orderAngle2[finder2]+360)),3)+30

            print('working on {}'.format(qq))
        #x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0],point1[1], point3x, point3y, point2[0],point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    elif ii == 6:
        goto = np.array((3,7))
        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder], orderAngle2[finder] + 0]
            point2 = [orderAngle[finder2], orderAngle2[finder2]]
            if qq == 2:
                point1 = [orderAngle[finder], orderAngle2[finder] + 360]
                point3x = point1[0] - 2 * np.divide(np.abs(point1[0] - orderAngle[finder2]), 4)
                point3y = (orderAngle2[finder] + 360) + 2 * np.divide(
                    np.abs((orderAngle2[finder] + 360) - orderAngle2[finder2]), 3)
            if qq == 3:
                point1 = [orderAngle[finder], orderAngle2[finder] + 0]
                point2 = [orderAngle[finder2]-360, orderAngle2[finder2]]

                point3x = -190#point1[0] + 1 * np.divide(np.abs(point1[0] - orderAngle[finder2]), 10)
                point3y = -60#(orderAngle2[finder]-360) - 1 * np.divide(np.abs((orderAngle2[finder] -360) - orderAngle2[finder2]), 3)
            elif qq == 6:
                point2 = [orderAngle[finder2] + 360, orderAngle2[finder2] + 360]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 360)), 3)
                point3y = (orderAngle[
                    finder]) + 80  # *np.divide(np.abs((orderAngle[finder]+360)-(orderAngle2[finder2]+360)),2)

            elif qq == 7:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30

            print('working on {}'.format(qq))
            # x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0], point1[1], point3x, point3y, point2[0], point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    elif ii == 7:
        goto = np.array((8,9))
        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder], orderAngle2[finder] + 0]
            point2 = [orderAngle[finder2], orderAngle2[finder2]]
            if qq == 8:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            elif qq == 9:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3)
            print('working on {}'.format(qq))
            # x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0], point1[1], point3x, point3y, point2[0], point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    elif ii == 8:
        goto = np.array((0,10))
        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder], orderAngle2[finder] + 0]
            point2 = [orderAngle[finder2], orderAngle2[finder2]]
            if qq == 0:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            elif qq == 10:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            print('working on {}'.format(qq))
            # x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0], point1[1], point3x, point3y, point2[0], point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    elif ii == 9:
        # goto = np.array((1,3,8))
        goto = np.array((1,8))

        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder], orderAngle2[finder] + 0]
            point2 = [orderAngle[finder2], orderAngle2[finder2]]
            if qq == 1:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            elif qq == 3:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            elif qq == 8:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            print('working on {}'.format(qq))
            # x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0], point1[1], point3x, point3y, point2[0], point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)

    elif ii == 10:
        # goto = np.array((0,1,2,8))
        goto = np.array((0,2))

        for qq in goto:
            finder = order[ii]
            finder2 = order[qq]
            point1 = [orderAngle[finder], orderAngle2[finder] + 0]
            point2 = [orderAngle[finder2], orderAngle2[finder2]]
            if qq == 0:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle2[finder]) + 2 * np.divide(
                    np.abs((orderAngle2[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            elif qq == 1:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            elif qq == 2:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            elif qq == 8:
                point2 = [orderAngle[finder2], orderAngle2[finder2]]
                point3x = point1[0] + np.divide(np.abs(point1[0] - (orderAngle[finder2] + 0)), 3)
                point3y = (orderAngle[finder]) + 2 * np.divide(
                    np.abs((orderAngle[finder]) - (orderAngle2[finder2] +0)), 3) + 30
            print('working on {}'.format(qq))
            # x,y = hanging_line(point1, point2)
            lista1 = curva_b(point1[0], point1[1], point3x, point3y, point2[0], point2[1])
            flat_list2 = [item for sublist in lista1 for item in sublist]
            x = np.array(flat_list2[::2])
            y = np.array(flat_list2[1::2])
            # xy = *zip(*lista1)
            if flatarray[ii,qq] < .1:
                sizeLine = smallSize
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                sizeLine = midSize
            else:
                sizeLine = largeSize
            line = ax10.plot(x, y, color=colors[ii, :], linewidth=sizeLine)
            if flatarray[ii,qq] < .1:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize1)
            elif flatarray[ii,qq] > .0999 and flatarray[ii,qq] < 0.19999:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize2)
            else:
                add_arrow(line[0], color=colors[ii, :], size=arrowSize3)
    else:
        print('skipping {}'.format(ii))

plt.tight_layout()


probMultiplier = 200
for ii in range(11): #range(numClusters):
    finder = order[ii]
    #index = np.where(yc == finder)
    if flatarray[ii,ii] < 0.55:
        sizeMarker = 50
    elif flatarray[ii,ii] < 0.70 and flatarray[ii,ii]>0.549:
        sizeMarker = 125
    else:
        sizeMarker = 200

    if ii == 5:
        ax10.scatter(orderAngle[finder],orderAngle2[finder]+360,sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)
    elif ii == 6:
        ax10.scatter(orderAngle[finder],orderAngle2[finder],sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)
        ax10.scatter(orderAngle[finder]+360,orderAngle2[finder]+360,sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)
    elif ii == 7:
        ax10.scatter(orderAngle[finder],orderAngle2[finder],sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)
        ax10.scatter(orderAngle[finder]+360,orderAngle2[finder]+360,sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)
    elif ii == 3:
        ax10.scatter(orderAngle[finder],orderAngle2[finder],sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)
        ax10.scatter(orderAngle[finder]-360,orderAngle2[finder],sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)
    elif ii == 1:
        ax10.scatter(orderAngle[finder],orderAngle2[finder],sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)
        ax10.scatter(orderAngle[finder],orderAngle2[finder]+360,sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)
    else:
        ax10.scatter(orderAngle[finder],orderAngle2[finder],sizeMarker,color=colors[ii,:],edgecolor='k',zorder=3)

#
# clusterPickle = 'moveInfoOver.pickle'
# output = {}
# output['orderAngle'] = orderAngle
# output['orderAngle2'] = orderAngle2
# output['colors'] = colors
# output['flatarray'] = flatarray
# output['order'] = order
#
# import pickle
# with open(clusterPickle,'wb') as f:
#     pickle.dump(output, f)


ax10.set_xlim([-400, 400])
ax10.set_ylim([-350, 480])

ax10.set_xticks([-360,-270,-180,-90,0,90,180,270,360])
ax10.set_xticklabels([r'-2$\pi$',r'-$\frac{3}{2}\pi$',r'-$\pi$',r'-$\frac{1}{2}\pi$','0',r'$\frac{1}{2}\pi$',r'$\pi$',r'$\frac{3}{2}\pi$',r'2$\pi$',])
ax10.set_xlabel(r'Offshore Propagation $\Longrightarrow$', fontsize=18)

ax10.set_yticks([-270,-180,-90,0,90,180,270,360])
ax10.set_yticklabels([r'-$\frac{3}{2}\pi$',r'-$\pi$',r'-$\frac{1}{2}\pi$','0',r'$\frac{1}{2}\pi$',r'$\pi$',r'$\frac{3}{2}\pi$',r'2$\pi$',])
ax10.set_ylabel(r'Onshore Propagation $\Longrightarrow$', fontsize=18)

plt.tight_layout()

left, bottom, width, height = [0.1, 0.13, 0.18, 0.18]
ax2sub = fig.add_axes([left, bottom, width, height])
ii = 6
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii,linewidth=2)
ii = 7
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='--')
ii = 3
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='--')
ax2sub.legend()

left, bottom, width, height = [0.5, 0.13, 0.18, 0.18]
ax2sub = fig.add_axes([left, bottom, width, height])
ii = 9
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii,linewidth=2)
ii = 8
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='--')
ii = 1
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='--')
ax2sub.legend()


left, bottom, width, height = [0.2, 0.73, 0.18, 0.18]
ax2sub = fig.add_axes([left, bottom, width, height])
ii = 8
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 10
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 0
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 2
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
# ii = 5
# finder = order[ii]
# ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
# ii = 7
# finder = order[ii]
# ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ax2sub.legend()




left, bottom, width, height = [0.78, 0.3, 0.18, 0.18]
ax2sub = fig.add_axes([left, bottom, width, height])
ii = 6
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 7
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 9
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 1
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 3
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ax2sub.legend()

left, bottom, width, height = [0.81, 0.8, 0.18, 0.18]
ax2sub = fig.add_axes([left, bottom, width, height])
ii = 4
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 5
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 6
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ii = 7
finder = order[ii]
ax2sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii, linestyle='-')
ax2sub.legend()


















left, bottom, width, height = [0.62, 0.83, 0.15, 0.15]
ax3sub = fig.add_axes([left, bottom, width, height])
ii = 4
finder = order[ii]
ax3sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)
ii = 5
finder = order[ii]
ax3sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)



left, bottom, width, height = [0.42, 0.8, 0.15, 0.15]
ax4sub = fig.add_axes([left, bottom, width, height])
ii = 0
finder = order[ii]
ax4sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)
ii = 2
finder = order[ii]
ax4sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)



left, bottom, width, height = [0.82, 0.44, 0.15, 0.15]
ax5sub = fig.add_axes([left, bottom, width, height])
ii = 5
finder = order[ii]
ax5sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)
ii = 6
finder = order[ii]
ax5sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)
# ii = 7
# finder = order[ii]
# ax5sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)


left, bottom, width, height = [0.06, 0.06, 0.15, 0.15]
ax5sub = fig.add_axes([left, bottom, width, height])
ii = 6
finder = order[ii]
ax5sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)
ii = 7
finder = order[ii]
ax5sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)



left, bottom, width, height = [0.12, 0.63, 0.15, 0.15]
ax5sub = fig.add_axes([left, bottom, width, height])
ii = 8
finder = order[ii]
ax5sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)
ii = 9
finder = order[ii]
ax5sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)


left, bottom, width, height = [0.3, 0.08, 0.15, 0.15]
ax7sub = fig.add_axes([left, bottom, width, height])
ii = 8
finder = order[ii]
ax7sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)
ii = 9
finder = order[ii]
ax7sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)


left, bottom, width, height = [0.22, 0.82, 0.15, 0.15]
ax6sub = fig.add_axes([left, bottom, width, height])
ii = 9
finder = order[ii]
ax6sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)
ii = 10
finder = order[ii]
ax6sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)
ii = 0
finder = order[ii]
ax6sub.plot(xinterpS, profiles[finder, :], color=colors[ii, :], label=ii)













bins = np.arange(-180,220,45)
avgMag1 = np.nan * np.ones((len(bins),len(bins)))
avgMag2 = np.nan * np.ones((len(bins),len(bins)))
stdMag1 = np.nan * np.ones((len(bins),len(bins)))
stdMag2 = np.nan * np.ones((len(bins),len(bins)))
numObs = np.nan * np.ones((len(bins),len(bins)))

for xx in range(len(bins)-1):
    for yy in range(len(bins)-1):
        index = np.where((phit2S[:,0] > bins[xx]) & (phit2S[:,0] < bins[xx+1]) & (phit2S[:,1] > bins[yy]) & (phit2S[:,1] < bins[yy+1]))
        avgMag1[yy,xx] = np.nanmean(RS[index[0],0])
        avgMag2[yy,xx] = np.nanmean(RS[index[0],1])
        stdMag1[yy,xx] = np.nanmean(RS[index[0],0])
        stdMag2[yy,xx] = np.nanmean(RS[index[0],1])
        numObs[yy,xx] = len(index[0])

meshbinsX,meshbinsY = np.meshgrid(bins,bins)

fig10 = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((2,2),(0,0),rowspan=1,colspan=1)
ax1.pcolor(meshbinsX,meshbinsY,avgMag1,vmin=0.5,vmax=1.3,cmap='bwr')
ax2 = plt.subplot2grid((2,2),(0,1),rowspan=1,colspan=1)
p1 = ax2.pcolor(meshbinsX,meshbinsY,avgMag2,vmin=0.5,vmax=1.3,cmap='bwr')
cb1 = plt.colorbar(p1,ax=ax2)

ax3 = plt.subplot2grid((2,2),(1,0),rowspan=1,colspan=1)
p3 = ax3.pcolor(meshbinsX,meshbinsY,stdMag1)#,vmin=0.3,vmax=1.3)
ax4 = plt.subplot2grid((2,2),(1,1),rowspan=1,colspan=1)
p4 = ax4.pcolor(meshbinsX,meshbinsY,stdMag2)#,vmin=0.3,vmax=1.3)
cb2 = plt.colorbar(p4,ax=ax4)

plt.figure()
plt.pcolor(meshbinsX,meshbinsY,numObs)
plt.colorbar()

plt.figure()
for ii in range(numClusters):
    finder = order[ii]
    index = np.where(yc == finder)
    avgProf = np.nanmean(alllines[index[0],:],axis=0)
    plt.plot(xinterpS, avgProf,color=colors[ii,:])
plt.title('Average of Profiles')
plt.show()


plt.figure()
ii = 5
finder = order[ii]
index = np.where(yc == finder)
avgProf = np.nanmean(alllines[index[0],:],axis=0)
plt.plot(xinterpS, avgProf,color=colors[ii,:])
plt.plot(xinterpS, profiles[finder, :], color=colors[ii, :])
plt.show()





fig3 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(4, 4, wspace=0.0, hspace=0.0)

ax = plt.subplot(gs[:, 0:2])
for ii in range(numClusters):
    finder = order[ii]
    ax.plot(xinterpS, profiles[finder,:],color=colors[ii,:],label=ii)
ax.legend()
ax.set_title('Centroids')

gr, gc = 0, 2
import matplotlib.cm as cm
# colors = cm.rainbow(np.linspace(0, 1, numClusters))
colors = cm.hsv(np.linspace(0, 1, (numClusters+1)))

for ii in range(numClusters):

    ax = plt.subplot(gs[gr, gc])

    finder = order[ii]

    index = np.where(yc == finder)
    # for qq in range(len(index[0])):
    #     ax.plot(xinterp,alllinesS[index[0][qq],:],color=[0.7, 0.7, 0.7])

    avgProf = np.nanmean(alllines[index[0],:],axis=0)
    stdProf = np.nanstd(alllines[index[0],:],axis=0)
    #ax.plot(xinterpS, avgProf,color=colors[ii,:])


    ax.plot(xinterpS, profiles[finder,:],color=colors[ii,:],label=ii,linewidth=2)
    # ax.plot(xinterpS, profiles[finder,:]+stdProf,'--',color=[0.5, 0.5, 0.5])
    # ax.plot(xinterpS, profiles[finder,:]-stdProf,'--',color=[0.5, 0.5, 0.5])
    #
    # ax.plot(xinterpS, profiles[finder,:]+2*stdProf,'--',color=[0.8, 0.8, 0.8])
    # ax.plot(xinterpS, profiles[finder,:]-2*stdProf,'--',color=[0.8, 0.8, 0.8])
    ax.plot(xinterpS, avgProf+stdProf,'--',color=[0.5, 0.5, 0.5])
    ax.plot(xinterpS, avgProf-stdProf,'--',color=[0.5, 0.5, 0.5])

    ax.plot(xinterpS, avgProf+2*stdProf,'--',color=[0.8, 0.8, 0.8])
    ax.plot(xinterpS, avgProf-2*stdProf,'--',color=[0.8, 0.8, 0.8])

    ax.set_xlim([0,500])
    ax.set_ylim([-7, 0.5])
    #ax.set_title('{}'.format(KMA.group_size.values[i]))
    #ax.text(400,0, group_size[sortedPeakInd[i]], fontweight='bold')
    #ax.text(400,0, sortedPhases[i], fontweight='bold')

    if gc > 0:
        ax.set_yticks([])

    if gr < (4):
        ax.set_xticks([])
    #  counter
    gr += 1
    if gr >= 4:
        gc += 1
        gr = 0



#
#
# madeUpOrder = np.array([ 9, 3, 0, 5, 14, 12, 7, 6, 1, 4, 13, 11,  8, 10,  2])
# madeUpOrderAngle = np.array([ 89.87321293, 244.83742932, 356.99748811,  75.14747898,
#        240.61349123, 142.77626464, 216.29014046, 170.16120042,
#        311.25367613,  43.16565455, 312.13700038, 268.6581865 ,
#        146.02703599, 332.01876024, 167.67728105])
#
# plt.figure()
# ax1 = plt.subplot2grid((6,12),(0,0),rowspan=3,colspan=6)
# for ii in range(numClusters):
#     finder = order[ii]
#     index = np.where(yc == finder)
#     avgProf = np.nanmean(alllines[index[0],:],axis=0)
#     # if orderMag[finder] > 0.67:
#
#     # ax1.plot(xinterpS, profiles[finder,:],color=colors[ii,:],label=np.round(orderAngle[finder]))
#     ax1.plot(xinterpS, avgProf,color=colors[ii,:],linewidth=orderMag[finder], label= "{} + {}".format(np.round(orderAngle[finder]), np.round(orderAngle2[finder])))
#
# ax1.legend()
# ax1.set_title('Centroids')
#
# ax2 = plt.subplot2grid((6,12),(3,0),rowspan=3,colspan=3)
# ids = np.array([6, 7, 8, 9, 10, 11, 12])
# for ii in ids:
#     finder = order[ii]
#     index = np.where(yc == finder)
#     avgProf = np.nanmean(alllines[index[0],:],axis=0)
#     if orderMag[finder] > 0.0:
#
#     # ax2.plot(xinterpS, profiles[finder,:],color=colors[ii,:],label=np.round(orderAngle[finder]))
#         ax2.plot(xinterpS, avgProf,color=colors[ii,:],linewidth=orderMag[finder], label= "{} + {}".format(np.round(orderAngle[finder]), np.round(orderAngle2[finder])))
#
# ax2.legend()
# ax2.set_title('Emergence of New Bar')
#
# ax3 = plt.subplot2grid((6,12),(3,3),rowspan=3,colspan=3)
# ids = np.array([0,1,2,3,4,5,6,7,8])
# for ii in ids:
#     finder = order[ii]
#     index = np.where(yc == finder)
#     avgProf = np.nanmean(alllines[index[0],:],axis=0)
#     # ax3.plot(xinterpS, profiles[finder,:],color=colors[ii,:],label=np.round(orderAngle[finder]))
#     if orderMag[finder] > 0.0:
#         ax3.plot(xinterpS, avgProf,color=colors[ii,:],linewidth=orderMag[finder], label= "{} + {}".format(np.round(orderAngle[finder]), np.round(orderAngle2[finder])))
#
# ax3.legend()
# ax3.set_title('Rapid Migration')
#
# ax = plt.subplot2grid((6,12),(0,7),colspan=3,rowspan=6)
#
# for ii in range(numClusters):
#     finder = order[ii]
#     index = np.where(yc == finder)
#     avgProf = np.nanmean(alllines[index[0],:],axis=0)
#     ax.plot((ii+1)*np.ones((len(index[0]))),timeS[index[0]],'o',color=colors[ii,:])
# ax.set_ylim([time[0], time[150]])
#
# ax4 = plt.subplot2grid((6,12),(0,10),colspan=3,rowspan=6)
#
# # ax4.set_cmap('RdBu_r')
#
# tg, xg = np.meshgrid(timeS, xinterpS)
# plt4 = ax4.pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8, cmap='RdBu_r')
# # plt.colorbar(plt4, ax=ax4, orientation='horizontal')
# ax4.set_ylim([time[0], time[150]])
# # ax[0].set_title('Surveys (dev.)')
#
#





#
#
# plt.figure()
# ax1 = plt.subplot2grid((3,3),(0,0),colspan=2,rowspan=3)
# for ii in range(numClusters):
#     finder = order[ii]
#     ax1.plot(xinterpS, profiles[finder,:],color=colors[ii,:],label=(ii+1))
# ax1.legend()
# ax1.set_title('Centroids')
# ax1.set_xlabel('Cross-shore (m)')
# ax1.set_ylabel('Depth (m relative to MHHW)')
#
# ax = plt.subplot2grid((3,3),(0,2),colspan=1,rowspan=3)
#
# for ii in range(numClusters):
#     finder = order[ii]
#     index = np.where(yc == finder)
#     avgProf = np.nanmean(alllines[index[0],:],axis=0)
#     ax.plot((ii+1)*np.ones((len(index[0]))),timeS[index[0]],'o',color=colors[ii,:])
#
# t1 = datetime.datetime(1980,5,1)
# t2 = datetime.datetime(1984,2,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.65, 0.65, 0.65])
# ax.add_patch(rect)
#
# t1 = datetime.datetime(1984,2,1)
# t2 = datetime.datetime(1988,2,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.8, 0.8, 0.8])
# ax.add_patch(rect)
#
# t1 = datetime.datetime(1988,2,1)
# t2 = datetime.datetime(1990,8,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.65, 0.65, 0.65])
# ax.add_patch(rect)
#
# t1 = datetime.datetime(1990,8,1)
# t2 = datetime.datetime(1993,2,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.8, 0.8, 0.8])
# ax.add_patch(rect)
#
# t1 = datetime.datetime(1993,2,1)
# t2 = datetime.datetime(1999,5,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.65, 0.65, 0.65])
# ax.add_patch(rect)
#
# t1 = datetime.datetime(1999,5,1)
# t2 = datetime.datetime(2004,10,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.8, 0.8, 0.8])
# ax.add_patch(rect)
#
# t1 = datetime.datetime(2004,10,1)
# t2 = datetime.datetime(2010,5,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.65, 0.65, 0.65])
# ax.add_patch(rect)
#
# t1 = datetime.datetime(2010,5,1)
# t2 = datetime.datetime(2014,10,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.8, 0.8, 0.8])
# ax.add_patch(rect)
#
# t1 = datetime.datetime(2014,10,1)
# t2 = datetime.datetime(2018,5,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.65, 0.65, 0.65])
# ax.add_patch(rect)
#
# t1 = datetime.datetime(2018,5,1)
# t2 = datetime.datetime(2020,9,1)
# start = mdates.date2num(t1)
# end = mdates.date2num(t2)
# width = end - start
# rect = Rectangle((0.5, start), 8.5, width, color=[0.8, 0.8, 0.8])
# ax.add_patch(rect)
# ax.set_xticks(np.arange(1, 9, 1))
# ax.set_xlabel('Cluster Number')
#
#
#
# t1 = datetime.datetime(1997,8,15)
# t2 = datetime.datetime(1997,10,30)
# sandyduck = np.where((timeS > t1) & (timeS < t2))
# xs = [] #np.empty(len(sandyduck[0],))
# ys = [] #np.empty(len(sandyduck[0],))
#
# for ii in range(len(sandyduck[0])):
#     finder = np.where((yc[sandyduck[0][ii]] == order))
#     ax.plot(finder[0]+1, timeS[sandyduck[0][ii]],'ko-',markersize=8,fillstyle='none') #, color="red", fillstyle="full") #, color=colors[ii, :])
#     xs.append(finder[0]+1)
#     ys.append(timeS[sandyduck[0][ii]])
# ax.plot(xs, ys,'k-') #, color="red", fillstyle="full") #, color=colors[ii, :])
# ax.text(10,ys[0],'SandyDuck97')
#
# t1 = datetime.datetime(1994,9,1)
# t2 = datetime.datetime(1994,10,30)
# sandyduck = np.where((timeS > t1) & (timeS < t2))
# xs = [] #np.empty(len(sandyduck[0],))
# ys = [] #np.empty(len(sandyduck[0],))
#
# for ii in range(len(sandyduck[0])):
#     finder = np.where((yc[sandyduck[0][ii]] == order))
#     ax.plot(finder[0]+1, timeS[sandyduck[0][ii]],'ko-',markersize=8,fillstyle='none') #, color="red", fillstyle="full") #, color=colors[ii, :])
#     xs.append(finder[0]+1)
#     ys.append(timeS[sandyduck[0][ii]])
# ax.plot(xs, ys,'k-') #, color="red", fillstyle="full") #, color=colors[ii, :])
# ax.text(10,ys[0],'Duck94')
#
#
# t1 = datetime.datetime(1990,9,20)
# t2 = datetime.datetime(1990,11,2)
# sandyduck = np.where((timeS > t1) & (timeS < t2))
# xs = [] #np.empty(len(sandyduck[0],))
# ys = [] #np.empty(len(sandyduck[0],))
#
# for ii in range(len(sandyduck[0])):
#     finder = np.where((yc[sandyduck[0][ii]] == order))
#     ax.plot(finder[0]+1, timeS[sandyduck[0][ii]],'ko-',markersize=8,fillstyle='none') #, color="red", fillstyle="full") #, color=colors[ii, :])
#     xs.append(finder[0]+1)
#     ys.append(timeS[sandyduck[0][ii]])
# ax.plot(xs, ys,'k-') #, color="red", fillstyle="full") #, color=colors[ii, :])
# ax.text(10,ys[0],'DELILAH')
#
#
# t1 = datetime.datetime(1986,9,1)
# t2 = datetime.datetime(1986,10,30)
# sandyduck = np.where((timeS > t1) & (timeS < t2))
# xs = [] #np.empty(len(sandyduck[0],))
# ys = [] #np.empty(len(sandyduck[0],))
#
# for ii in range(len(sandyduck[0])):
#     finder = np.where((yc[sandyduck[0][ii]] == order))
#     ax.plot(finder[0]+1, timeS[sandyduck[0][ii]],'ko-',markersize=8,fillstyle='none') #, color="red", fillstyle="full") #, color=colors[ii, :])
#     xs.append(finder[0]+1)
#     ys.append(timeS[sandyduck[0][ii]])
# ax.plot(xs, ys,'k-') #, color="red", fillstyle="full") #, color=colors[ii, :])
# ax.text(10,ys[0],'SUPERDUCK')
#
# t1 = datetime.datetime(2019,9,15)
# t2 = datetime.datetime(2019,10,30)
# sandyduck = np.where((timeS > t1) & (timeS < t2))
# xs = [] #np.empty(len(sandyduck[0],))
# ys = [] #np.empty(len(sandyduck[0],))
#
# for ii in range(len(sandyduck[0])):
#     finder = np.where((yc[sandyduck[0][ii]] == order))
#     ax.plot(finder[0]+1, timeS[sandyduck[0][ii]],'ko-',markersize=8,fillstyle='none') #, color="red", fillstyle="full") #, color=colors[ii, :])
#     xs.append(finder[0]+1)
#     ys.append(timeS[sandyduck[0][ii]])
# ax.plot(xs, ys,'k-') #, color="red", fillstyle="full") #, color=colors[ii, :])
# # ax.arrow(xs[1][0], ys[1], (xs[2][0]-xs[1][0]), (ys[2]-ys[1])) #, fc="k", ec="k",head_width=0.05, head_length=0.1 )
# ax.text(10,ys[0],'DUNEX Pilot')
#
# # t1 = datetime.datetime(1997,9,1)
# # t2 = datetime.datetime(1997,10,30)
# # start = mdates.date2num(t1)
# # end = mdates.date2num(t2)
# # width = end - start
# # rect = Rectangle((6.5, start), 8.5, width, color=[1, 0, 0])
# # ax.add_patch(rect)
# #
# # t1 = datetime.datetime(1994,9,1)
# # t2 = datetime.datetime(1994,10,30)
# # start = mdates.date2num(t1)
# # end = mdates.date2num(t2)
# # width = end - start
# # rect = Rectangle((6.5, start), 8.5, width, color=[1, 0, 0])
# # ax.add_patch(rect)
# #
# # t1 = datetime.datetime(1990,10,1)
# # t2 = datetime.datetime(1990,10,30)
# # start = mdates.date2num(t1)
# # end = mdates.date2num(t2)
# # width = end - start
# # rect = Rectangle((6.5, start), 8.5, width, color=[1, 0, 0])
# # ax.add_patch(rect)
# #
# # t1 = datetime.datetime(1986,9,1)
# # t2 = datetime.datetime(1986,10,30)
# # start = mdates.date2num(t1)
# # end = mdates.date2num(t2)
# # width = end - start
# # rect = Rectangle((6.5, start), 8.5, width, color=[1, 0, 0])
# # ax.add_patch(rect)
# # for xx in indices:
# #     t1 = timeList[xx][0]
# #     t2 = timeList[xx][-1]
# #     # # convert to matplotlib date representation
# #
# #
# #     # Plot rectangle
# #     qq = dwtList[xx]
# #     if qq > 30:
# #
# #     else:
# #         rect = Rectangle((start, 0), width, 1, color=dwtcolors[qq-1])
# #
#
# # ax4.set_xlim([datetime(2011,6,1), datetime(2011,11,1)])
#
#

asdfg







clusterPickle = 'sandbarsSouthernTransectPythonClusters8.pickle'
output = {}
output['time'] = timeS
output['alllines'] = alllines
output['xinterp'] = xinterp
output['order'] = order
output['profiles'] = profiles
output['bmu'] = yc
output['numCenters'] = numClusters
import pickle
with open(clusterPickle,'wb') as f:
    pickle.dump(output, f)






#
# fig1 = plt.figure()
# ax1 = plt.subplot2grid((5,5),(0,2),rowspan=1,colspan=1)
# ax1.plot(xinterpS,np.nanmean(alllines,axis=0))
#
# ax2 = plt.subplot2grid((5,5),(1,1),rowspan=1,colspan=1)
# finder = np.where((bmu2 ==0))
# ax2.plot(xinterpS,np.nanmean(alllines[finder[0],:],axis=0))
# ax3 = plt.subplot2grid((5,5),(1,3),rowspan=1,colspan=1)
# finder = np.where((bmu2 ==1))
# ax3.plot(xinterpS,np.nanmean(alllines[finder[0],:],axis=0))
#
# ax4 = plt.subplot2grid((5,5),(2,1),rowspan=1,colspan=1)
# finder = np.where((bmu3 ==0))
# ax4.plot(xinterpS,np.nanmean(alllines[finder[0],:],axis=0))
# ax5 = plt.subplot2grid((5,5),(2,3),rowspan=1,colspan=1)
# finder = np.where((bmu3 ==1))
# ax5.plot(xinterpS,np.nanmean(alllines[finder[0],:],axis=0))
# ax6 = plt.subplot2grid((5,5),(2,4),rowspan=1,colspan=1)
# finder = np.where((bmu3 ==2))
# ax6.plot(xinterpS,np.nanmean(alllines[finder[0],:],axis=0))
#
# ax7 = plt.subplot2grid((5,5),(3,1),rowspan=1,colspan=1)
# finder = np.where((bmu4 ==0))
# ax7.plot(xinterpS,np.nanmean(alllines[finder[0],:],axis=0))
# ax8 = plt.subplot2grid((5,5),(3,2),rowspan=1,colspan=1)
# finder = np.where((bmu4 ==1))
# ax8.plot(xinterpS,np.nanmean(alllines[finder[0],:],axis=0))
# ax9 = plt.subplot2grid((5,5),(3,3),rowspan=1,colspan=1)
# finder = np.where((bmu4 ==2))
# ax9.plot(xinterpS,np.nanmean(alllines[finder[0],:],axis=0))
# ax10 = plt.subplot2grid((5,5),(3,4),rowspan=1,colspan=1)
# finder = np.where((bmu4 ==3))
# ax10.plot(xinterpS,np.nanmean(alllines[finder[0],:],axis=0))
#


# negindex = np.where((phitS<0))
# phitS[negindex] = phitS[negindex]+2*np.pi
# color=colors[ii,:]
fig1 = plt.figure(figsize=(18,10))

# 1190 is all members (596)
ax1 = plt.subplot2grid((8,16),(0,7),rowspan=2,colspan=2)
ax1.plot(xinterpS,np.nanmean(alllines[members[1190]],axis=0),'k-')
ax1.plot(xinterpS,np.nanmean(alllines[members[1190]],axis=0)+np.nanstd(alllines[members[1190]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax1.plot(xinterpS,np.nanmean(alllines[members[1190]],axis=0)-np.nanstd(alllines[members[1190]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax1.set_title('All Profiles',size=14)
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.text(400,-1,len(members[1190]))

## 1189 = 334
## 1187 = 262
##  334 + 262 = 596
ax2 = plt.subplot2grid((8,16),(2,3),rowspan=2,colspan=2)
ax2.plot(xinterpS,np.nanmean(alllines[members[1189]],axis=0),'k-')
ax2.plot(xinterpS,np.nanmean(alllines[members[1189]],axis=0)+np.nanstd(alllines[members[1189]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax2.plot(xinterpS,np.nanmean(alllines[members[1189]],axis=0)-np.nanstd(alllines[members[1189]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax2.set_title('Inner-Bar',size=14)
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.text(400,-1,len(members[1189]))

bbox_args = dict(boxstyle="round", fc="0.8")
arrow_args = dict(arrowstyle="->")

# Here we'll demonstrate the extents of the coordinate system and how
# we place annotating text.

# ax1.annotate('figure fraction : 0, 0', xy=(0.5, 0.5), xycoords='figure fraction',
#              xytext=(20, 20), textcoords='offset points',
#              ha="left", va="bottom",
#              bbox=bbox_args,
#              arrowprops=arrow_args)




ax3 = plt.subplot2grid((8,16),(2,11),rowspan=2,colspan=2)
ax3.plot(xinterpS,np.nanmean(alllines[members[1187]],axis=0),'k-')
ax3.plot(xinterpS,np.nanmean(alllines[members[1187]],axis=0)+np.nanstd(alllines[members[1187]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax3.plot(xinterpS,np.nanmean(alllines[members[1187]],axis=0)-np.nanstd(alllines[members[1187]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax3.set_title('Outer-Bar',size=14)
ax3.set_xticklabels([])
ax3.set_yticklabels([])
ax3.text(400,-1,len(members[1187]))

## 1188 = 225
## 1184 = 109
# 225 + 109 = 334
ax4 = plt.subplot2grid((8,16),(4,5),rowspan=2,colspan=2)
ax4.plot(xinterpS,np.nanmean(alllines[members[1188]],axis=0),'k-')
ax4.plot(xinterpS,np.nanmean(alllines[members[1188]],axis=0)+np.nanstd(alllines[members[1188]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax4.plot(xinterpS,np.nanmean(alllines[members[1188]],axis=0)-np.nanstd(alllines[members[1188]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax4.set_title('Prominent Inner-Bar',size=14)
ax4.set_xticklabels([])
ax4.set_yticklabels([])
ax4.text(400,-1,len(members[1188]))

ax5 = plt.subplot2grid((8,16),(4,1),rowspan=2,colspan=2)
ax5.plot(xinterpS,np.nanmean(alllines[members[1184]],axis=0),'k-')
ax5.plot(xinterpS,np.nanmean(alllines[members[1184]],axis=0)+np.nanstd(alllines[members[1184]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax5.plot(xinterpS,np.nanmean(alllines[members[1184]],axis=0)-np.nanstd(alllines[members[1184]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax5.set_xticklabels([])
ax5.set_yticklabels([])
ax5.set_title('Inner-Berm',size=14)
ax5.text(400,-1,len(members[1184]))

## 1185 = 131
## 1183 = 131
# 131 + 131 = 262
ax6 = plt.subplot2grid((8,16),(4,9),rowspan=2,colspan=2)
ax6.plot(xinterpS,np.nanmean(alllines[members[1185]],axis=0),'k-')
ax6.plot(xinterpS,np.nanmean(alllines[members[1185]],axis=0)+np.nanstd(alllines[members[1185]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax6.plot(xinterpS,np.nanmean(alllines[members[1185]],axis=0)-np.nanstd(alllines[members[1185]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax6.set_xticklabels([])
ax6.set_yticklabels([])
ax6.set_title('Prominent Outer-Bar',size=14)
ax6.text(400,-1,len(members[1185]))

ax7 = plt.subplot2grid((8,16),(4,13),rowspan=2,colspan=2)
ax7.plot(xinterpS,np.nanmean(alllines[members[1183]],axis=0),'k-')
ax7.plot(xinterpS,np.nanmean(alllines[members[1183]],axis=0)+np.nanstd(alllines[members[1183]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax7.plot(xinterpS,np.nanmean(alllines[members[1183]],axis=0)-np.nanstd(alllines[members[1183]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax7.set_xticklabels([])
ax7.set_yticklabels([])
ax7.set_title('Double-Bar',size=14)
ax7.text(400,-1,len(members[1183]))


## 1186 = 129
## 1182 = 96
# 129 + 96 = 225
ax8 = plt.subplot2grid((8,16),(6,4),rowspan=2,colspan=2)
ax8.plot(xinterpS,np.nanmean(alllines[members[1186]],axis=0))
ax8.plot(xinterpS,np.nanmean(alllines[members[1186]],axis=0)+np.nanstd(alllines[members[1186]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax8.plot(xinterpS,np.nanmean(alllines[members[1186]],axis=0)-np.nanstd(alllines[members[1186]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax8.set_xticklabels([])
ax8.set_yticklabels([])
temp8 = phitS[members[1186],0]
indT = np.where((temp8<-(np.pi/2+.3)))
temp8[indT] = temp8[indT]+2*np.pi
temp8b = phitS[members[1186],1]
indT = np.where((temp8b<0))
temp8b[indT] = temp8b[indT]+2*np.pi
# ax8.set_title(['EOF1 = ', np.round(100*np.median(np.real(RS[members[1186],0])))/100, ' ', np.round(100*np.median(phitS[members[1186],0])*180/np.pi)/100])
# ax8.set_title('EOF1 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1186],0])))/100, np.round(np.median(temp8)*180/np.pi)))
# ax8.set_title('EOF1 = {:.2f} @ {:.0f} \nEOF2 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1186],0])))/100, np.round(np.median(temp8)*180/np.pi), np.round(100*np.median(np.real(RS[members[1186],1])))/100, np.round(np.median(temp8b)*180/np.pi)))

ax8.text(400,-1,len(members[1186]))
# phitS[:,0],RS[:,0]
ax9 = plt.subplot2grid((8,16),(6,6),rowspan=2,colspan=2)
ax9.plot(xinterpS,np.nanmean(alllines[members[1182]],axis=0))
ax9.plot(xinterpS,np.nanmean(alllines[members[1182]],axis=0)+np.nanstd(alllines[members[1182]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax9.plot(xinterpS,np.nanmean(alllines[members[1182]],axis=0)-np.nanstd(alllines[members[1182]],axis=0),'--',color=[0.5, 0.5, 0.5])
temp7 = phitS[members[1182],0]
temp7b = phitS[members[1182],1]
indT = np.where((temp7b<-(np.pi/2-.2)))
temp7b[indT] = temp7b[indT]+2*np.pi
# ax9.set_title('EOF1 = {:.2f} @ {:.0f} \nEOF2 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1182],0])))/100, np.round(np.median(temp7)*180/np.pi), np.round(100*np.median(np.real(RS[members[1182],1])))/100, np.round(np.median(temp7b)*180/np.pi)))

# ax9.set_title('EOF1 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1182],0])))/100, np.round(np.median(temp7)*180/np.pi)))
ax9.set_xticklabels([])
ax9.set_yticklabels([])
ax9.text(400,-1,len(members[1182]))

## 1179 = 40
## 1177 = 69
# 40 + 69 = 109
ax10 = plt.subplot2grid((8,16),(6,0),rowspan=2,colspan=2)
ax10.plot(xinterpS,np.nanmean(alllines[members[1179]],axis=0))
ax10.plot(xinterpS,np.nanmean(alllines[members[1179]],axis=0)+np.nanstd(alllines[members[1179]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax10.plot(xinterpS,np.nanmean(alllines[members[1179]],axis=0)-np.nanstd(alllines[members[1179]],axis=0),'--',color=[0.5, 0.5, 0.5])
temp6 = phitS[members[1179],0]
indT = np.where((temp6<0))
temp6[indT] = temp6[indT]+2*np.pi
# ax10.set_title('EOF1 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1179],0])))/100, np.round(np.median(temp6)*180/np.pi)))
temp6b = phitS[members[1179],1]
indTb = np.where((temp6b<0))
temp6b[indTb] = temp6b[indTb]+2*np.pi
# ax10.set_title('EOF1 = {:.2f} @ {:.0f} \nEOF2 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1179],0])))/100, np.round(np.median(temp6)*180/np.pi), np.round(100*np.median(np.real(RS[members[1179],1])))/100, np.round(np.median(temp6b)*180/np.pi)))

ax10.set_xticklabels([])
ax10.set_yticklabels([])
ax10.text(400,-1,len(members[1179]))

ax11 = plt.subplot2grid((8,16),(6,2),rowspan=2,colspan=2)
ax11.plot(xinterpS,np.nanmean(alllines[members[1177]],axis=0))
ax11.plot(xinterpS,np.nanmean(alllines[members[1177]],axis=0)+np.nanstd(alllines[members[1177]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax11.plot(xinterpS,np.nanmean(alllines[members[1177]],axis=0)-np.nanstd(alllines[members[1177]],axis=0),'--',color=[0.5, 0.5, 0.5])
temp5 = phitS[members[1177],0]
indT = np.where((temp5<0))
temp5[indT] = temp5[indT]+2*np.pi
temp5b = phitS[members[1177],1]
temp5b = temp5b+2*np.pi
# ax11.set_title('EOF1 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1177],0])))/100, np.round(np.median(temp5)*180/np.pi)))
# ax11.set_title('EOF1 = {:.2f} @ {:.0f} \nEOF2 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1177],0])))/100, np.round(np.median(temp5)*180/np.pi), np.round(100*np.median(np.real(RS[members[1177],1])))/100, np.round(np.median(temp5b)*180/np.pi)))

ax11.set_xticklabels([])
ax11.set_yticklabels([])
ax11.text(400,-1,len(members[1177]))

## 1181 = 109
## 1180 = 61
## 1178 = 62
## 1176 = 67
## 1175 = 65
## 1174 = 46
## 1173 = 44
## 1172 = 39
## 1171 = 70
## 1165 = 22

ax12 = plt.subplot2grid((8,16),(6,8),rowspan=2,colspan=2)
ax12.plot(xinterpS,np.nanmean(alllines[members[1181]],axis=0))
ax12.plot(xinterpS,np.nanmean(alllines[members[1181]],axis=0)+np.nanstd(alllines[members[1181]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax12.plot(xinterpS,np.nanmean(alllines[members[1181]],axis=0)-np.nanstd(alllines[members[1181]],axis=0),'--',color=[0.5, 0.5, 0.5])
temp = phitS[members[1181],0]
indT = np.where((temp<0))
temp[indT] = temp[indT]+2*np.pi
tempb = phitS[members[1181],1]
indTb = np.where((tempb<0))
tempb[indTb] = tempb[indTb]+2*np.pi

# ax12.set_title('EOF1 = {:.2f} @ {:.0f} \nEOF2 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1181],0])))/100, np.round(np.median(temp)*180/np.pi), np.round(100*np.median(np.real(RS[members[1181],1])))/100, np.round(np.median(tempb)*180/np.pi)))

# ax12.set_title('EOF1 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1181],0])))/100, np.round(np.median(temp)*180/np.pi)))
ax12.set_xticklabels([])
ax12.set_yticklabels([])
ax12.text(400,-1,len(members[1181]))


ax13 = plt.subplot2grid((8,16),(6,10),rowspan=2,colspan=2)
ax13.plot(xinterpS,np.nanmean(alllines[members[1165]],axis=0))
ax13.plot(xinterpS,np.nanmean(alllines[members[1165]],axis=0)+np.nanstd(alllines[members[1165]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax13.plot(xinterpS,np.nanmean(alllines[members[1165]],axis=0)-np.nanstd(alllines[members[1165]],axis=0),'--',color=[0.5, 0.5, 0.5])
temp2 = phitS[members[1165],0]
indT = np.where((temp2<0))
temp2[indT] = temp2[indT]+2*np.pi
temp2b = phitS[members[1165],1]
indTb = np.where((temp2b<0))
temp2b[indTb] = temp2b[indTb]+2*np.pi
# ax13.set_title('EOF1 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1165],0])))/100, np.round(np.median(temp2)*180/np.pi)))
# ax13.set_title('EOF1 = {:.2f} @ {:.0f} \nEOF2 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1165],0])))/100, np.round(np.median(temp2)*180/np.pi), np.round(100*np.median(np.real(RS[members[1171],1])))/100, np.round(np.median(temp2b)*180/np.pi)))
ax13.set_xticklabels([])
ax13.set_yticklabels([])
ax13.text(400,-1,len(members[1165]))

ax14 = plt.subplot2grid((8,16),(6,12),rowspan=2,colspan=2)
ax14.plot(xinterpS,np.nanmean(alllines[members[1180]],axis=0))
ax14.plot(xinterpS,np.nanmean(alllines[members[1180]],axis=0)+np.nanstd(alllines[members[1180]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax14.plot(xinterpS,np.nanmean(alllines[members[1180]],axis=0)-np.nanstd(alllines[members[1180]],axis=0),'--',color=[0.5, 0.5, 0.5])
temp4 = phitS[members[1180],0]
indT = np.where((temp4<0))
temp4[indT] = temp4[indT]+2*np.pi
temp4b = phitS[members[1180],1]
indT = np.where((temp4b<0))
temp4b[indT] = temp4b[indT]+2*np.pi
indT = np.where((temp4b<(np.pi/2)))
temp4b[indT] = temp4b[indT]+2*np.pi

# ax14.set_title('EOF1 = {:.2f} @ {:.0f} \nEOF2 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1180],0])))/100, np.round(np.median(temp4)*180/np.pi), np.round(100*np.median(np.real(RS[members[1180],1])))/100, np.round(np.median(temp4b)*180/np.pi)))
ax14.set_xticklabels([])
ax14.set_yticklabels([])
ax14.text(400,-1,len(members[1180]))

ax15 = plt.subplot2grid((8,16),(6,14),rowspan=2,colspan=2)
ax15.plot(xinterpS,np.nanmean(alllines[members[1171]],axis=0))
ax15.plot(xinterpS,np.nanmean(alllines[members[1171]],axis=0)+np.nanstd(alllines[members[1171]],axis=0),'--',color=[0.5, 0.5, 0.5])
ax15.plot(xinterpS,np.nanmean(alllines[members[1171]],axis=0)-np.nanstd(alllines[members[1171]],axis=0),'--',color=[0.5, 0.5, 0.5])
temp3 = phitS[members[1171],0]
indT = np.where((temp3<0))
temp3[indT] = temp3[indT]+2*np.pi
indT2 = np.where((temp3<np.pi/2))
temp3[indT2] = temp3[indT2]+2*np.pi
temp3b = phitS[members[1171],1]
# indT = np.where((temp3b<0))
# temp3b[indT] = temp3b[indT]+2*np.pi
# indT = np.where((temp4b<(np.pi/2)))
# temp4b[indT] = temp4b[indT]+2*np.pi
# ax15.set_title('EOF1 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1171],0])))/100, np.round(np.median(temp3)*180/np.pi)))
# ax15.set_title('EOF1 = {:.2f} @ {:.0f} \nEOF2 = {:.2f} @ {:.0f}'.format(np.round(100*np.median(np.real(RS[members[1171],0])))/100, np.round(np.median(temp3)*180/np.pi), np.round(100*np.median(np.real(RS[members[1171],1])))/100, np.round(np.median(temp3b)*180/np.pi)))

ax15.set_xticklabels([])
ax15.set_yticklabels([])
ax15.text(400,-1,len(members[1171]))

plt.tight_layout()

from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
    BboxConnectorPatch
from matplotlib.transforms import Bbox, TransformedBbox, \
    blended_transform_factory


# xmin = 0
# xmax = 500
# trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
# trans2 = blended_transform_factory(ax3.transData, ax3.transAxes)
# bbox = Bbox.from_extents(xmin, 0, xmax, 1)
# mybbox1 = TransformedBbox(bbox, trans1)
# mybbox2 = TransformedBbox(bbox, trans2)
mybbox1 = ax1.bbox
mybbox2 = ax2.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=3, loc2=1, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)
mybbox1 = ax1.bbox
mybbox2 = ax3.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=4, loc2=2, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)

mybbox1 = ax2.bbox
mybbox2 = ax4.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=4, loc2=2, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)
mybbox1 = ax2.bbox
mybbox2 = ax5.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=3, loc2=1, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)

mybbox1 = ax3.bbox
mybbox2 = ax6.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=3, loc2=1, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)
mybbox1 = ax3.bbox
mybbox2 = ax7.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=4, loc2=2, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)

mybbox1 = ax4.bbox
mybbox2 = ax9.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=4, loc2=1, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)

mybbox1 = ax4.bbox
mybbox2 = ax8.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=3, loc2=2, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)

mybbox1 = ax5.bbox
mybbox2 = ax10.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=3, loc2=2, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)

mybbox1 = ax5.bbox
mybbox2 = ax11.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=4, loc2=1, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)

mybbox1 = ax6.bbox
mybbox2 = ax12.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=3, loc2=2, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)

mybbox1 = ax6.bbox
mybbox2 = ax13.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=4, loc2=1, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)
mybbox1 = ax7.bbox
mybbox2 = ax14.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=3, loc2=2, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)
mybbox1 = ax7.bbox
mybbox2 = ax15.bbox
ctest = BboxConnector(mybbox1, mybbox2, loc1=4, loc2=1, linewidth=2, capstyle='projecting',joinstyle='miter')
ctest.set_clip_on(False)
ax3.add_patch(ctest)

# import itertools
# ii = itertools.count(data.shape[0])
# clusters = [{'node_id': next(ii), 'left': x[0], 'right':x[1]} for x in model.children_]
#
# import copy
# n_points = data.shape[0]
# members = {i:[i] for i in range(n_points)}
# for cluster in clusters:
#     node_id = cluster["node_id"]
#     members[node_id] = copy.deepcopy(members[cluster["left"]])
#     members[node_id].extend(copy.deepcopy(members[cluster["right"]]))
#
# on_split = {c["node_id"]: [c["left"], c["right"]] for c in clusters}
# up_merge = {c["left"]: {"into": c["node_id"], "with": c["right"]} for c in clusters}
# # up_merge.update({c["right"]: {"into": c["node_id"], "with": c["left"]} for c in clusters})
# #
#
# from matplotlib.transforms import Bbox, TransformedBbox, \
#     blended_transform_factory
#
# from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
#     BboxConnectorPatch
#
#
# def connect_bbox(bbox1, bbox2,
#                  loc1a, loc2a, loc1b, loc2b,
#                  prop_lines, prop_patches=None):
#     if prop_patches is None:
#         prop_patches = prop_lines.copy()
#         prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2
#
#     c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
#     c1.set_clip_on(False)
#     c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
#     c2.set_clip_on(False)
#
#     bbox_patch1 = BboxPatch(bbox1, **prop_patches)
#     bbox_patch2 = BboxPatch(bbox2, **prop_patches)
#
#     p = BboxConnectorPatch(bbox1, bbox2,
#                            # loc1a=3, loc2a=2, loc1b=4, loc2b=1,
#                            loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
#                            **prop_patches)
#     p.set_clip_on(False)
#
#     return c1, c2, bbox_patch1, bbox_patch2, p
#
#
# def zoom_effect01(ax1, ax2, xmin, xmax, **kwargs):
#     """
#     ax1 : the main axes
#     ax1 : the zoomed axes
#     (xmin,xmax) : the limits of the colored area in both plot axes.
#
#     connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
#     be marked.  The keywords parameters will be used ti create
#     patches.
#
#     """
#
#     trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
#     trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)
#
#     bbox = Bbox.from_extents(xmin, 0, xmax, 1)
#
#     mybbox1 = TransformedBbox(bbox, trans1)
#     mybbox2 = TransformedBbox(bbox, trans2)
#
#     prop_patches = kwargs.copy()
#     prop_patches["ec"] = "none"
#     prop_patches["alpha"] = 0.2
#
#     c1, c2, bbox_patch1, bbox_patch2, p = \
#         connect_bbox(mybbox1, mybbox2,
#                      loc1a=3, loc2a=2, loc1b=4, loc2b=1,
#                      prop_lines=kwargs, prop_patches=prop_patches)
#
#     ax1.add_patch(bbox_patch1)
#     ax2.add_patch(bbox_patch2)
#     ax2.add_patch(c1)
#     ax2.add_patch(c2)
#     ax2.add_patch(p)
#
#     return c1, c2, bbox_patch1, bbox_patch2, p
#
#
# def zoom_effect02(ax1, ax2, **kwargs):
#     """
#     ax1 : the main axes
#     ax1 : the zoomed axes
#
#     Similar to zoom_effect01.  The xmin & xmax will be taken from the
#     ax1.viewLim.
#     """
#
#     tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
#     trans = blended_transform_factory(ax2.transData, tt)
#
#     mybbox1 = ax1.bbox
#
#     mybbox2 = TransformedBbox(ax1.viewLim, trans)
#
#     prop_patches = kwargs.copy()
#     prop_patches["ec"] = "none"
#     prop_patches["alpha"] = 0.2
#
#     c1, c2, bbox_patch1, bbox_patch2, p = \
#         connect_bbox(mybbox1, mybbox2,
#                      loc1a=3, loc2a=2, loc1b=4, loc2b=1,
#                      prop_lines=kwargs, prop_patches=prop_patches)
#
#     ctest = BboxConnector(mybbox1, mybbox2, loc1=3, loc2=2)
#     ctest.set_clip_on(False)
#     ax2.add_patch(ctest)
#     # ax1.add_patch(bbox_patch1)
#     # ax2.add_patch(bbox_patch2)
#     # ax2.add_patch(c1)
#     # ax2.add_patch(c2)
#     # ax2.add_patch(p)
#
#     return c1, c2, bbox_patch1, bbox_patch2, p, mybbox1, mybbox2
#
#
# import matplotlib.pyplot as plt
#
# plt.figure(1, figsize=(5, 5))
# ax1 = plt.subplot(221)
# ax2 = plt.subplot(212)
# ax2.set_xlim(0, 1)
# ax2.set_xlim(0, 5)
# zoom_effect01(ax1, ax2, 0.2, 0.8)
#
#
# ax1 = plt.subplot(222)
# ax1.set_xlim(2, 3)
# ax2.set_xlim(0, 5)
# c1, c2, bbox_patch1, bbox_patch2, p, mybbox1, mybbox2 = zoom_effect02(ax1, ax2)
#
# plt.show()
