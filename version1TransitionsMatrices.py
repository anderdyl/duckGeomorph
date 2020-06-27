
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

import pickle
dbfile = open('sandbarsSouthernTransect.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()
alllines = data['alllines']
xinterp = data['xinterp']
time = data['time']


# Ok, first lets load a DWT struct

import scipy.io as sio
from scipy.io.matlab.mio5_params import mat_struct
import numpy as np

def ReadMatfile(p_mfile):
    'Parse .mat file to nested python dictionaries'

    def RecursiveMatExplorer(mstruct_data):
        # Recursive function to extrat mat_struct nested contents

        if isinstance(mstruct_data, mat_struct):
            # mstruct_data is a matlab structure object, go deeper
            d_rc = {}
            for fn in mstruct_data._fieldnames:
                d_rc[fn] = RecursiveMatExplorer(getattr(mstruct_data, fn))
            return d_rc

        else:
            # mstruct_data is a numpy.ndarray, return value
            return mstruct_data

    # base matlab data will be in a dict
    mdata = sio.loadmat(p_mfile, squeeze_me=True, struct_as_record=False)
    mdata_keys = [x for x in mdata.keys() if x not in
                  ['__header__','__version__','__globals__']]

    # use recursive function
    dout = {}
    for k in mdata_keys:
        dout[k] = RecursiveMatExplorer(mdata[k])
    return dout



mat = scipy.io.loadmat('matlabsSouthernTransectClusters.mat')

numCenters = int(mat['clusters']['NumberCenters'][0][0].flatten())
group = mat['clusters']['Group'][0][0]
bmus = mat['clusters']['PatternsGroup'][0][0].flatten()
order = mat['order'].flatten()
#sorted_bmus = np.zeros((len(bmus),))
sorted_bmus = np.tile(0,(len(bmus),), )

#sorted_time = np.tile(0,(len(kma.labels_),), )
numClusters = numCenters

for i in range(numCenters):
    posc = np.where(bmus == order[i])
    sorted_bmus[posc] = int(i)
    #sorted_time[posc] = time[posc]


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




# def transitionDictionary(bmus,dates):

bins = dict()
date = dict()
nextbin = dict()
nextdate = dict()
prevbin = dict()
prevdate = dict()
for xx in range(numCenters):
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



m,m2 = transition_matrix(sorted_bmus)
for row in m: print(' '.join('{0:.2f}'.format(x) for x in row))

flat_list = [item for sublist in m for item in sublist]
flatarray = np.asarray(flat_list)
flatarray.resize(numClusters, numClusters)

flat_list2 = [item for sublist in m2 for item in sublist]
flatarray2 = np.asarray(flat_list2)
flatarray2.resize(numClusters, numClusters)

matrix = flatarray

import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, numClusters))
colors2 = cm.gray(np.linspace(0, 1, numClusters))


## Chord option #1

# chord diagram
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

import numpy as np

LW = 0.3

def polar2xy(r, theta):
    return np.array([r*np.cos(theta), r*np.sin(theta)])

def hex2rgb(c):
    return tuple(int(c[i:i+2], 16)/256.0 for i in (1, 3 ,5))

def IdeogramArc(start=0, end=60, radius=1.0, width=0.2, ax=None, color=(1,0,0)):
    # start, end should be in [0, 360)
    if start > end:
        start, end = end, start
    start *= np.pi/180.
    end *= np.pi/180.
    # optimal distance to the control points
    # https://stackoverflow.com/questions/1734745/how-to-create-circle-with-b%C3%A9zier-curves
    opt = 4./3. * np.tan((end-start)/ 4.) * radius
    inner = radius*(1-width)
    verts = [
        polar2xy(radius, start),
        polar2xy(radius, start) + polar2xy(opt, start+0.5*np.pi),
        polar2xy(radius, end) + polar2xy(opt, end-0.5*np.pi),
        polar2xy(radius, end),
        polar2xy(inner, end),
        polar2xy(inner, end) + polar2xy(opt*(1-width), end-0.5*np.pi),
        polar2xy(inner, start) + polar2xy(opt*(1-width), start+0.5*np.pi),
        polar2xy(inner, start),
        polar2xy(radius, start),
        ]

    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.LINETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CLOSEPOLY,
             ]

    if ax == None:
        return verts, codes
    else:
        path = Path(verts, codes)
        print(color)
        patch = patches.PathPatch(path, facecolor=color, edgecolor=color, lw=LW, alpha=0.5)
        #patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.4,), lw=LW)

        ax.add_patch(patch)


def ChordArc(start1=0, end1=60, start2=180, end2=240, radius=1.0, chordwidth=0.7, ax=None, color=(1,0,0)):
    # start, end should be in [0, 360)
    if start1 > end1:
        start1, end1 = end1, start1
    if start2 > end2:
        start2, end2 = end2, start2
    start1 *= np.pi/180.
    end1 *= np.pi/180.
    start2 *= np.pi/180.
    end2 *= np.pi/180.
    opt1 = 4./3. * np.tan((end1-start1)/ 4.) * radius
    opt2 = 4./3. * np.tan((end2-start2)/ 4.) * radius
    rchord = radius * (1-chordwidth)
    verts = [
        polar2xy(radius, start1),
        polar2xy(radius, start1) + polar2xy(opt1, start1+0.5*np.pi),
        polar2xy(radius, end1) + polar2xy(opt1, end1-0.5*np.pi),
        polar2xy(radius, end1),
        polar2xy(rchord, end1),
        polar2xy(rchord, start2),
        polar2xy(radius, start2),
        polar2xy(radius, start2) + polar2xy(opt2, start2+0.5*np.pi),
        polar2xy(radius, end2) + polar2xy(opt2, end2-0.5*np.pi),
        polar2xy(radius, end2),
        polar2xy(rchord, end2),
        polar2xy(rchord, start1),
        polar2xy(radius, start1),
        ]

    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             ]

    if ax == None:
        return verts, codes
    else:
        path = Path(verts, codes)
        #patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.4,), lw=LW)
        patch = patches.PathPatch(path, facecolor=color, edgecolor=color, lw=LW, alpha=0.5)

        ax.add_patch(patch)

def selfChordArc(start=0, end=60, radius=1.0, chordwidth=0.7, ax=None, color=(1,0,0)):
    # start, end should be in [0, 360)
    if start > end:
        start, end = end, start
    start *= np.pi/180.
    end *= np.pi/180.
    opt = 4./3. * np.tan((end-start)/ 4.) * radius
    rchord = radius * (1-chordwidth)
    verts = [
        polar2xy(radius, start),
        polar2xy(radius, start) + polar2xy(opt, start+0.5*np.pi),
        polar2xy(radius, end) + polar2xy(opt, end-0.5*np.pi),
        polar2xy(radius, end),
        polar2xy(rchord, end),
        polar2xy(rchord, start),
        polar2xy(radius, start),
        ]

    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             ]

    if ax == None:
        return verts, codes
    else:
        path = Path(verts, codes)
        #patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.4,), lw=LW)
        patch = patches.PathPatch(path, facecolor=color, edgecolor=color, lw=LW, alpha=0.5)

        ax.add_patch(patch)

def chordDiagram(X, ax, colors=None, width=0.1, pad=2, chordwidth=0.7):
    """Plot a chord diagram
    Parameters
    ----------
    X :
        flux data, X[i, j] is the flux from i to j
    ax :
        matplotlib `axes` to show the plot
    colors : optional
        user defined colors in rgb format. Use function hex2rgb() to convert hex color to rgb color. Default: d3.js category10
    width : optional
        width/thickness of the ideogram arc
    pad : optional
        gap pad between two neighboring ideogram arcs, unit: degree, default: 2 degree
    chordwidth : optional
        position of the control points for the chords, controlling the shape of the chords
    """
    # X[i, j]:  i -> j
    x = X.sum(axis = 1) # sum over rows
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)

    if colors is None:
    # use d3.js category10 https://github.com/d3/d3-3.x-api-reference/blob/master/Ordinal-Scales.md#category10
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                  '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        if len(x) > 10:
            print('x is too large! Use x smaller than 10')
        colors = [hex2rgb(colors[i]) for i in range(len(x))]

    # find position for each start and end
    y = x/np.sum(x).astype(float) * (360 - pad*len(x))

    pos = {}
    arc = []
    nodePos = []
    start = 0
    for i in range(len(x)):
        end = start + y[i]
        arc.append((start, end))
        angle = 0.5*(start+end)
        #print(start, end, angle)
        if -30 <= angle <= 210:
            angle -= 90
        else:
            angle -= 270
        nodePos.append(tuple(polar2xy(1.1, 0.5*(start+end)*np.pi/180.)) + (angle,))
        z = (X[i, :]/x[i].astype(float)) * (end - start)
        ids = np.argsort(z)
        z0 = start
        for j in ids:
            pos[(i, j)] = (z0, z0+z[j])
            z0 += z[j]
        start = end + pad

    for i in range(len(x)):
        start, end = arc[i]
        print(colors[i])
        IdeogramArc(start=start, end=end, radius=1.0, ax=ax, color=colors[i], width=width)
        start, end = pos[(i,i)]
        selfChordArc(start, end, radius=1.-width, color=colors[i], chordwidth=chordwidth*0.7, ax=ax)
        for j in range(i):
            color = colors[i]
            if X[i, j] > X[j, i]:
                color = colors[j]
            start1, end1 = pos[(i,j)]
            start2, end2 = pos[(j,i)]
            ChordArc(start1, end1, start2, end2,
                     radius=1.-width, color=colors[i], chordwidth=chordwidth, ax=ax)

    #print(nodePos)
    return nodePos

def matplotlib_to_plotly(cmap, pl_entries):
    h = 1.0/(pl_entries-1)
    pl_colorscale = []

    for k in range(pl_entries):
        #C = map(np.uint8, np.array(cmap(k*h)[:3])*255)
        C = np.uint8(np.array(cmap(k*h)[:3])*255)
        #pl_colorscale.append([k*h, 'rgb'+str((C[0], C[1], C[2]))])
        pl_colorscale.append('rgba'+str((C[0], C[1], C[2], 0.6)))
    return pl_colorscale


##################################
#if __name__ == "__main__":
fig = plt.figure(figsize=(6,6))
#flux = matrix
flux = matrix
ax = plt.axes([0,0,1,1])

    #nodePos = chordDiagram(flux, ax, colors=[hex2rgb(x) for x in ['#666666', '#66ff66', '#ff6666', '#6666ff']])
#nodePos = chordDiagram(flux, ax)
tempColorsInd = np.where((colors==1))
tempColorsInd2 = np.where((colors==0))

tempColors = colors
tempColors[tempColorsInd] = 0.9999
tempColors[tempColorsInd2] = 0.0001
tempColors = tempColors
nodePos = chordDiagram(flux,ax,colors=tempColors[:,0:3])
ax.axis('off')
prop = dict(fontsize=16*0.8, ha='center', va='center')
#plt.title('Winter')
plt.show()


onlyTransitions = flatarray
onlyTransitions[0][0] = 0
onlyTransitions[1][1] = 0
onlyTransitions[2][2] = 0
onlyTransitions[3][3] = 0
onlyTransitions[4][4] = 0
onlyTransitions[5][5] = 0
onlyTransitions[6][6] = 0
onlyTransitions[7][7] = 0
onlyTransitions[8][8] = 0
onlyTransitions[9][9] = 0
onlyTransitions[10][10] = 0
onlyTransitions[11][11] = 0
onlyTransitions[12][12] = 0
onlyTransitions[13][13] = 0
onlyTransitions[14][14] = 0

sources = list()
targets = list()
values = list()
for ff in range(numCenters):

    tempSource = ff*np.ones((numCenters,))
    tempTarget = np.arange(0,numCenters,)

    sources.append(tempSource)
    targets.append(tempTarget)
    values.append(onlyTransitions[ff])

sourcesSankey = np.asarray(sources).flatten()
valuesSankey = np.asarray(values).flatten()
targetsSankey = np.asarray(targets).flatten()

import matplotlib.colors
rainbow_cmap = matplotlib.cm.get_cmap('rainbow')
rainbow_rgb = []
norm = matplotlib.colors.Normalize(vmin=0, vmax=255)

for i in range(0, 255):
       k = matplotlib.colors.colorConverter.to_rgb(rainbow_cmap(norm(i)))
       rainbow_rgb.append(k)

colorsPyplot = matplotlib_to_plotly(rainbow_cmap,numCenters)
colorsNodes = matplotlib_to_plotly(rainbow_cmap,numCenters)
from itertools import repeat
for ff in range(numClusters-1):
    colorsPyplot.extend(colorsPyplot)

import plotly.graph_objects as go

fig = go.Figure(go.Sankey(
    arrangement = "snap",
    node = {
        "label": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"],
        "x": [0.1, 0.1, 0.1, 0.25, 0.25, 0.25, 0.4, 0.4, 0.4, 0.55, 0.55, 0.55, 0.7, 0.7, 0.7],
        "y": [0.1, 0.3, 0.5, 0.15, 0.35, 0.55, 0.2, 0.4, 0.6, 0.25, 0.45, 0.65, 0.3, 0.5, 0.7],
        "color": colorsNodes,
        'pad':10},  # 10 Pixels
    link = {
        "source": sourcesSankey,
        "target": targetsSankey,
        "value": valuesSankey,
        "color": colorsPyplot}))
fig.write_image("fig2.png")


#
# fig = go.Figure(go.Sankey(
#     arrangement = "snap",
#     node = {
#         "label": ["A", "B", "C", "D", "E", "F"],
#         "x": [0.2, 0.1, 0.5, 0.7, 0.3, 0.5],
#         "y": [0.7, 0.5, 0.2, 0.4, 0.2, 0.3],
#         'pad':10},  # 10 Pixels
#     link = {
#         "source": [0, 0, 1, 2, 5, 4, 3, 5],
#         "target": [5, 3, 4, 3, 0, 2, 2, 3],
#         "value": [1, 5, 1, 1, 1, 1, 1, 2]}))
#
# fig.write_image("fig2.png")
