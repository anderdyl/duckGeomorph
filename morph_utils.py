
import numpy as np
from sklearn.cluster import KMeans
import xarray as xr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from matplotlib import cm


def KMA_simple(xds_PCA, num_clusters, repres=0.95):
    '''
    KMeans Classification for PCA data

    xds_PCA:
        (n_components, n_components) PCs
        (n_components, n_features) EOFs
        (n_components, ) variance
    num_clusters
    repres

    returns a xarray.Dataset containing KMA data
    '''

    # PCA data
    variance = xds_PCA.variance.values[:]
    EOFs = xds_PCA.EOFs.values[:]
    PCs = xds_PCA.PCs.values[:]

    var_anom_std = xds_PCA.pred_std.values[:]
    var_anom_mean = xds_PCA.pred_mean.values[:]
    time = xds_PCA.pred_time.values[:]

    # APEV: the cummulative proportion of explained variance by ith PC
    APEV = np.cumsum(variance) / np.sum(variance)*100.0
    nterm = np.where(APEV <= repres*100)[0][-1]

    PCsub = PCs[:, :nterm+1]
    data = PCsub
    data_std = np.std(data, axis=0)
    data_mean = np.mean(data, axis=0)

    # normalize but keep PCs weigth
    data_norm = np.ones(data.shape)*np.nan
    for i in range(PCsub.shape[1]):
        data_norm[:,i] = np.divide(data[:,i]-data_mean[i], data_std[0])


    EOFsub = EOFs[:nterm+1, :]

    # KMEANS
    kma = KMeans(n_clusters=num_clusters, n_init=2000).fit(PCsub)

    # groupsize
    _, group_size = np.unique(kma.labels_, return_counts=True)

    # groups
    d_groups = {}
    for k in range(num_clusters):
        d_groups['{0}'.format(k)] = np.where(kma.labels_==k)
    # TODO: STORE GROUPS WITHIN OUTPUT DATASET

    # centroids
    centroids = np.dot(kma.cluster_centers_, EOFsub)
    # centroids
    # centroids = np.zeros((num_clusters, data.shape[1]))
    # for k in range(num_clusters):
    #     centroids[k,:] = np.mean(data[d_groups['{0}'.format(k)],:], axis=1)
    # km, x and var_centers
    km = np.multiply(
        centroids,
        np.tile(var_anom_std[0], (num_clusters, 1))
    ) + np.tile(var_anom_mean, (num_clusters, 1))

    # sort kmeans
    kma_order = np.argsort(np.mean(-km, axis=1))

    # reorder clusters: bmus, km, cenEOFs, centroids, group_size
    sorted_bmus = np.zeros((len(kma.labels_),),)*np.nan
    for i in range(num_clusters):
        posc = np.where(kma.labels_ == kma_order[i])
        sorted_bmus[posc] = i
    sorted_km = km[kma_order]
    sorted_cenEOFs = kma.cluster_centers_[kma_order]
    sorted_centroids = centroids[kma_order]
    sorted_group_size = group_size[kma_order]

    return xr.Dataset(
        {
            'bmus': (('n_pcacomp'), sorted_bmus.astype(int)),
            'cenEOFs': (('n_clusters', 'n_features'), sorted_cenEOFs),
            'centroids': (('n_clusters','n_pcafeat'), sorted_centroids),
            'Km': (('n_clusters','n_pcafeat'), sorted_km),
            'group_size': (('n_clusters'), sorted_group_size),

            # PCA data
            'PCs': (('n_pcacomp','n_features'), PCsub),
            'variance': (('n_pcafeat',), variance),
            'time': (('n_pcacomp',), time),
        }
    )


def sort_cluster_gen_corr_end(centers, dimdim):
    '''
    SOMs alternative
    '''
    # TODO: DOCUMENTAR.

    # get dimx, dimy
    dimy = np.floor(np.sqrt(dimdim)).astype(int)
    dimx = np.ceil(np.sqrt(dimdim)).astype(int)

    if not np.equal(dimx*dimy, dimdim):
        # TODO: RAISE ERROR
        pass

    dd = distance_matrix(centers, centers)
    qx = 0
    sc = np.random.permutation(dimdim).reshape(dimy, dimx)

    # get qx
    for i in range(dimy):
        for j in range(dimx):

            # row F-1
            if not i==0:
                qx += dd[sc[i-1,j], sc[i,j]]

                if not j==0:
                    qx += dd[sc[i-1,j-1], sc[i,j]]

                if not j+1==dimx:
                    qx += dd[sc[i-1,j+1], sc[i,j]]

            # row F
            if not j==0:
                qx += dd[sc[i,j-1], sc[i,j]]

            if not j+1==dimx:
                qx += dd[sc[i,j+1], sc[i,j]]

            # row F+1
            if not i+1==dimy:
                qx += dd[sc[i+1,j], sc[i,j]]

                if not j==0:
                    qx += dd[sc[i+1,j-1], sc[i,j]]

                if not j+1==dimx:
                    qx += dd[sc[i+1,j+1], sc[i,j]]

    # test permutations
    q=np.inf
    go_out = False
    for i in range(dimdim):
        if go_out:
            break

        go_out = True

        for j in range(dimdim):
            for k in range(dimdim):
                if len(np.unique([i,j,k]))==3:

                    u = sc.flatten('F')
                    u[i] = sc.flatten('F')[j]
                    u[j] = sc.flatten('F')[k]
                    u[k] = sc.flatten('F')[i]
                    u = u.reshape(dimy, dimx, order='F')

                    f=0
                    for ix in range(dimy):
                        for jx in range(dimx):

                            # row F-1
                            if not ix==0:
                                f += dd[u[ix-1,jx], u[ix,jx]]

                                if not jx==0:
                                    f += dd[u[ix-1,jx-1], u[ix,jx]]

                                if not jx+1==dimx:
                                    f += dd[u[ix-1,jx+1], u[ix,jx]]

                            # row F
                            if not jx==0:
                                f += dd[u[ix,jx-1], u[ix,jx]]

                            if not jx+1==dimx:
                                f += dd[u[ix,jx+1], u[ix,jx]]

                            # row F+1
                            if not ix+1==dimy:
                                f += dd[u[ix+1,jx], u[ix,jx]]

                                if not jx==0:
                                    f += dd[u[ix+1,jx-1], u[ix,jx]]

                                if not jx+1==dimx:
                                    f += dd[u[ix+1,jx+1], u[ix,jx]]

                    if f<=q:
                        q = f
                        sc = u

                        if q<=qx:
                            qx=q
                            go_out=False

    return sc.flatten('F')



def Plot_DWTs_Mean_Anom(xds_KMA, xds_var, kind='mean', mask_land=None, p_export=None):
    '''
    Plot Daily Weather Types (bmus mean)
    kind - mean/anom
    '''

    bmus = xds_KMA['sorted_bmus'].values[:]
    n_clusters = len(xds_KMA.n_clusters.values[:])

    var_max = np.max(xds_var.values)
    var_min = np.min(xds_var.values)
    scale = 1/100.0  # scale from Pa to mbar

    # get number of rows and cols for gridplot
    n_rows, n_cols = GetBestRowsCols(n_clusters)

    # get cluster colors
    cs_dwt = colors_dwt(n_clusters)

    # plot figure
    #fig = plt.figure(figsize=(_faspect*_fsize, _fsize))
    fig = plt.figure(figsize=(10, 10))

    gs = gridspec.GridSpec(n_rows, n_cols, wspace=0.0, hspace=0.0)
    gr, gc = 0, 0

    for ic in range(n_clusters):

        if kind=='mean':
            # data mean
            it = np.where(bmus==ic)[0][:]
            c_mean = xds_var.isel(time=it).mean(dim='time')
            c_plot = np.multiply(c_mean, scale)  # apply scale

        elif kind=='anom':
            # data anomally
            it = np.where(bmus==ic)[0][:]
            t_mean = xds_var.mean(dim='time')
            c_mean = xds_var.isel(time=it).mean(dim='time')
            c_anom = c_mean - t_mean
            c_plot = np.multiply(c_anom, scale)  # apply scale

        # dwt color
        clr = cs_dwt[ic]

        # axes plot
        ax = plt.subplot(gs[gr, gc])
        pc = axplot_DWT(
            ax, np.flipud(c_plot),
            vmin = var_min, vmax = var_max,
            wt_num = ic+1,
            land = mask_land, wt_color = clr,
        )

        # anomalies colorbar center at 0
        if kind == 'anom':
            pc.set_clim(-6,6)

        # get lower positions
        if gr==n_rows-1 and gc==0:
            pax_l = ax.get_position()
        elif gr==n_rows-1 and gc==n_cols-1:
            pax_r = ax.get_position()

        # counter
        gc += 1
        if gc >= n_cols:
            gc = 0
            gr += 1

    # add a colorbar
    cbar_ax = fig.add_axes([pax_l.x0, pax_l.y0-0.05, pax_r.x1 - pax_l.x0, 0.02])
    cb = fig.colorbar(pc, cax=cbar_ax, orientation='horizontal')
    if kind=='mean':
        cb.set_label('Pressure (mbar)')
    elif kind=='anom':
        cb.set_label('Pressure anomalies (mbar)')

    # show / export
#    if not p_export
#     plt.show()
#    else:
#        fig.savefig(p_export, dpi=_fdpi)
#        plt.close()




def axplot_DWT(ax, dwt, vmin, vmax, wt_num, land=None, wt_color=None):
    'axes plot EOFs 2d map'

    cmap = copy.deepcopy(cm.get_cmap('RdBu_r'))

    # EOF pcolormesh
    pc = ax.pcolormesh(
        dwt,
        cmap = cmap, shading = 'gouraud',
        clim = (vmin, vmax),
    )

    # optional mask land
    if type(land).__module__ == np.__name__:
        landc = land.copy()
        landc[np.isnan(land)]=1
        landc[land==1]=np.nan
        ax.pcolormesh(
            np.flipud(landc),
            cmap=colors.ListedColormap([wt_color]), shading='gouraud',
        )

    # wt text
    ax.text(0.87, 0.85, wt_num, transform=ax.transAxes, fontweight='bold')

    # customize axis
    ax.set_xticks([])
    ax.set_yticks([])

    return pc




def GetBestRowsCols(n):
    'try to square number n, used at gridspec plots'
    from math import sqrt

    sqrt_n = sqrt(n)
    if sqrt_n.is_integer():
        n_r = int(sqrt_n)
        n_c = int(sqrt_n)
    else:
        l_div = GetDivisors(n)
        n_c = l_div[int(len(l_div)/2)]
        n_r = int(n/n_c)

    return n_r, n_c

def GetDivisors(x):
    l_div = []
    i = 1
    while i<x:
        if x%i == 0:
            l_div.append(i)
        i = i + 1
    return l_div




def colors_dwt(num_clusters):

    # 42 DWT colors
    l_colors_dwt = [
        (1.0000, 0.1344, 0.0021),
        (1.0000, 0.2669, 0.0022),
        (1.0000, 0.5317, 0.0024),
        (1.0000, 0.6641, 0.0025),
        (1.0000, 0.9287, 0.0028),
        (0.9430, 1.0000, 0.0029),
        (0.6785, 1.0000, 0.0031),
        (0.5463, 1.0000, 0.0032),
        (0.2821, 1.0000, 0.0035),
        (0.1500, 1.0000, 0.0036),
        (0.0038, 1.0000, 0.1217),
        (0.0039, 1.0000, 0.2539),
        (0.0039, 1.0000, 0.4901),
        (0.0039, 1.0000, 0.6082),
        (0.0039, 1.0000, 0.8444),
        (0.0039, 1.0000, 0.9625),
        (0.0039, 0.8052, 1.0000),
        (0.0039, 0.6872, 1.0000),
        (0.0040, 0.4510, 1.0000),
        (0.0040, 0.3329, 1.0000),
        (0.0040, 0.0967, 1.0000),
        (0.1474, 0.0040, 1.0000),
        (0.2655, 0.0040, 1.0000),
        (0.5017, 0.0040, 1.0000),
        (0.6198, 0.0040, 1.0000),
        (0.7965, 0.0040, 1.0000),
        (0.8848, 0.0040, 1.0000),
        (1.0000, 0.0040, 0.9424),
        (1.0000, 0.0040, 0.8541),
        (1.0000, 0.0040, 0.6774),
        (1.0000, 0.0040, 0.5890),
        (1.0000, 0.0040, 0.4124),
        (1.0000, 0.0040, 0.3240),
        (1.0000, 0.0040, 0.1473),
        (0.9190, 0.1564, 0.2476),
        (0.7529, 0.3782, 0.4051),
        (0.6699, 0.4477, 0.4584),
        (0.5200, 0.5200, 0.5200),
        (0.4595, 0.4595, 0.4595),
        (0.4100, 0.4100, 0.4100),
        (0.3706, 0.3706, 0.3706),
        (0.2000, 0.2000, 0.2000),
        (     0, 0, 0),
    ]

    # get first N colors
    np_colors_base = np.array(l_colors_dwt)
    np_colors_rgb = np_colors_base[:num_clusters]

    return np_colors_rgb

