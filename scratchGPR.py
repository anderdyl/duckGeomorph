

import pickle
import numpy as np
from matplotlib import pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, Matern, WhiteKernel
import random

pickleName = 'eofsForHypotheticalProfiles23Sept2020.pickle'
dbfile = open(pickleName, 'rb')
eofs = pickle.load(dbfile)
dbfile.close()

pickleName = 'hypotheticalStorms.pickle'
dbfile = open(pickleName, 'rb')
storms = pickle.load(dbfile)
dbfile.close()
predictors = np.vstack((storms['hypotheticalStorms'][0:14,:],storms['hypotheticalStorms'][15:43,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][44:45,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][46:69,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][70:97,:]))
predictors = np.vstack((predictors,storms['hypotheticalStorms'][98:100,:]))

predictors = np.hstack((predictors,eofs['oldEOFscores']))
predictands = eofs['newEOFscores']

output = dict()
output['predictors'] = predictors
output['predictands'] = predictands
import scipy.io
scipy.io.savemat('profileTrainingTemp.mat',output)

randomIndices = np.arange(0, 95)
np.random.shuffle(randomIndices)

tot = 95
sub = 19

plt.figure()
ax1 = plt.subplot2grid((1, 5), (0, 0), rowspan=1, colspan=1)
ax2 = plt.subplot2grid((1, 5), (0, 1), rowspan=1, colspan=1)
ax3 = plt.subplot2grid((1, 5), (0, 2), rowspan=1, colspan=1)
ax4 = plt.subplot2grid((1, 5), (0, 3), rowspan=1, colspan=1)
ax5 = plt.subplot2grid((1, 5), (0, 4), rowspan=1, colspan=1)
counter = 0
for xx in range(1):

    valIndices = randomIndices[counter:counter+sub]

    select = np.in1d(range(tot), valIndices)

    mode = 0
    trainX1 = predictors[select,:]
    trainY1 = predictands[select,mode]
    valX1 = predictors[~select,:]
    valY1 = predictands[~select,mode]
    mode = 1
    trainX2 = predictors[select,:]
    trainY2 = predictands[select,mode]
    valX2 = predictors[~select,:]
    valY2 = predictands[~select,mode]
    mode = 2
    trainX3 = predictors[select,:]
    trainY3 = predictands[select,mode]
    valX3 = predictors[~select,:]
    valY3 = predictands[~select,mode]
    mode = 3
    trainX4 = predictors[select,:]
    trainY4 = predictands[select,mode]
    valX4 = predictors[~select,:]
    valY4 = predictands[~select,mode]
    mode = 4
    trainX5 = predictors[select,:]
    trainY5 = predictands[select,mode]
    valX5 = predictors[~select,:]
    valY5 = predictands[~select,mode]

    # Instantiate a Gaussian Process model
    kernel1 = 300 * RBF(length_scale=10,length_scale_bounds=(1e-2, 1e4)) + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1e+1))
    kernel2 = 200 * RBF(length_scale=10,length_scale_bounds=(1e1, 1e4)) + WhiteKernel(noise_level=1e-1, noise_level_bounds=(1e-5, 1e+2))
    kernel3 = 200 * RBF(length_scale=100,length_scale_bounds=(1e1, 1e5)) + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1e+2))
    kernel4 = 200 * RBF(length_scale=10,length_scale_bounds=(1e-1, 1e4)) + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1e+1))
    kernel5 = 200 * RBF(length_scale=100,length_scale_bounds=(1e-1, 1e4)) + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1e+1))
    # kernel1 = 100.0 * Matern(length_scale=200, nu=15)#0.1 * RBF(length_scale=1,length_scale_bounds=(1e-3, 1e2))# + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1e+1))
    # kernel2 = 100.0 * Matern(length_scale=100, nu=1.5)#0.2 * RBF(length_scale=1,length_scale_bounds=(1e-3, 1e2))# + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1e+1))
    # kernel3 = 10.0 * Matern(length_scale=1, nu=1.5)#0.3 * RBF(length_scale=1,length_scale_bounds=(1e-3, 1e2))# + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1e+1))
    # kernel4 = 10.0 * Matern(length_scale=1, nu=1.5)#0.4 * RBF(length_scale=1,length_scale_bounds=(1e-3, 1e2))# + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1e+1))
    # kernel5 = 10.0 * Matern(length_scale=1, nu=1.5)#0.5 * RBF(length_scale=1,length_scale_bounds=(1e-3, 1e2))# + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1e+1))

    # kernel = 1.0 * Matern(length_scale=1.0, nu=1.5)#
    #
    #
    gp1 = GaussianProcessRegressor(kernel=kernel1, optimizer = 'fmin_l_bfgs_b', n_restarts_optimizer=50,normalize_y=False,random_state=5,alpha=.1)
    gp2 = GaussianProcessRegressor(kernel=kernel2, optimizer = 'fmin_l_bfgs_b', n_restarts_optimizer=50,normalize_y=False,random_state=5,alpha=.3)
    gp3 = GaussianProcessRegressor(kernel=kernel3, optimizer = 'fmin_l_bfgs_b', n_restarts_optimizer=50,normalize_y=False,random_state=5,alpha=.3)
    gp4 = GaussianProcessRegressor(kernel=kernel4, optimizer = 'fmin_l_bfgs_b', n_restarts_optimizer=50,normalize_y=False,random_state=5,alpha=.1)
    gp5 = GaussianProcessRegressor(kernel=kernel5, optimizer = 'fmin_l_bfgs_b', n_restarts_optimizer=50,normalize_y=False,random_state=5,alpha=.1)
    # gp1 = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=30,alpha=.1,normalize_y=True,random_state=1)
    # gp2 = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=30,alpha=.1,normalize_y=True,random_state=1)
    # gp3 = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=30,alpha=.3,normalize_y=True,random_state=1)
    # gp4 = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=30,alpha=.1,normalize_y=True,random_state=1)
    # gp5 = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=30,alpha=.3,normalize_y=True,random_state=1)

    # Fit to data using Maximum Likelihood Estimation of the parameters
    gp1.fit(trainX1, trainY1)
    gp2.fit(trainX2, trainY2)
    gp3.fit(trainX3, trainY3)
    gp4.fit(trainX4, trainY4)
    gp5.fit(trainX5, trainY5)

    y_pred1, sigma1 = gp1.predict(valX1, return_std=True)
    y_pred2, sigma2 = gp2.predict(valX2, return_std=True)
    y_pred3, sigma3 = gp3.predict(valX3, return_std=True)
    y_pred4, sigma4 = gp4.predict(valX4, return_std=True)
    y_pred5, sigma5 = gp5.predict(valX5, return_std=True)

    ax1.plot(valY1,y_pred1,'.')
    ax1.set_title('EOF 1')
    ax1.plot([-12,12],[-12,12],'k--')
    ax1.set_xlim([-12, 12])
    ax1.set_ylabel('GPR PCs')
    ax1.set_xlabel('Post-Storm PCs')

    ax2.plot(valY2,y_pred2,'.')
    ax2.set_title('EOF 2')
    ax2.plot([-9,9],[-9,9],'k--')
    ax2.set_xlim([-9, 9])
    ax2.set_xlabel('Post-Storm PCs')

    ax3.plot(valY3,y_pred3,'.')
    ax3.set_title('EOF 3')
    ax3.plot([-6,6],[-6,6],'k--')
    ax3.set_xlim([-6, 6])
    ax3.set_xlabel('Post-Storm PCs')

    ax4.plot(valY4,y_pred4,'.')
    ax4.set_title('EOF 4')
    ax4.plot([-6,6],[-6,6],'k--')
    ax4.set_xlim([-6, 6])
    ax4.set_xlabel('Post-Storm PCs')

    ax5.plot(valY5,y_pred5,'.')
    ax5.set_title('EOF 5')
    ax5.plot([-4.5,4.5],[-4.5,4.5],'k--')
    ax5.set_xlim([-4.5, 4.5])
    ax5.set_xlabel('Post-Storm PCs')

    counter = counter + sub





