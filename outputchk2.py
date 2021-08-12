# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 06:00:23 2021

@author: agharag
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pickle as pkl
import os

def onclick(event):
    if event.button==3:
        ix, iy = round(event.xdata,2), round(event.ydata,2)
#        print ('x = {}, y = {}'.format(ix,iy))
        coords.append((ix, iy))
    elif event.button==1:
        print("use rightclick to pick points")


def fit(dist_arr,Z_arr):
    def Equilibrium(x, a, b):
        return a*np.power(x, b)


    pars, cov = curve_fit(f=Equilibrium, xdata=dist_arr, ydata=Z_arr, p0=[0, 0], bounds=(-np.inf, np.inf))
    
    stdevs = np.sqrt(np.diag(cov))
    res = Z_arr - Equilibrium(dist_arr, *pars)
    return pars, stdevs,res


# def makefit(R1,R2):
#     ind_s=np.where(R>=R1)[0][0]
#     ind_e=np.where(R>=R2)[0][0]
#     pars, stdevs,res=fit(R[ind_s:ind_e+1]-R[ind_s],np.insert(post[i][j][ind_s+1:ind_e+1],0,pre[i][j][ind_s])-pre[i][j][ind_s])
#     equZ= pars[0]*np.power(R[ind_s:ind_e+1]-R[ind_s], pars[1])+pre[i][j][ind_s]
#     return equZ, pars, ind_s, ind_e
def makefit(R1,R2):
    ind_s=np.where(X>=R1)[0][0]
    ind_e=np.where(X>=R2)[0][0]
    pars, stdevs,res=fit(X[ind_s:ind_e+1]-X[ind_s],np.insert(newPost[ind_s+1:ind_e+1],0,newPre[ind_s])-newPre[ind_s])
    equZ= pars[0]*np.power(X[ind_s:ind_e+1]-X[ind_s], pars[1])+newPre[ind_s]
    return equZ, pars, ind_s, ind_e



geomorphdir = '/home/dylananderson/projects/duckGeomorph/scarp_data_full_DAversion/'
files = os.listdir(geomorphdir)
files.sort()
files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]



for i in range(25):
    dbfile = open(os.path.join(files_path[i+3], files[i+3]), 'rb')
    scarpPointsTemp = pkl.load(dbfile)
    dbfile.close()
    if i == 0:
        scarpPointsMid = scarpPointsTemp[:,:]
    else:
        scarpPointsMid = np.vstack((scarpPointsMid, scarpPointsTemp[:,:]))


scarpPointsMid = scarpPointsMid[:,0,:]


with open(r'/home/dylananderson/projects/duckGeomorph/scarp_data_full_DAversion/newcase_profiles_0000_0500_revised.pkl', "rb") as input_file:
    outputProfiles = pkl.load(input_file)

with open(r'/home/dylananderson/projects/duckGeomorph/scarp_data_full_DAversion/newcase_profiles_0500_1000.pkl', "rb") as input_file:
    outputProfiles2 = pkl.load(input_file)

with open(r'/home/dylananderson/projects/duckGeomorph/scarp_data_full_DAversion/newcase_profiles_1000_1250.pkl', "rb") as input_file:
    outputProfiles3 = pkl.load(input_file)

preProfiles1 = np.vstack((outputProfiles['pre'],outputProfiles2['pre']))
postProfiles1 = np.vstack((outputProfiles['post'],outputProfiles2['post']))
preProfiles = np.vstack((preProfiles1,outputProfiles3['pre']))
postProfiles = np.vstack((postProfiles1,outputProfiles3['post']))
dist = outputProfiles['dist']
pre = preProfiles
post = postProfiles
R = dist
X = np.arange(0,np.max(R),0.5)

#
# # data1=pkl.load(open('case_0_499_profile_result.pickle','rb'))
# # R=data1['dist']
# # pre=data1['pre']
# # post=data1['post']
# with open(r'/home/dylananderson/projects/duckGeomorph/scarp_data_final/newcase_profiles_0_500.pkl', "rb") as input_file:
#     outputProfiles = pkl.load(input_file)
#
# with open(r'/home/dylananderson/projects/duckGeomorph/scarp_data_final/newcase_profiles_500_1000.pkl', "rb") as input_file:
#     outputProfiles2 = pkl.load(input_file)
#
# preProfiles = np.vstack((outputProfiles['pre'],outputProfiles2['pre']))
# postProfiles = np.vstack((outputProfiles['post'],outputProfiles2['post']))
# dist = outputProfiles['dist']
# R = dist
# pre = preProfiles
# post = postProfiles
#
#
# X = np.arange(0,np.max(R),0.5)
#
# import os
# geomorphdir = '/home/dylananderson/projects/duckGeomorph/scarp_data_final/'
# files = os.listdir(geomorphdir)
# files.sort()
# files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]
# for i in range(6):
#     dbfile = open(os.path.join(files_path[i+2], files[i+2]), 'rb')
#     scarpPointsTemp = pkl.load(dbfile)
#     dbfile.close()
#     if i == 0:
#         scarpPointsMid = scarpPointsTemp[:,:]
#     else:
#         scarpPointsMid = np.vstack((scarpPointsMid, scarpPointsTemp[:,:]))


#
# Data=[]
# for i in range(len(files)):
#     dbfile = open(os.path.join(files_path[i], files[i]), 'rb')
#     data2 = pkl.load(dbfile)
#     dbfile.close()
#     Data.append(data2)

#    if i == 0:
#        scarpPointsMid = scarpPointsTemp[:,1,:]
#        scarpPointsLeft = scarpPointsTemp[:,0,:]
#        scarpPointsRight = scarpPointsTemp[:,2,:]
#    else:
#        scarpPointsMid = np.vstack((scarpPointsMid, scarpPointsTemp[:,1,:]))
#        scarpPointsLeft = np.vstack((scarpPointsLeft, scarpPointsTemp[:,0,:]))
#        scarpPointsRight = np.vstack((scarpPointsRight, scarpPointsTemp[:,2,:]))



#Data=[]
#for i in range (5):
#    filename='scarp_{}_{}.pkl'.format(i*100, i*100 +99)
#    data2=pkl.load(open(filename,'rb'))
#    Data.append(data2)
    
Data=np.array(scarpPointsMid)
Data2 = Data#.reshape(1000, 3, 6)


fig1=plt.figure()

Coord=[]

for ii in range(1):
    i=ii+802 ##change to i=ii+(start)
    temp=Data2[i,:]
    Arrj=[]
    pcoords=[]
    #for j in range(3):
    j = 1
    coords=[]
    plt.cla()
    plt.xlim(100,300)
    plt.ylim(-5,6)
    plt.title('case {}'.format(i))
    from scipy.interpolate import interp1d
    postFit = interp1d(R,post[i,:],kind='cubic')
    preFit = interp1d(R,pre[i,:],kind='cubic')
    newPost = postFit(X)
    newPre = preFit(X)
    # plt.plot(R,post[i][j])
    # plt.plot(R,pre[i][j],'--')
    plt.plot(X,newPost)
    plt.plot(X,newPre,'--')
    if np.all(np.isnan(temp)):# or (i==180 and j==1):
        ind_s=0
        ind_e=1
        equZ=[0,0]
        pars=[0,0]
            
    else:
        #equZ,pars,ind_s,ind_e=makefit(temp[j][0],temp[j][2])
        equZ,pars,ind_s,ind_e=makefit(temp[0],temp[4])

    ln,=plt.plot(X[ind_s:ind_e+1],equZ,'--',color='hotpink',linewidth=2)

    # ln,=plt.plot(R[ind_s:ind_e+1],equZ,'--',color='hotpink',linewidth=2)
    plt.plot(temp[0],temp[1],'o',color='magenta')
    plt.plot(temp[2],temp[3],'o',color='cyan')
                
    plt.pause(0.25)
    cid = fig1.canvas.mpl_connect('button_press_event', onclick)
                
    while not plt.waitforbuttonpress():
        try:
            if len(coords)>=2:

                equZ,pars,ind_s,ind_e=makefit(coords[0][0],coords[-1][0])
                #print('equZ = {}'.format(equZ))
                print('pars = {}'.format(pars))
                print('X = {}'.format(X[ind_s]))
                ln.remove()
                # ln,=plt.plot(R[ind_s:ind_e+1],equZ,'--',color='hotpink',linewidth=2)
                ln,=plt.plot(X[ind_s:ind_e+1],equZ,'--',color='hotpink',linewidth=2)

                plt.show()
        except  TypeError:
            coords=[]
            print('reset :: pick 1st and 2nd points')
        except  ValueError:
            coords=[]
            print('reset :: pick 1st and 2nd points')

        pass
                
        # alldata=[R[ind_s],pre[i][j][ind_s],pars[0],pars[1],R[ind_e],post[i][j][ind_e]]
    alldata=[X[ind_s],newPre[ind_s],pars[0],pars[1],X[ind_e],newPost[ind_e]]

    if alldata[0]<1:
        alldata=[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
    pcoords.append(alldata)

    Coord.append(pcoords) #np.round(pcoords,3))
    pklname='scarp_data_{}_{}_DAversionFinal.pkl'.format(i-ii,i)

Coord=np.array(Coord)


# with open(pklname,'wb') as f:
#     pkl.dump(Coord, f)

