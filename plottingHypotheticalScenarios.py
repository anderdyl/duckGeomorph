

import numpy as np
import pickle
import matplotlib.pyplot as plt


# Python function to print common elements in three sorted arrays
def findCommon(ar1, ar2, ar3, n1, n2, n3):
    # Initialize starting indexes for ar1[], ar2[] and ar3[]
    i, j, k = 0, 0, 0

    # Iterate through three arrays while all arrays have elements
    while (i < n1 and j < n2 and k < n3):

        # If x = y and y = z, print any of them and move ahead
        # in all arrays
        if (ar1[i] == ar2[j] and ar2[j] == ar3[k]):
            print(ar1[i])
            output = ar1[i]
            i += 1
            j += 1
            k += 1

        # x < y
        elif ar1[i] < ar2[j]:
            i += 1

        # y < z
        elif ar2[j] < ar3[k]:
            j += 1

        # We reach here when x > y and z < y, i.e., z is smallest
        else:
            k += 1
    return output


with open(r'hypoNourishStormTides2.pickle', "rb") as input_file:
    outputMDA = pickle.load(input_file)

mdaSub = outputMDA['hypo']

outHypo = dict()
outHypo['hypo'] = mdaSub
import scipy.io
scipy.io.savemat('hypoFinal500.mat',outHypo)




with open(r'hypotheticalTrainingConditions2.pickle', "rb") as input_file:
    outputHypo = pickle.load(input_file)

PCs1 = outputHypo['PCs1']
PCs2 = outputHypo['PCs2']
PCs3 = outputHypo['PCs3']


with open(r'nourishmentProfileEOFs2.pickle', "rb") as input_file:
    outputEOFs = pickle.load(input_file)

meanZ = outputEOFs['meanZ']
EOFs = outputEOFs['EOFs']
x = outputEOFs['x']




plt.figure()
nr = 0
nc = 0
for num in range(100):

        ax = plt.subplot2grid((5,5),(nr,nc),rowspan=1,colspan=1)
        finder = np.where(((mdaSub[num,0]-.001)<PCs1) & ((mdaSub[num,0]+.001)>PCs1))
        finder2 = np.where(((mdaSub[num,1]-.001)<PCs2) & ((mdaSub[num,1]+.001)>PCs2))
        finder3 = np.where(((mdaSub[num,2]-.001)<PCs3) & ((mdaSub[num,2]+.001)>PCs3))

        ztemp = meanZ + EOFs[0, :]*mdaSub[num,0] + EOFs[1, :]*mdaSub[num,1] + EOFs[2, :] * mdaSub[num,2] + EOFs[3, :] * \
                mdaSub[num,3] + EOFs[4, :]*mdaSub[num,4] + EOFs[5, :]*mdaSub[num,5] + EOFs[6, :]*mdaSub[num,6] + EOFs[7,:]*mdaSub[num,7] + EOFs[8,:]*mdaSub[num,8]

        output = findCommon(finder[0],finder2[0],finder3[0],len(finder[0]),len(finder2[0]),len(finder3[0]))
        ax.plot(x, ztemp, '-',label='Xbeach Profile')
        ax.plot((x[0],x[-1]),(0,0),'.-',label='zero water level')
        ax.plot((x[0],x[-1]),(mdaSub[num,9],mdaSub[num,9]),label='MSL')
        # ax.plot(xShort, zInt[output, :],'--')
        # ax.set_title('profile{}'.format(num))
        ax.text(250,2,'profile{}'.format(num))
        ax.text(250,-1,'MSL={}m'.format(np.round(mdaSub[num,9]*1000)/1000))
        ax.text(250,-2,'HS={}m'.format(np.round(mdaSub[num,10]*10)/10))
        ax.text(250,-3,'Surge={}m'.format(np.round(mdaSub[num,12]*10)/10))

        ax.set_xlim([0, 400])
        ax.set_ylim([-4, 4])
        if nc < 4:
            nc = nc+1
        else:
            nc = 0
            nr = nr+1
        # plt.savefig('/home/dylananderson/projects/duckGeomorph/hypoNourishments/profile{}.png'.format(num))
        # plt.close()




plt.figure()
plt.hist(mdaSub[0:1000,9],20)
plt.xlabel('MSL (m)')
plt.title('First 1000 scenarios')



plt.figure()
plt.subplot2grid((1,4),(0,0),rowspan=1,colspan=1)
plt.hist(mdaSub[0:1000,0],20)
plt.xlabel('EOF #1')
plt.title('First 1000 scenarios')
plt.subplot2grid((1,4),(0,1),rowspan=1,colspan=1)
plt.hist(mdaSub[0:1000,1],20)
plt.xlabel('EOF #2')
plt.title('First 1000 scenarios')
plt.subplot2grid((1,4),(0,2),rowspan=1,colspan=1)
plt.hist(mdaSub[0:1000,2],20)
plt.xlabel('EOF #3')
plt.title('First 1000 scenarios')
plt.subplot2grid((1,4),(0,3),rowspan=1,colspan=1)
plt.hist(mdaSub[0:1000,3],20)
plt.xlabel('EOF #4')
plt.title('First 1000 scenarios')


