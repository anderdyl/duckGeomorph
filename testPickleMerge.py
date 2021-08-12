# import pickle
# pkl1 = open('testPickles/2010_JulAug_L5_UPPER_output.pkl', 'rb')
# data = pickle.load(pkl1)
# pkl1.close()
#
# pkl2 = open('testPickles/2010_SepOct_L5_UPPER_output.pkl', 'rb')
# data2 = pickle.load(pkl2)
# pkl2.close()
#
# dates1 = data['dates']
# dates2 = data2['dates']
# shorelines1 = data['shorelines']
# shorelines2 = data2['shorelines']
#
# mergedDates = []
# mergedDates.extend(dates1)
# mergedDates.extend(dates2)
#
# mergedShorelines = []
# mergedShorelines.extend(shorelines1)
# mergedShorelines.extend(shorelines2)
#
# mergedDict = dict()
# mergedDict['dates'] = mergedDates
# mergedDict['shorelines'] = mergedShorelines

import pickle
import os

mergedPickle = 'mergedPickle.pickle'


dir = '/home/dylananderson/projects/duckGeomorph/testPickles'
# Need to sort the files to ensure correct temporal order...
files = os.listdir(dir)
files.sort()
files_path = [os.path.join(os.path.abspath(dir), x) for x in files]

mergedData = dict()
mergedData['dates'] = []
mergedData['shorelines'] = []
mergedData['filename'] = []
mergedData['cloud_cover'] = []

for i in files_path:
    # function in here that appends the latest dictionary
    pkl = open(i, 'rb')
    data = pickle.load(pkl)
    pkl.close()
    datesNew = data['dates']
    shorelinesNew = data['shorelines']
    mergedData['dates'].extend(datesNew)
    mergedData['shorelines'].extend(shorelinesNew)


with open(mergedPickle,'wb') as f:
    pickle.dump(mergedData, f)
