



import pickle

with open(r'hypotheticalTrainingConditions2.pickle', "rb") as input_file:
    outputData = pickle.load(input_file)

dataN = outputData['dataN']
dataT = outputData['dataT']
dataS = outputData['dataS']

M2amp = outputData['M2amp']
S2amp = outputData['S2amp']
K1amp = outputData['K1amp']
N2amp = outputData['N2amp']

M2period = outputData['M2period']
N2period = outputData['N2period']
K1period = outputData['K1period']
S2period = outputData['S2period']

M2speed = outputData['M2speed']
N2speed = outputData['N2speed']
K1speed = outputData['K1speed']
S2speed = outputData['S2speed']

PCs1 = outputData['PCs1']
PCs2 = outputData['PCs2']
PCs3 = outputData['PCs3']
PCs4 = outputData['PCs4']
PCs5 = outputData['PCs5']
PCs6 = outputData['PCs6']
PCs7 = outputData['PCs7']
PCs8 = outputData['PCs8']
PCs9 = outputData['PCs9']






