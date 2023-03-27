import adios2
import sys
import numpy as np

if len(sys.argv) != 3 :
    print('usage: ', sys.argv[0], ' infile outfile')
    sys.exit()

inFile = sys.argv[1]
out = sys.argv[2]

outName = out + '.txt'

##Read in data.
f = adios2.open(inFile, 'r')
ID = f.read('ID')
R = f.read('R')
Z = f.read('Z')
THETA = f.read('Theta')
PSI = f.read('Psi')

fOut = open(outName, 'w')
fOut.write('ID, R, Z, THETA, PSI\n')
fOut = open(outName, 'a') ##open in append mode now.

n = len(ID)
for i in range(0,n) :
    if i % 500 == 0 : print('Processing field line ', i)

    id = ID[i]
    validIdx = np.where(id >= 0)
    validIds = id[validIdx]
    r = R[i][validIdx]
    z = Z[i][validIdx]

    theta = THETA[i][validIdx]
    psi = PSI[i][validIdx]
    print(len(id), len(psi), len(theta))


    data = list(zip(validIds, r, z, theta, psi))
    np.savetxt(fOut, data, delimiter=',', fmt='%d, %12.10lf, %12.10lf, %12.10lf, %12.10lf')
