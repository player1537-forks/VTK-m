import adios2
import sys
import matplotlib.pyplot as plt
import numpy as np

print(sys.argv)

if len(sys.argv) != 3 :
    print('usage: ', sys.argv[0], ' infile skip')
    sys.exit()

inFile = sys.argv[1]
skip = int(sys.argv[2])

f=adios2.open(inFile,'r')
r=f.read('r')
z=f.read('z')
psi=f.read('psi')
theta=f.read('theta')
f.close()

#f = adios2.open('xgc.particle.0000000.init.bp', 'r')
f = adios2.open('xgc.particle.init.bp', 'r')
seedIds = f.read('egid')
ephase = f.read('ephase')

def dumpPoints(ID, R,Z) :
    if ID == 0 :
        f = open('XGC.RZ.txt', 'w')
        f.write('ID, R, Z, PUNC\n')
    else:
        f = open('XGC.RZ.txt', 'a')

    n = len(R)
    for i in range(0,n) :
        f.write('%d, %15.14f, %15.14f, %d\n' % (ID, R[i], Z[i], i))
    f.close()


ns = np.shape(r)[0]



param = {'c':'0.2', 'alpha':0.5}

plt.figure(figsize=[10,10])
print('ns= ', ns)
i0 = 7300
i1 = 7301

print('i0/i1= ', i0, i1, skip)
R = range(0, ns, skip)
#R = [100, 500, 4000, 6000, 7300]
#R = range(0, ns, skip)
R = [7300]

for i in R :
    tmp=r[i,:]
    msk=np.nonzero(tmp)
    # print (i, len(msk[0]))
    rt=r[i,msk]
    zt=z[i,msk]
    p = np.where(seedIds == i)[0][0]
    print('p= ', p)
    seed = (ephase[p][0], ephase[p][1], ephase[p][2] )

    #print(i, j, ': R,Z=', rt[0][0], zt[0][0])
    #print('%16.14f, %16.14f, 0' %(rt[0][0], zt[0][0]))
    print('Seed= %16.14f, %16.14f, %16.14f' %(seed[0], seed[1], seed[2]))
    dumpPoints(i, rt[0], zt[0])

    plt.scatter(rt, zt, s=3, edgecolor='none', **param)
    plt.axis('scaled');

plt.show()
