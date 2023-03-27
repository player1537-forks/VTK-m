import sys
import numpy as np
import adios2 as ad
import random

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import mpl_scatter_density
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

norm = ImageNormalize(vmin=0.0, vmax=200, stretch=LogStretch())

from matplotlib.colors import ListedColormap, LinearSegmentedColormap


inFiles = []
for f in sys.argv[1:-1] :
    if '.bp' not in f :
        print('Error: BP file not specified.')
    inFiles.append(f)

if '.png' not in sys.argv[-1]  :
    print('Error! Last argument must be an image file')
    sys.exit()


def readFile(inFile, rIn=[], zIn=[]) :
    print('reading file:', inFile)
    f = ad.open(inFile, "r")
    r = f.read("R")
    z = f.read("Z")
    # psi = f.read("psi")
    # theta = f.read("theta")
    gid = f.read("ID")
    time = f.read("TimeStep")
    f.close()
    print('Read is done!')

    gid += 1
    validMask = np.nonzero(gid)
    R = r[validMask]
    Z = z[validMask]
    ID = gid[validMask]

    N = len(R)
    print('N=', N)
    idx = list(range(N))
    random.seed(43)
    random.shuffle(idx)

    R0 = R
    Z0 = Z
    for i in range(N) :
        R[i] = R0[idx[i]]
        Z[i] = Z0[idx[i]]
    if len(rIn) > 0 :
        R = np.concatenate((rIn, R), axis=None)
        Z = np.concatenate((zIn, Z), axis=None)
    return (R,Z)


def myplot(Rvals, Zvals, cmapName, outFile) :
    ns = np.shape(Rvals)[0]
    param = {"alpha": 0.5}

    cmap = plt.get_cmap(cmapName)
    print("colormap.N:", cmap.N, ns, cmapName)

    newcolors = cmap(np.linspace(0, 1, cmap.N))
    empty = np.array([0, 0, 0, 0])
    newcolors[:1, :] = empty
    cmap2 = ListedColormap(newcolors)

    fig = plt.figure(figsize=[12, 10])
    ax = fig.add_subplot(1, 1, 1, projection="scatter_density")
    density = ax.scatter_density(Rvals, Zvals, norm=norm, cmap=cmap, zorder=10)

    density = ax.scatter_density(Rvals, Zvals, norm=norm, cmap=cmap2, zorder=30)
    fig.colorbar(density, label="Number of points per pixel")
    plt.axis("scaled")
    plt.xlim([4.9, 5.4])
    plt.ylim([-3.6, -3.1])

    #plt.show()
    plt.savefig(outFile, format='png')


R = Z = []
for f in sys.argv[1:-1] :
    (R,Z) = readFile(f, R, Z)

print('len=', len(R), len(Z))
myplot(R, Z, 'hot', sys.argv[-1])
