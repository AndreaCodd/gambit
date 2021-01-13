import importlib, os, sys
sys.path.append(os.getcwd())
import argparse
import numpy as np
from esys.escript import length

parser = argparse.ArgumentParser(description='this computes nearest neighbour distance for observation points', epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')

args = parser.parse_args()
config = importlib.import_module(args.config)

print('computing closest neighbour distance for observation points from file ',config.obsPts_file)


# import measuring locations
pts=np.loadtxt(config.obsPts_file, delimiter=',')
mdist=[]
for ind1 in range(len(pts)):
    dist = config.maxDist
    for ind2 in range(len(pts)):
        if ind1 != ind2:
            dist2= length(pts[ind1]-pts[ind2])
            if dist2 < dist :
                dist = dist2
    mdist.append(dist)
mdist=np.array(mdist)
np.savetxt(config.minDist_file, mdist, delimiter=",")
print('nearest neighbour distance saved in file ',config.minDist_file)
print('minimum distance ', min(mdist), ' maximum distance', max(mdist))
