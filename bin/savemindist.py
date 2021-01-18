__copyright__ = "Copyright (c) 2021 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Andrea Codd"

import importlib, os, sys
sys.path.append(os.getcwd())
import argparse
import numpy as np
from esys.escript import length

parser = argparse.ArgumentParser(description='This computes nearest neighbour distance for observation points and saves it as a csv file.  It is used to make the variable ground mesh.', epilog="version 01/2021 by a.codd@uq.edu.au")
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
