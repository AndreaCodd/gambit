__copyright__ = "Copyright (c) 2021 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Andrea Codd"

import importlib, sys, os
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='this program makes a subset of the data points', epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')

args = parser.parse_args()
config = importlib.import_module(args.config)

# import measuring locations
pts = np.loadtxt(config.big_obsPts_file, delimiter=',')
gravs = np.loadtxt(config.big_gravity_data_file, delimiter=',')
accs = np.loadtxt(config.big_acc_data_file, delimiter=',')
pick = config.pick

print('The orginal data sets were for ',len(pts),' observation points.')

smallpts = np.array(pts[::pick])
smallgravs = np.array(gravs[::pick])
smallaccs = np.array(accs[::pick])

np.savetxt("Grav_small_MeasEs.csv",smallpts, delimiter = ',')
np.savetxt("Grav_small_gz.csv",smallgravs, delimiter = ',')
np.savetxt("Grav_small_acc.csv",smallaccs, delimiter = ',')

print('The smaller data sets have ',len(smallpts), ' observation points.')
