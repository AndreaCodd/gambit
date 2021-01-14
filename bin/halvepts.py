import importlib, sys, os
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='this generates gmsh geo file for synthetic application', epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')

args = parser.parse_args()
config = importlib.import_module(args.config)

# import measuring locations
pts = np.loadtxt(config.obsPts_file, delimiter=',')
gravs = np.loadtxt(config.gravity_data_file, delimiter=',')
accs = np.loadtxt(config.accuracy_data_file, delimiter=',')



smallpts = np.array(pts[::4])
smallgravs = np.array(gravs[::4])
smallaccs = np.array(accs[::4])

np.savetxt("Grav_small_MeasEs.csv",smallpts, delimiter = ',')
np.savetxt("Grav_small_gz.csv",smallgravs, delimiter = ',')
np.savetxt("Grav_small_acc.csv",smallaccs, delimiter = ',')
