#!/usr/bin/python3
import esys.escript.unitsSI as U
import importlib, sys, os
sys.path.insert(0, os.getcwd())
import argparse
parser = argparse.ArgumentParser(description='Inputs needed to run gravity inversion for plane data.',epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')

mu       = 1.e-8
k_0    = 1.*U.kg/U.m**3
rho_r    = 2650.*U.kg/U.m**3
atol     = 0.        # absolute tolerance  
rtol     = 1.e-4     # relative tolerance   *|r| <= atol+rtol*|r0|*  (energy norm)
pdetol   = 1.e-10    # make sure this is not more than the square of rtol
iter_max = 500 
data_scale = 1.

Magmag = 45000.0
Incl = 0.0
Decl = -90.

mesh_name = "syth_testmesh.fly"
data_file = "test_mag.nc"
output_name = "M_synthetic_test_k0_{0:1.3e}_mu_{1:1.3e}".format(k_0,mu)


#VerboseLevel = "low"
#VerboseLevel = "medium"
VerboseLevel = "high"
