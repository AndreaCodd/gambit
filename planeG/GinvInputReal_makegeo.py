#!/usr/bin/python3
import esys.escript.unitsSI as U
import importlib, sys, os
sys.path.insert(0, os.getcwd())
import argparse
parser = argparse.ArgumentParser(description='Inputs needed to run gravity inversion for plane data.',epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')

mu       = 1.e-14
rho_0    = 1.*U.kg/U.m**3
rho_r    = 2650.*U.kg/U.m**3
atol     = 0.        # absolute tolerance  
rtol     = 1.e-3     # relative tolerance   *|r| <= atol+rtol*|r0|*  (energy norm)
pdetol   = 1.e-10    # make sure this is not more than the square of rtol
iter_max = 500 
data_scale = 1.e-6

mesh_name = "G_201x338test.fly"
data_file = "Gravity_201x338.nc"
output_name = "G_test_201x338_rho0_{0:1.3e}_mu_{1:1.3e}".format(rho_0,mu)


#VerboseLevel = "low"
#VerboseLevel = "medium"
VerboseLevel = "high"
