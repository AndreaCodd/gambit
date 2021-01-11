#!/usr/bin/python3
import esys.escript.unitsSI as U
import importlib, sys, os
sys.path.insert(0, os.getcwd())
import argparse
parser = argparse.ArgumentParser(description='Inputs needed to run gravity inversion for plane data.',epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')

mu       = 1.e-6
k_0      = 1.*U.kg/U.m**3
rho_r    = 2650.*U.kg/U.m**3
atol     = 0.        # absolute tolerance  
rtol     = 1.e-3     # relative tolerance   *|r| <= atol+rtol*|r0|*  (energy norm)
pdetol   = 1.e-10    # make sure this is not more than the square of rtol
iter_max = 500 
data_scale = 1.

Magmag =  58014.0333    #nT
Incl   = -64.68784390   # positive is “down”
Decl   =  0.86171459    # clockwise from true north

mesh_name = "Gravity_201x338.fly"
data_file = "Magnetic_201x338.nc"
output_name = "M_201x338_k0_{0:1.3e}_mu_{1:1.3e}".format(k_0,mu)


VerboseLevel = "low"
#VerboseLevel = "medium"
#VerboseLevel = "high"
