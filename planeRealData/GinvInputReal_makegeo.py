#
# gambit
#
# Configuration file for gravity inversion for use by planeGravInv.py 
# mesh has been made with mkGeoWithData2D.py
#
# Inversion constants:
#
# scale between misfit and regularization
mu       = 1.e-14
#
# used to scale computed density.  kg/m^3    
rho_0    = 1.
#
# IPCG tolerance *|r| <= atol+rtol*|r0|*  (energy norm)
# absolute tolerance for IPCG interations
atol     = 0.        
#
# relative tolerance for IPCG iterations   
rtol     = 1.e-3
#
# tolerance for solving PDEs 
# make sure this is not more than the square of rtol     
pdetol   = 1.e-10    
#
# maximum number of IPCG iterations
iter_max = 500 
#
# data scale.  Program assumes m/s^2, 
# converts micrometres/s^2 to m/s^2
data_scale = 1.e-6
#
#
# File names
# mesh file name.  This needs to be in msh or fly format
mesh_name = "G_201x338test.fly"
#
# data file name in netcdf format.  See readme.md for more details.
data_file = "Gravity_201x338.nc"
#
# output file name for .csv output and silo output
output_name = "G_test_201x338_rho0_{0:1.3e}_mu_{1:1.3e}".format(rho_0,mu)
#
#
# Level for the verbosity of the output, "low", "medium" or "high".
# low: 
#   screen outputs:
#      data range, 
#      summaries of gravity data and final gravity
#      initial, final and difference misfits
#   file output:
#      silo of final solution
# medium: low outputs + 
#   screen outputs:
#      residual norm from the IPCG iterations
# high: medium outputs + 
#   screen outputs:
#      misfit and smoothing value at each iteration step
#   file outputs:
#      csv files for misfit and smoothing at each IPCG iteration
#      silos at misfit values of 0.05, 0.01, 0.008 and 0.005. (Initial misfit is 0.5.)
VerboseLevel = "low"
#VerboseLevel = "medium"
#VerboseLevel = "high"
