

input_file = 

mu       = 1.e-14
rho_0    = 1.*U.kg/U.m**3
rho_r    = 2650.*U.kg/U.m**3
atol     = 0.        # absolute tolerance  
rtol     = 1.e-2     # relative tolerance   *|r| <= atol+rtol*|r0|*  (energy norm)
pdetol   = 1.e-10    # make sure this is not more than the square of rtol
iter_max = 500 
data_scale = 1.e-6

mesh_name = "Gravity_201x338.fly"
data_file = "Gravity_201x338.nc"
output_name = "G_201x338_rho0_{0:1.3e}_mu_{1:1.3e}".format(rho_0,mu)


VerboseLevel = "low"
#VerboseLevel = "medium"
#VerboseLevel = "high"
