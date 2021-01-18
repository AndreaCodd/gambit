#
# gambit
#
gravity_data_file = "Grav_gz.csv"
acc_data_file = "Grav_acc.csv"
obsPts_file = "Grav_MeasEs.csv"
mesh_name = "smallPointMesh.fly"

maxDist = 5000

rho_0  = 1.
atol   = 0.
rtol   = 1.e-2
pdetol = 1.e-8
iter_max = 500
s = 900.
a = 72000.
b = 8.    

data_scale = 1.e-6

output_name = 'pointG_big_s_'+str(int(s))+'_a_'+str(int(a))+'_b_'+str(int(b))+'_'

###  dataweighting # equal, relative, accuracy
dataWt = "accuracy" 

# assumes ground level is horizontal at z=0 
###  depthWeight    # noWt, coreWt, baseWt, updown
depthWeight = "noWt" 
#coreD = -20000.

VerboseLevel = "low"
#VerboseLevel = "medium"
#VerboseLevel = "high"






