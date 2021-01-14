#
# gambit
#

gravity_data_file = "Grav_small_gz.csv"
accuracy_data_file = "Grav_small_acc.csv"
obsPts_file = "Grav_small_MeasEs.csv"
mesh_name = "smallPointMesh.fly"


rho_0  = 1.
atol   = 0.
rtol   = 1.e-3
pdetol = 1.e-8
iter_max = 500
s = 900.
a = 72000.
b = 8.    

data_scale = 1.e-6

save_name = 'pointG_s_'+str(int(s))+'_a_'+str(int(a))+'_b_'+str(int(b))+'_'

###  dataweighting # equal, relative, accuracy
dataWt = "accuracy" 

# assumes ground level is horizontal at z=0 
###  depthWeight    # noWt, coreWt, baseWt, updown
depthWeight = "noWt" 
#coreD = -20000.




