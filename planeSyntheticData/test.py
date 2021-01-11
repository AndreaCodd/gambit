# 
#.. create syth_test.geo file form this file
#
#     python3  mkSyntheticGeoData2D.py test
#
#.. create 3D mesh syth_testmesh.msh
#
#  gmsh -3 -format msh2 -o syth_testmesh.msh syth_test.geo
#
#.. run for synthetic data 
#
#    python3 mkSyntheticData2D.py -s synth -g -m test
#
#  creating test_grav.nc and test_mag.nc
#
#  next step is start the inversion process:
#
#    python3 mkGeoFromNc.py -c 60 -p 60 -h 1. -P benchmark1 -g test_grav.nc -m test_mag.nc
#
project="test"
km=1000.

gravfile=project+"_grav.nc"
magfile=project+"_mag.nc"

#
# It is assumed that the XY-data arrays are flat and parallel to the surface at a given height and
# corresponding data are not changing vertically across a thin layer. 
#
# .... these are the horizontal coordinates of the lower-left (south-west) end of the data array in the mesh [m]:
DataRefX=0.0
DataRefY=0.0

# ... this is the height of the grav and magnetic data above ground [m] (can be zero)
DataHeightAboveGround=0.#1*km

# .... this total extent of the data array [m]:(total length of the array)
#DataSpacingX=1*km
#DataSpacingY=1*km
LDataY=40*km
LDataX=70*km

# ... number of data points in east-west (X) and north-south (Y) direction:
DataNumX=71
DataNumY=41

# Note: the resolution specified here should roughly match the resolution of the actual data as input data are interpolated to the resolution in the mesh

# ... this is the "thickness" of the data array = the thickness of the vertical layer. 
DataMeshSizeVertical=1*km

# ... this is the thickness of region below the data area. In essence it defines the depth of the inversion
CoreThickness=60*km

# ... there is also an air layer and this is its thickness [m] (no updates for density and magnetization here)
AirLayerThickness=30*km

# ... there is padding around the core and air layer. For the subsurface there will be updates in padding region but not for the air layer  
PaddingAir=60000.0
PaddingX=60000.0
PaddingY=60000.0
PaddingZ=60000.0

# ... these are factors by which the DataMeshSizeVertical is raised in the air layer and in the core. 

MeshSizeAirFactor=10
MeshSizeCoreFactor=5

#  ... these are factors by which the core and air layer mesh size are raised for the padding zone. 
MeshSizePaddingFactor=5


# name of the mesh file (gmsh 3 file format)
meshfile=project+'mesh.msh'


B_hx = 45000.0 # background magnetic field in nT x direction
B_hz = 0.0     # background magnetic field in nT y direction
B_hy = 0.0     # background magnetic field in nT z direction
#
# this defines the assumed true density and magnetization:
# 
s1={ 'xc' : LDataX/3, 'yc' : LDataY/3, 'zc' : -CoreThickness*0.2, 'r' : 10*km }
s2={ 'xc' : 2*LDataX/3, 'yc' : 2*LDataY/3, 'zc' : -CoreThickness*0.1, 'r' : 8*km }

# ... 500 kg/m^3 over the union of sphere 1 and 2
true_density = [(-320, [s1]),(500, [s2]) ] 

# ... 0.1 on sphere 1 and 0.03 on sphere 2:
true_magnetization= [ ( 0.1, [s1]), (0.03, [s2])]
