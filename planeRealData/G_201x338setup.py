# 
#.. create <project>.geo file form this file
#
#     mkGeoWithData2D.py test
#
#.. create 3D mesh <meshfile>:
#
#  gmsh -3 -o testmesh.msh test.geo
#
project="G_201x338test"
#
# It is assumed that the XY-data arrays are flat and parallel to the surface at a given height and
# corresponding data do not change vertically across a thin layer. 
#
# .... these are the horizontal coordinates of the lower-left (south-west) end of the data array in the mesh [m]:
DataRefX=0.0
DataRefY=0.0

# ... this is the height of the grav data above ground [m] (can be zero)
DataHeightAboveGround=0

# ... number of data points in east-west (X) and north-south (Y) direction:
DataNumX=201
DataNumY=338

# .... this spacing of the data array [m]:
LDataX=800*(DataNumX-1)
LDataY=800*(DataNumY-1)

# Note: the resolution specified here should roughly match the resolution of the actual data as input data are interpolated to the resolution in the mesh

# ... this is the "thickness" of the data array = the thickness of the vertical layer. 
DataMeshSizeVertical=800

# ... this is the thickness of region below the data area. In essence it defines the depth of the inversion
CoreThickness=100000

# ... there is also an air layer and this is its thickness [m] (no updates for density and magnetization here)
AirLayerThickness=100000

# ... there is padding around the core and air layer. For the subsurface there will be updates in padding region but not for the air layer  
PaddingAir=100000.0
PaddingX=80000.0
PaddingY=135000.0
PaddingZ=100000.0

# ... these are factors by which the DataMeshSizeVertical is raised in the air layer and in the core. 

MeshSizeAirFactor=10
MeshSizeCoreFactor=5

#  ... these are factors by which the core and air layer mesh size are raised for the padding zone. 
MeshSizePaddingFactor=5


# name of the mesh file (gmsh 3 file format)
meshfile='G_201x338test.msh'


