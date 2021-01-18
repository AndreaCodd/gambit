#
# gambit
#
# This file contains the information to make the mesh for the original point data set.
#
#
# The original data sets are in csv format and are for 6519 observation points.
# 1. `Grav_MeasEs.csv` contains the cartesian coordinates of the observation points (m).
# 2. `Grav_gz.csv` contains the Bouger corrected gravity measurements ( micro m/s^2).
# 3. `Grav_acc.csv` cotains the measurement accuracy (micro m/s^2).
#
gravity_data_file = "Grav_gz.csv"
accuracy_data_file = "Grav_acc.csv"
obsPts_file = "Grav_MeasEs.csv"
#
# The inversion code has a measurement scale factor so it is not absolutely necessary
# to have the same units for the gravity measurments.  The scale factor conversts the 
# measurements to m/s^2.  It is assumed that the spatial measurements are in m.
#
# Element size at observation points depends on the distance to its nearest neighbour. 
# To make the ground variable mesh, the minimum nearest neighbour distance is computed for each observation point
# and recorded in the file 
minDist_file = "Grav_minDist.csv"
#
# GMSH mesh making parameters
#
# A bound for the maximum nearest neighbour distance.  This value will be used to compute element size for any point
# that has no nearest neighbour closer.  It does not need to be the maximum distance.
maxDist = 5000
#
#  For each observatio point, element length is, in general
#      nearest neighbour distance / spacing0 
# 
#  For nearest neighbour distance < mindist1, element length is, 
#      nearest neighbour distance / spacing1 
#
#  For nearest neighbour distance < mindist1, element length is, 
#      nearest neighbour distance / spacing2 
#
mindist1 = 100
mindist2 = 20
spacing0 = 6
spacing1 = 4
spacing2 = 2
#
# The core area is 1.01*span of observation points in the x and y directions.
# Ground level is at 0m and is horizontal.  There is a core ground region and core air region.
groundLevel = 0
coreDepth = - 20000
coreAir = 10000
#
# mesh sizes vary.
# Mesh size at the top of the core air region
MCtop = 3000
# Mesh size at the bottom of the core ground region
MCbase = 2000
# mesh size at the core ground region.  This is probably over ruled by the mesh size of the observation points.
MCground = 2000
#
# The buffer region has different element lengths.
# Mesh size at the tob of the buffer in the air
MBtop = 20000
# Mesh size at the edge of the buffer region at ground level
MBbase = 20000
# Mesh size at the base of the buffer region.
MBground = 10000
# Name for the geo file.
geo_name = "PointMesh.geo"









