# gambit - gravity inversion for point and plane data and magnetic inversion for plane data
A package for geophysical inversion of point or plane gravity data and plane magnetic data.  It is based on [esys-escript](https://github.com/esys-escript/esys-escript.github.io) with python3 and supports parallel execution with domain decomposition for threading as well as MPI. The inversion code is for plane data in netcdf format and point data in csv format.  Code is provided to assist in constructing a mesh using [GMSH](https://gitlab.onelab.info/gmsh/gmsh).

# Examples
[planeRealData](https://github.com/AndreaCodd/gambit/blob/main/planeRealData/readme.md)

[planeSyntheticData](https://github.com/AndreaCodd/gambit/blob/main/planeSyntheticData/readme.md)

[pointData](https://github.com/AndreaCodd/gambit/blob/main/pointData/readme.md)
