### Real point gravity data example
Files:
  * Grav\_MeasEs.csv
  * Grav\_gz.csv
  * Gravacc.csv
  * pointG\_Meshconfig.py
  * makeVariableGroundMesh.py
  * GPointinvInputReal.py
  * mkfly.py
  * pointGravityInversion.py
  
Steps:
1. compute minimum distance between observation points
`python3 ~/gambit/bin/savemindist.py pointG\_Meshconfig`
2. generate geo file
`python3 ~/gambit/bin/makeVariableGroundMesh.py pointG\_Meshconfig`
3. generate mesh from geo file
`gmsh -3 General.ExpertMode = 0 -o PointMesh.msh PointMesh.geo`
`gmsh -3 -format msh2 -o PointMesh.msh PointMesh.geo`
4. convert msh file to fly file (optional)
`run-escript ~/gambit/bin/mkfly.py PointMesh`
5. run inversion
`run-escript ~/gambit/bin/planeGravInv.py GinvInputReal_withgeo`
  
### Synthetic point gravity data example using observation points and mesh from real example
Files:
  * Grav\_MeasEs.csv
  * Grav\_gz.csv
  * Gravacc.csv
  * pointG\_Meshconfig.py
  * makeVariableGroundMesh.py
  * GPointinvInputReal.py
  * mkfly.py
  * pointGravityInversion.py
  
Steps:
1. compute minimum distance between observation points
`python3 ~/gambit/bin/savemindist.py pointG\_Meshconfig`
2. generate geo file
`python3 ~/gambit/bin/makeVariableGroundMesh.py pointG\_Meshconfig`
3. generate mesh from geo file
`gmsh -3 General.ExpertMode = 0 -o PointMesh.msh PointMesh.geo`
`gmsh -3 -format msh2 -o PointMesh.msh PointMesh.geo`
4. convert msh file to fly file (optional)
`run-escript ~/gambit/bin/mkfly.py PointMesh`
5. run inversion
`run-escript ~/gambit/bin/planeGravInv.py GinvInputReal_withgeo`

## Comments
1. An older output format for gmsh needs to be used so that tagging works.
2. It is not necessary to convert msh files to fly files but in some cases it is faster to load fly files in the inversion code.  The inversion code is generally run multiple times with different scaling factors.
3. Three levels of output are possible:
    + low: outputs data range, summaries of gravity data and final gravity, and initial, final and difference misfits as well as a silo of the final solution.
    + medium: low outputs + outputs residual norm from the PCG iterations
    + high: medium outputs + outputs misfit and smoothing value at each iteration step.  Also saves silos at misfit values of 0.05, 0.01, 0.008 and 0.005.  (Initial misfit is 0.5.)
4. run-escript can be run with threads.


