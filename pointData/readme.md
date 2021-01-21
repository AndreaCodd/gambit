### Point gravity inversion - real data sets
The original data sets are in csv format and are for 6519 observation points.
1. `Grav_MeasEs.csv` contains the cartesian coordinates of the observation points.
2. `Grav_gz.csv` contains the Bouger corrected gravity measurements in micro m/s^2.
3. `Grav_acc.csv` cotains the measurement accuracy in micro m/s^2.

It is possible to run the inversion on a smaller data set.  The smaller data set has 1630 data points and can be created using halvepts.py.
 - `python3 ~/gambit/bin/halvepts.py pointG_small_Meshconfig`
    * `Grav_small_MeasEs.csv` contains the cartesian coordinates of the observation points.
    * `Grav_small_gz.csv` contains the Bouger corrected gravity measurements in micro m/s^2.
    * `Grav_small_acc.csv` cotains the measurement accuracy in micro m/s^2.

Two configuration files are needed, one for the mesh and one for the inversion
- pointG_small_Meshconfig.py or pointG_Meshconfig.py
- pointG_small_Realconfig.py or pointG_Realconfig.py

The inversion can either be run with the small data sets or the large data sets. 

Steps:
1. Compute minimum distance between observation points. Element length at observations is determined by the distance to its nearest neighbour.
 - `python3 ~/gambit/bin/savemindist.py pointG_small_Meshconfig`
      * Grav_small_minDist.py
 - `python3 ~/gambit/bin/savemindist.py pointG_Meshconfig`
      * Grav_minDist.py
2. Generate geo file.  This uses mesh config files pointG_small_Meshconfig.py or pointG_Meshconfig.py
 - `python3 ~/gambit/bin/makeVariableGroundMesh.py pointG_small_Meshconfig`
 - `python3 ~/gambit/bin/makeVariableGroundMesh.py pointG_Meshconfig`
3. Generate mesh from geo file.
 - `gmsh -3 -o smallPointMesh.msh smallPointMesh.geo`
 - `gmsh -3 -format msh2 -o PointMesh.msh PointMesh.geo`
4. Convert msh file to fly file (optional).
 - `run-escript ~/gambit/bin/mkfly.py smallPointMesh`
 - `run-escript ~/gambit/bin/mkfly.py PointMesh`
5. Run inversion.
 - `run-escript ~/gambit/bin/pointGravityInversion.py pointG_small_Realconfig`
 - `run-escript ~/gambit/bin/pointGravityInversion.py pointG_Realconfig`
  
### Synthetic point gravity data example using observation points, accuracy measurements and mesh from real example.
A synthetic example for the observation locations and accuracies from the real example but with 4 regions with non zero density, 1 sphere, 2 hemispheres and one large rectangular prism.  The mean density is 0kg/m^3.  Make the mesh as for the real example
  
Steps:
1. compute gravity measurements
`python3 ~/gambit/bin/savemindist.py pointG_Meshconfig`
2. run inversion
`run-escript ~/gambit/bin/planeGravInv.py GinvInputReal_withgeo`

## Comments
1. For the element lengths defined in the given geo codes, the small mesh has 205 thousand vertices and 1.4 million elements and the large mesh has 2650 thousand vertices and 18 million elements, a little over 10 times bigger.
2. It is not necessary to convert msh files to fly files but in some cases it is faster to load fly files in the inversion code.  The inversion code is generally run multiple times with different scaling factors.
3. Three levels of output are possible:
    + low: outputs data range, summaries of gravity data and final gravity, and initial, final and difference misfits as well as a silo of the final solution.
    + medium: low outputs + outputs residual norm from the PCG iterations
    + high: medium outputs + outputs misfit and smoothing value at each iteration step.  Also saves silos at misfit values of 0.05, 0.01, 0.008 and 0.005.  (Initial misfit is 0.5.)
4. run-escript can be run with threads.   `run-escript -t10 ~/gambit/bin/pointGravityInversion.py pointG_Realconfig`


