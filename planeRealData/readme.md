### Real example with existing geo file
Files:
  * Gravity\_201x338.nc
  * Gravity\_201x338.geo
  * GinvInputReal\_withgeo.py
  * mkfly.py
  * planeGravInv.py
Steps:
1. generate mesh from geo file  
`gmsh -3 -format msh2 -o Gravity_201x338.msh Gravity_201x338.geo`
2. convert msh file to fly file (optional)
`run-escript ~/gambit/bin/mkfly.py Gravity_201x338`
3. run inversion
`run-escript ~/gambit/bin/planeGravInv.py GinvInputReal_withgeo`

### Real example make geo file
Files:
  * Gravity_201x338.nc
  * G_201x338setup.py
  * OnGroundTemplate.geo
  * mkGeoWithData2D.py
  * GinvInputReal_makegeo.py
  * mkfly.py
  * planeGravInv.py
Steps:
1. generate geo
`python3 mkGeoWithData2D.py G_201x338setup`
2. generate mesh
`gmsh -3 -format msh2 -o G_201x338test.msh G_201x338test.geo`
3. convert .msh format to .fly (optional)
`run-escript mkfly.py G_201x338test`
4. run inversion
`run-escript planeGravInv.py GinvInputReal_makegeo`
  
### Synthetic example
Files:
  * test.py 
  * OnGroundTemplate.geo 
  * AboveGroundTemplate.geo
  * mkSyntheticGeoData2D.py
  * mkSyntheticData2D.py
  * GinvInput\_synth.py
  * mkfly.py
  * planeGravInv.py
Steps:
1. generate geo
`python3  mkSyntheticGeoData2D.py test`
2. generate mesh
`gmsh -3 -format msh2 -o syth_testmesh.msh syth_test.geo`
3. convert msh format to fly (optional)
`run-escript mkfly.py syth_testmesh`
4. make synthetic data
`python3 mkSyntheticData2D.py -s synth -g -m test`
5. run inversion
`run-escript planeGravInv.py GinvInput_synth.py`

## Comments
1. An older output format for gmsh needs to be used so that tagging works.
2. It is not necessary to convert msh files to fly files but in some cases it is faster to load fly files in the inversion code.  The inversion code is generally run multiple times with different scaling factors.
3. Three levels of output are possible:
    + low: outputs data range, summaries of gravity data and final gravity, and initial, final and difference misfits as well as a silo of the final solution.
    + medium: low outputs + outputs residual norm from the PCG iterations
    + high: medium outputs + outputs misfit and smoothing value at each iteration step.  Also saves silos at misfit values of 0.05, 0.01, 0.008 and 0.005.  (Initial misfit is 0.5.)
4. run-escript can be run with threads.
