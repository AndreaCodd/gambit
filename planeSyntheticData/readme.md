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
`python3 ~/gambit/bin/mkSyntheticGeoData2D.py test`
2. generate mesh
`gmsh -3 -format msh2 -o syth_testmesh.msh syth_test.geo`
3. convert msh format to fly (optional)
`run-escript ~/gambit/bin/mkfly.py syth_testmesh`
4. make synthetic data
`python3 ~/gambit/bin/mkSyntheticData2D.py -s synth -g -m test`
5. run gravity inversion
`run-escript ~/gambit/bin/planeGravInv.py GinvInput_synth`

Output:
Three levels of output are possible:
+ low: outputs to the screen: data range, summaries of gravity data and final gravity, and initial, final and difference misfits.  Outputs final solution in a silo file that can be viewed in VisIt.   
+ medium: low outputs + outputs to the screen residual norm from the PCG iterations.
+ high: medium outputs + outputs to the screen misfit and smoothing values at each iteration step.  Also saves silos at misfit values of 0.05, 0.01, 0.008 and 0.005.  (Initial misfit is 0.5.)

## Comments
1. An older output format for gmsh needs to be used so that tagging works.
2. It is not necessary to convert msh files to fly files but in some cases it is faster to load fly files in the inversion code.  The inversion code is generally run multiple times with different scaling factors.
