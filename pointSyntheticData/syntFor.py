
from esys.escript import *
import esys.escript as escript
import esys.finley as finley
import esys.escript.unitsSI as U
import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDE,  SolverOptions
from esys.escript.pdetools import PCG
from esys.downunder import *
from esys.weipa import *
from esys.escript import integrate, wherePositive, length
#from esys.escript import saveDataCSV, loadDataCSV
#from esys.escript.pdetools import MaskFromTag
import sys
#from scipy.io import netcdf_file
from time import time
import pickle as pickle



mesh_file = 'variable_groundfiner.fly' #'vg_fine9.fly'
#mesh_file='MtIsaCloncurryMesh3.fly'

MeasEs = np.loadtxt('MtIsaCloncurry_MeasEs.csv',delimiter=',')
recorders=[]
for bob in MeasEs:
    recorders.append((bob[0],bob[1],bob[2]))

# build domain
print("Loading the mesh..."+mesh_file)
dom = finley.ReadMesh(mesh_file, numDim=3,
                         diracPoints = [sp for sp in recorders],
                         diracTags = [st for st in range(len(recorders))])

coord=dom.getX()
depthWt = whereNegative(coord[2])
minZ=inf(coord[2])

d1 = -40
d2 = -145.239923051
d3 = 130
d4 = 50

R0 = 10000
Cx0 = -60000
Cy0 = 15000
Cz0 = -1000

R1 = 10000
Cx1 = 50000
Cy1 = 20000
Cz1 = -10000

R2 = 15000
Cx2 = 60000
Cy2 = -15000
Cz2 = 0

dens1 = d1*(whereNegative(coord[0]-20000)*wherePositive(coord[0]+20000)*
               whereNegative(coord[1]-20000)*wherePositive(coord[1]+20000)*
               wherePositive(coord[2]+10000)*whereNegative(coord[2]+1000))
dens2 = d2*(whereNegative((coord[0]-Cx0)**2 + (coord[1]-Cy0)**2 + (coord[2]-Cz0)**2- R0**2)*
               whereNegative(coord[2]-Cz0))
dens3 = d3*(whereNegative((coord[0]-Cx1)**2 + (coord[1]-Cy1)**2 + (coord[2]-Cz1)**2- R1**2)*
                whereNegative(coord[2]+2000))
dens4 = d4*whereNegative((coord[0]-Cx2)**2 + (coord[1]-Cy2)**2 + (coord[2]-Cz2)**2- R2**2)*depthWt

density=dens1+dens2+dens3+dens4
totdens=integrate(density)
print(d2,totdens)

#saveSilo("syntheticDens", density=density)

FOSLSpde = LinearPDE(dom,numEquations=3,numSolutions=3)
FOSLSpde.setSymmetryOn()
qtop = whereZero(coord[2]-sup(coord[2]))
qbottom = whereZero(coord[2]-inf(coord[2]))
qleft = whereZero(coord[0]-inf(coord[0]))
qright = whereZero(coord[0]-sup(coord[0]))
qfront = whereZero(coord[1]-inf(coord[1]))
qback = whereZero(coord[1]-sup(coord[1]))

q=Data(0, (3,), Solution(dom))
q[0]=qleft+qright+qtop
q[1]=qfront+qback+qtop
q[2]=qbottom
FOSLSpde.setValue(q=q)
A=Data(0,(3,3,3,3),Function(dom))
X = Data(0,(3,3),Function(dom))
for jj in range(3):
    X[jj,jj] = 4.0*np.pi*U.Gravitational_Constant*density
    for kk in range(3):
       A[jj,jj,kk,kk] = Scalar(1.,Function(dom))
       if kk < jj:
          A[kk,jj,kk,jj] = Scalar(1.,Function(dom))
          A[jj,kk,jj,kk] = Scalar(1.,Function(dom))
          A[kk,jj,jj,kk] = -Scalar(1.,Function(dom))
          A[jj,kk,kk,jj] = -Scalar(1.,Function(dom)) 


FOSLSpde.setValue(X=X)
FOSLSpde.setValue(A=A)
Foptions=FOSLSpde.getSolverOptions()
Foptions.setPackage(SolverOptions.TRILINOS)
Foptions.setSolverMethod(SolverOptions.PCG)
Foptions.setPreconditioner(SolverOptions.AMG)       
Foptions.setTolerance(1e-8)
Foptions.setTrilinosParameter("number of equations",3)

U = FOSLSpde.getSolution()
Uz = -np.array(Locator(ContinuousFunction(dom),recorders)(U[2]))*1.e6   # add -ve for direction of measurement
np.savetxt("synthgz.csv", Uz, delimiter=",")
saveSilo("syntheticDens", density=density, gravity=-U[2])


gd = 0.0
xmin = inf(coord[0])
xmax = sup(coord[0])
ymin = inf(coord[1])
ymax = sup(coord[1])
gVol = (xmax-xmin)*(ymax-ymin)*(gd-minZ)


rho = density
rhobar = (1./gVol)*integrate(rho)
sigma = (1./gVol)*integrate((rho-rhobar)**2)
stddev=np.sqrt(sigma)
print('Min density  ', inf(rho))
print('Max density  ', sup(rho))
print('mean density ', rhobar)
print('variance     ', sigma)
print('stddev       ', stddev)  


