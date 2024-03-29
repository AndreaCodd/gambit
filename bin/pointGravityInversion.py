#!/usr/bin/python3
__copyright__ = "Copyright (c) 2021 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Andrea Codd"

import importlib, sys, os
sys.path.insert(0, os.getcwd())
import argparse

from esys.escript import *
from esys.finley import ReadGmsh, ReadMesh
import esys.escript.unitsSI as U
import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDE,  SolverOptions
from esys.escript.pdetools import PCG
from esys.downunder import *
from esys.weipa import *


class FOSLSGravity(object):
    def __init__(self, domain, gz, recorders, rho_0, P0, wdsq, 
                       mu, a=0., b=1., atol=1.0, rtol=1.0, iter_max=100, 
                       pde_tol=1e-8, name='bob', verboseLevel="low"):
        self.domain = domain
        self.gz = np.array(gz)
        self.w = - kronecker(3)[2]
        self.locG = Locator(ContinuousFunction(self.domain),recorders)
        self.rho_0 = rho_0  
        self.P0 =  P0   
        self.wdsq = np.array(wdsq)
        self.mu = mu
        self.a = a    
        self.b = b
        self.atol = atol
        self.rtol = rtol
        self.iter_max = iter_max
        self.pdetol = pdetol
        self.name = name
        self.numes = len(self.gz)
        self.verboseLevel = verboseLevel
        self.beta = -4.0*np.pi*U.Gravitational_Constant

        #boundaries       
        coord=self.domain.getX()
        self.qtop = whereZero(coord[2]-sup(coord[2]))
        self.qbottom = whereZero(coord[2]-inf(coord[2]))
        self.qleft = whereZero(coord[0]-inf(coord[0]))
        self.qright = whereZero(coord[0]-sup(coord[0]))
        self.qfront = whereZero(coord[1]-inf(coord[1]))
        self.qback = whereZero(coord[1]-sup(coord[1]))

        # pdes
        self.dPpde = self.setupdPpde()
        self.dppde = self.setupdppde()
        self.FOSLSpde = self.setupFOSLSpde()

    def setupFOSLSpde(self):
        FOSLSpde = LinearPDE(self.domain,numEquations=3,numSolutions=3)
        FOSLSpde.setSymmetryOn()
        q=Data(0, (3,), Solution(self.domain))
        q[0]=self.qleft+self.qright+self.qtop
        q[1]=self.qfront+self.qback+self.qtop
        q[2]=self.qbottom
        FOSLSpde.setValue(q=q)
        A=Data(0,(3,3,3,3),Function(self.domain))
        for jj in range(3):
            for kk in range(3):
                A[jj,jj,kk,kk] = Scalar(1.,Function(self.domain))
                if kk < jj:
                    A[kk,jj,kk,jj] = Scalar(1.,Function(self.domain))
                    A[jj,kk,jj,kk] = Scalar(1.,Function(self.domain))
                    A[kk,jj,jj,kk] = -Scalar(1.,Function(self.domain))
                    A[jj,kk,kk,jj] = -Scalar(1.,Function(self.domain)) 
        FOSLSpde.setValue(A=A)
        Foptions=FOSLSpde.getSolverOptions()
        Foptions.setPackage(SolverOptions.TRILINOS)
        Foptions.setSolverMethod(SolverOptions.PCG)
        Foptions.setPreconditioner(SolverOptions.AMG)       
        Foptions.setTolerance(self.pdetol)
        Foptions.setTrilinosParameter("number of equations",3)
        Foptions.setTrilinosParameter("reuse: type","full")
        return FOSLSpde

    def setupdPpde(self):
        aa = self.a*self.a
        bb = self.b*self.b
        pde=LinearSinglePDE(self.domain, isComplex=False)
        pde.setValue(A = aa*kronecker(3))
        pde.setValue(D = Scalar(bb  , Function(self.domain)))
        pde.setSymmetryOn()
        q=self.qleft+self.qright+self.qfront+self.qback+self.qtop+self.qbottom
        pde.setValue(q=q)
        Foptions=pde.getSolverOptions()
        Foptions.setPackage(SolverOptions.TRILINOS)
        Foptions.setSolverMethod(SolverOptions.PCG)
        Foptions.setPreconditioner(SolverOptions.AMG)       
        Foptions.setTolerance(self.pdetol)
        Foptions.setTrilinosParameter("reuse: type","full")
        return pde


    def setupdppde(self):
        aabb = self.a*self.a*self.b*self.b
        pde=LinearSinglePDE(self.domain, isComplex=False)
        pde.setValue(A = aabb*kronecker(3))
        pde.setValue(D = Scalar(1., Function(self.domain)))
        q=self.qleft+self.qright+self.qfront+self.qback+self.qtop+self.qbottom
        pde.setValue(q=q)
        Foptions=pde.getSolverOptions()
        Foptions.setPackage(SolverOptions.TRILINOS)
        Foptions.setSolverMethod(SolverOptions.PCG)
        Foptions.setPreconditioner(SolverOptions.AMG)       
        Foptions.setTolerance(self.pdetol)
        Foptions.setTrilinosParameter("reuse: type","full")
        return pde

    def getRHSsolve(self, f):
        Y = Data(0,(3,),Function(self.domain))
        X = Data(0,(3,3),Function(self.domain))
        wtWd = Data(0,(3,),DiracDeltaFunctions(self.domain))
        for e in range(self.numes):
            bob = self.w*wdsq[e]*f[e] 
            wtWd.setTaggedValue(e, bob)
        self.FOSLSpde.setValue(X = X, Y = Y , y_dirac=wtWd )
        return self.FOSLSpde.getSolution()        
        

    def RHS(self): 
        U = self.getRHSsolve(self.gz) 
        NewY =Data(0.,(4,),Function(self.domain))
        NewY[3] =  (self.beta*self.rho_0/self.mu)*div(U)
        NewX = Data(0.,(4,3),Function(self.domain))
        return ArithmeticTuple (NewY, NewX)         

    def Aprod(self,P):
        # left hand side <AP,Q> = (SP, SQ)+(DP, Q) = (S1 P + D P, Q) + (S2 P, grad(Q)) 
        # returns tuple(S1 P + D P , S2 P)
        Udown = - self.getGravity(P)[2]
        U2pts = np.array(self.locG(Udown))
        cU = self.locG(U2pts)
        U =  self.getRHSsolve(cU) 

        a=self.a
        aa=a*a
        bb=self.b*self.b
        NewY =Data(0.,(4,),Function(self.domain))
        NewX =Data(0.,(4,3),Function(self.domain))

        gradp0 = grad(P[0])
        gradp1 = grad(P[1])
        gradp2 = grad(P[2])
        gradp3 = grad(P[3])
        curl2 = gradp1[0]-gradp0[1]
        curl1 = gradp0[2]-gradp2[0]
        curl0 = gradp2[1]-gradp1[2]
        y3 = -a*(gradp0[0] + gradp1[1] + gradp2[2]) + P[3]
        NewY[0] = bb*(P[0]-a*gradp3[0])
        NewY[1] = bb*(P[1]-a*gradp3[1])
        NewY[2] = bb*(P[2]-a*gradp3[2])
        NewY[3] = (self.beta*self.rho_0/self.mu)*div(U) + y3
        for kk in range(3):
            NewX[kk,kk] = -a*y3
            NewX[3,kk]= -a*NewY[kk]              # aa*gradp3[kk]-a*P[kk] 
        NewX[0,1] = -aa*curl2
        NewX[0,2] = -aa*curl1
        NewX[1,0] = aa*curl2            
        NewX[1,2] = -aa*curl0
        NewX[2,0] = aa*curl1
        NewX[2,1] = aa*curl0
        return ArithmeticTuple(NewY, NewX)     
        
    def Msolve(self,R):
        # solve for U, (S*S U,V) = (R[0],V) + (R[1],grad(V))
        # (U_i,V_i)+a^2(grad U_i,grad V_i) = (R[0]_i,V_i)+(R[1]_i,grad (V)_i)
        U =Data(0.,(4,),Solution(self.domain))
        for ind1 in range(3):
            Y=R[0][ind1]
            X=R[1][ind1]
            self.dPpde.setValue(Y=Y, X=X)
            U[ind1] = self.dPpde.getSolution()
        Y=R[0][3]
        X=R[1][3]
        self.dppde.setValue(Y=Y, X=X)
        U[3] = self.dppde.getSolution()
        return U       
  

    def bilinearform(self, P, R): 
        #print("bilinear form")  
        # R = (NewY , NewX)
        # returns (P , NewY ) + (grad(P) , NewX)   
        PR0= inner(P,R[0])
        PR1 = inner(grad(P),R[1])    # grad(a)[i,j,k] = partial a[i,j]/partial k  
        return integrate(PR0+PR1)
    
    def getGravity(self,P):
        Y = Data(0,(3,),Function(self.domain))
        X = Data(0,(3,3),Function(self.domain))
        wtWd = Data(0,(3,),DiracDeltaFunctions(self.domain))
        for jj in range(3):
            X[jj,jj] = self.beta*self.rho_0*P[3]
        self.FOSLSpde.setValue(Y = Y, X=X, y_dirac = wtWd)
        U = self.FOSLSpde.getSolution()
        return U  

    def getSPSP(self,P):
        a=self.a
        b=self.b
        gradp0 = grad(P[0])
        gradp1 = grad(P[1])
        gradp2 = grad(P[2])
        gradp3 = grad(P[3])
        SP0 = integrate((b*(P[0]-a*gradp3[0]))**2)
        SP1 = integrate((b*(P[1]-a*gradp3[1]))**2)
        SP2 = integrate((b*(P[2]-a*gradp3[2]))**2)
        SP3 = integrate((-a*(gradp0[0] + gradp1[1]+gradp2[2])+P[3])**2)
        SP4 = integrate((a*(-gradp1[2] + gradp2[1]))**2)
        SP5 = integrate((a*(gradp0[2] - gradp2[0]))**2)
        SP6 = integrate((a*(-gradp0[1] + gradp1[0]))**2)
        SPSP=SP0+SP1+SP2+SP3+SP4+SP5 +SP6
        return SPSP


    def myPCG(self, x,r,itermax,rtol):
        # x intial approximation P0 (4 elements)
        # r initial residual (4 elements but first 3 are zero)
        piter=0   # iteration count
        mfs = []
        smooths = []
        rzrzs = []
        rhat = self.Msolve(r)  
        d = rhat                                          
        rhat_dot_r = self.bilinearform(rhat, r)             
        if rhat_dot_r<0: print("negative norm.")
        rzrz0 = rhat_dot_r
        norm_r0=np.sqrt(rhat_dot_r)
        atol2=self.rtol*norm_r0
        if atol2<=0:
           print("Non-positive tolarance.")
        print(("PCG: initial residual norm = %e (absolute tolerance = %e)"%(norm_r0, atol2)))
        # this bit not actually needed, just initial output for csvs
        # will need to fix for varying down
        Udown = - self.getGravity(x)[2]
        U2pts = np.array(self.locG(Udown))
        diffG= U2pts - self.gz 
        mf=np.inner(diffG, diffG*self.wdsq )
        smooth=self.getSPSP(x)
        smooths.append(smooth)
        mfs.append(mf)
        rzrzs.append(1.0)
        print(piter,'mf',mf,'smooth',smooth, 'rzrz 1.0' )

        while not np.sqrt(rhat_dot_r) <= atol2:
           piter+=1
           if piter  >= iter_max: 
               print("maximum number of %s steps reached."%iter_max)
               break 
           q=self.Aprod(d)
           alpha = rhat_dot_r / self.bilinearform(d, q)
           x += alpha * d
           r += q * (-alpha)      
           rhat=self.Msolve(r)
           rhat_dot_r_new = self.bilinearform(rhat, r)
           beta = rhat_dot_r_new / rhat_dot_r
           rhat+=beta * d
           d=rhat
           rhat_dot_r = rhat_dot_r_new
           if rhat_dot_r<0: print("negative norm.")
           U = - self.getGravity(x)
           U2data = np.array(self.locG(U[2]))
           diffG=U2data-self.gz
           mf=np.inner(diffG, diffG*self.wdsq )
           smooth=self.getSPSP(x)
           print(piter, 'mf',mf,'smooth',smooth, 'rzrz', rhat_dot_r/rzrz0)
           mfs.append(mf)
           smooths.append(smooth)
           rzrzs.append(np.single(rhat_dot_r/rzrz0))
        print(("PCG: tolerance reached after %s steps."%piter))
        smooths=np.array(smooths)
        mfs=np.array(mfs)
        rzrzs=np.array(rzrzs)
        np.savetxt(self.name+'smooths.csv', smooths, delimiter=",")
        np.savetxt(self.name+'mfs.csv', mfs,delimiter=",")
        np.savetxt(self.name+'rzrzs.csv', rzrzs,delimiter=",")
        np.savetxt(self.name+'compg.csv',U2data,delimiter=",")
        np.savetxt(self.name+'diffG.csv',diffG,delimiter=",")
        return x#, smooths, mfs, rzrzs

    def solve(self):
        r = self.RHS()
        if self.verboseLevel=="low":
            P,r,rhatr = PCG(r, self.Aprod, self.P0, self.Msolve, self.bilinearform,
                 atol=self.atol, rtol=self.rtol, iter_max=self.iter_max, 
                 initial_guess=True, verbose=False)
        elif self.verboseLevel=="medium":
            P,r,rhatr = PCG(r, self.Aprod, self.P0, self.Msolve, self.bilinearform,
                 atol=self.atol, rtol=self.rtol, iter_max=self.iter_max, 
                 initial_guess=True, verbose=True)
        elif self.verboseLevel == "high":
            P = self.myPCG(self.P0, r, self.iter_max, self.rtol)
        U = self.getGravity(P)

        pdeG = LinearSinglePDE(self.domain)
        pdeG.setSymmetryOn()
        pdeG.setValue(A = kronecker(3))
        pdeG.setValue(q = self.qtop)
        pdeG.setValue(Y = - self.beta*self.rho_0*P[3])
        optionsG=pdeG.getSolverOptions()
        optionsG.setPackage(SolverOptions.TRILINOS)
        optionsG.setSolverMethod(SolverOptions.PCG)
        optionsG.setPreconditioner(SolverOptions.AMG)

        u = pdeG.getSolution()
        gradu= grad(u, ReducedFunction(dom))

        saveSilo(self.name+"_final", gravity = - U[2], grav2 = - gradu[2], rho=P[3]*self.rho_0)
        print('results silo saved to '+self.name+"_final"+'.silo')
        return P[3]



########################################################################
### Input files and variables from file 
parser = argparse.ArgumentParser(description='Gravity inversion for point data in csv format.', epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')
args = parser.parse_args()
config = importlib.import_module(args.config)
print("Configuration "+args.config+".py imported.")

rho_0    = config.rho_0
atol     = config.atol  
rtol     = config.rtol
pdetol   = config.pdetol
iter_max = config.iter_max 
data_scale = config.data_scale

s = config.s 
a = config.a
b = config.b
mu=1./(8*np.pi*s*a**3)

gz = np.loadtxt(config.gravity_data_file, delimiter=',')*data_scale
acc = np.array(np.loadtxt(config.acc_data_file, delimiter=','))*data_scale
MeasEs = np.loadtxt(config.obsPts_file, delimiter=',')
recorders=[]
for bob in MeasEs:
    recorders.append((bob[0],bob[1],bob[2]))

gz=np.array(gz)
            
measnum = len(gz)
norm_data_sq = np.inner(gz,gz)
print(measnum)


dataWt = config.dataWt
depthWeight = config.depthWeight 
# data weighting
wdsq = 1./(2.*norm_data_sq)*np.ones(measnum)
if dataWt =='relative':
    wdsq=1./(2.*measnum*gz**2)
if dataWt =='accuracy':
    wdsq = np.array(1./(measnum*acc**2))
    
# build domain 
filename, file_extension = os.path.splitext(config.mesh_name)
if file_extension == ".msh":
    dom=ReadGmsh(config.mesh_name, numDim = 3,
                         diracPoints = [sp for sp in recorders], 
                         diracTags = [st for st in range(len(recorders))])
else:
    dom=ReadMesh(config.mesh_name, numDim = 3,
                         diracPoints = [sp for sp in recorders], 
                         diracTags = [st for st in range(len(recorders))])
print("Mesh read from "+config.mesh_name) 

coord=dom.getX()
coord2=ReducedFunction(dom).getX()
depthWt = whereNegative(coord2[2])

minZ =inf(coord[2])
gd = 0.0
xmin = inf(coord[0])
xmax = sup(coord[0])
ymin = inf(coord[1])
ymax = sup(coord[1])
gVol = (xmax-xmin)*(ymax-ymin)*(gd-minZ)


if depthWeight == "coreWt":
    coreD = config.coreD
    factor = 1.+coord2[2]/coreD*wherePositive(coord2[2]-coreD)+(1.-wherePositive(coord2[2]-coreD))
    depthWt = depthWt*(factor)
if depthWeight == "baseWt":
    depthWt = depthWt*coord2[2]/inf(coord2[2])

rho_e = Scalar(rho_0,ReducedFunction(dom))*depthWt

#saveSilo("bobish", rhoref=rho_e)

print("saved bobbish")
P0 = Data(0., (4,), Solution(dom)) 
grav = FOSLSGravity(dom, gz=gz, recorders=recorders, rho_0=rho_e, P0=P0, 
                    wdsq=wdsq, mu=mu, a = a, b = b, atol=atol, rtol=rtol, 
                    iter_max = iter_max, pde_tol=pdetol, name = config.output_name,
                    verboseLevel = config.VerboseLevel)
p = grav.solve()


rho = p* rho_e
rhobar = (1./gVol)*integrate(rho)
sigma = (1./gVol)*integrate((rho-rhobar)**2)
stddev=np.sqrt(sigma)
print('Min density  ', inf(rho))
print('Max density  ', sup(rho))
print('mean density ', rhobar)
print('variance     ', sigma)
print('stddev       ', stddev)  

print("finished")




    
