#!/usr/bin/python3
__copyright__ = "Copyright (c) 2021 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Andrea Codd"

import importlib, sys, os
sys.path.insert(0, os.getcwd())
import argparse

from esys.escript import *
from esys.finley import ReadGmsh, ReadMesh
#import esys.escript as escript
#import esys.finley as finley
import esys.escript.unitsSI as U
import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import PCG
from esys.downunder import *
from esys.weipa import * 
#from esys.escript import integrate, wherePositive, length
#from esys.escript import saveDataCSV
#from esys.escript.pdetools import MaskFromTag
#import sys
from scipy.io import netcdf_file
#from time import time 

class ACmag(object):
    def __init__(self, domain, w_e, bBx, bBy, bBz, dmag_e, k_0, mu, m0,
         atol=1.0, rtol=1.0, iter_max=100, pde_tol=1e-8, output_name='solutions', verboseLevel="low"):

        self.domain = domain
        self.w_e = w_e
        self.bBv = (bBx*kronecker(3)[0]+bBy*kronecker(3)[1]+bBz*kronecker(3)[2])
        self.magbB = length(self.bBv)
        self.wbBv = -self.w_e*self.bBv/self.magbB
        self.dmag_e = dmag_e
        self.dd   = integrate(self.dmag_e*self.dmag_e)                                 
        self.mu   = mu*self.dd
        self.m0   = m0       
        self.atol = atol
        self.rtol = rtol
        self.iter_max = iter_max
        self.pdetol = pdetol
        self.verboseLevel = verboseLevel
        self.k_0 = k_0                          
        self.zr_zr = []
        self.mfs=[]
        self.smooths=[]
        
        #boundaries       
        z=self.domain.getX()[2]
        top = whereZero(z-sup(z))
        self.pdeu = self.setupPDE()
        self.pdeu.setValue(q=top)
       
    def setupPDE(self):
        pde=LinearSinglePDE(self.domain, isComplex=False)
        pde.setValue(A=kronecker(3))  
        pde.setSymmetryOn()
        options=pde.getSolverOptions()
        options.setTolerance(self.pdetol)
        options.setSolverMethod(SolverOptions.PCG)
        options.setPackage(SolverOptions.TRILINOS)
        options.setPreconditioner(SolverOptions.AMG)
        options.setTrilinosParameter("max levels", 10)  
        options.setTrilinosParameter("problem: symmetric", True)
        options.setTrilinosParameter("reuse: type", "full")
        return pde

    def RHS(self):                  
        self.pdeu.setValue(X=self.dmag_e*self.wbBv,Y=Scalar(0., Solution(self.domain)))
        u = self.pdeu.getSolution()
        return ArithmeticTuple (-self.k_0*inner(self.bBv,grad(u))/self.mu,Vector(0., Solution(self.domain)))

    def Aprod(self,m):
        self.pdeu.setValue(Y=Scalar(0., Solution(self.domain)), X = self.k_0*self.bBv*m)      
        u = self.pdeu.getSolution()         
        self.pdeu.setValue(Y=Scalar(0., Solution(self.domain)), X = self.wbBv*inner(self.wbBv,grad(u))) 
        u = self.pdeu.getSolution()
        return ArithmeticTuple(self.k_0*inner(self.bBv,grad(u))/self.mu ,grad(m))    
        
    def Msolve(self,r):
        self.pdeu.setValue(X=r[1], Y=r[0])
        return self.pdeu.getSolution()        

    def bilinearform(self, m, r):
        return integrate(inner(grad(m), r[1])+ m*r[0])
    
    def compCmag(self,x):
        self.pdeu.setValue(X = self.k_0*x*self.bBv, Y = self.zero)
        u = self.pdeu.getSolution()       
        return inner(self.wbBv,-grad(u))


    def myPCG(self, x,r):
       piter=0
       rhat=self.Msolve(r)
       d = rhat
       rhat_dot_r = self.bilinearform(rhat, r)
       if rhat_dot_r<0: print("negative norm.")
       rzrz0=rhat_dot_r
       norm_r0=np.sqrt(rhat_dot_r)
       atol2=self.rtol*norm_r0
       if atol2<=0:
          print("Non-positive tolarance.")
       mf = 0.5*integrate((self.compCmag(x)-self.dmag_e)**2)/self.dd
       mfold = 0.5
       self.smooths.append(0.)
       self.mfs.append(0.5)
       self.zr_zr.append(1.0)
       donep01=False
       donep05=False
       donep005=False
       donep003=False
       donep008=False
       while not np.sqrt(rhat_dot_r) <= atol2:
           piter += 1
           if piter  >= iter_max: 
               print("maximum number of %s steps reached."%iter_max)
               break 
           q = self.Aprod(d)
           alpha = rhat_dot_r / self.bilinearform(d, q)
           x += alpha * d
           r += q * (-alpha)      
           rhat = self.Msolve(r)
           rhat_dot_r_new = self.bilinearform(rhat, r)
           beta = rhat_dot_r_new / rhat_dot_r
           rhat += beta * d
           d=rhat
           rhat_dot_r = rhat_dot_r_new
           if rhat_dot_r<0: print("negative norm.")
           cmag= self.compCmag(x)
           mf=0.5*integrate((cmag-self.dmag_e)**2)/self.dd
           if mf<0.05 and not donep05:
               print("0.05")
               saveSilo(self.output_name+"p05", cmag=cmag, m=x*self.k_0)
               donep05=True
           if mf<0.01 and not donep01:
               print("0.01")
               saveSilo(self.output_name+"p01", cmag=cmag, m=x*self.k_0)
               donep01=True
           if mf<0.005 and not donep005:
               print("0.005")
               saveSilo(self.output_name+"p005", cmag=cmag, m=x*self.k_0)
               donep005=True
           if mf<0.003 and not donep003:
               print("0.003")
               saveSilo(self.output_name+"p003", cmag=cmag, m=x*self.k_0)
               donep003=True
           if mf<0.008 and not donep008:
               print("0.008")
               saveSilo(self.output_name+"p008", cmag=cmag, m=x*self.k_0)
               donep008=True
           self.mfs.append(mf)
           smooth=integrate(inner(grad(x),grad(x)))
           self.smooths.append(np.double(smooth))
           self.zr_zr.append(np.single(rhat_dot_r/rzrz0))
           if mf>mfold:
               print("disaster disaster disaster disaster disaster disaster disaster") 
           print('mf ', mf)
           print(("PCG: iteration step %s: residual norm = %e"%(piter, np.sqrt(rhat_dot_r))))
       print(("PCG: tolerance reached after %s steps."%piter))
       allsmooths=np.array(self.smooths)
       allmfs=np.array(self.mfs)
       allzr_zr=np.array(self.zr_zr)
       np.savetxt(self.output_name+"smooths.csv",allsmooths, delimiter=",")  
       np.savetxt(self.output_name+"mfs.csv",allmfs,delimiter=",")
       np.savetxt(self.output_name+"zrzr.csv",allzr_zr, delimiter=",")  
       return x,
    

    def solve(self):
        tstart=time()
        r = self.RHS()
        if self.verboseLevel == "low":
            m = PCG(r, self.Aprod, self.m0, self.Msolve, self.bilinearform,
                 atol=self.atol, rtol=self.rtol, iter_max=self.iter_max, 
                 initial_guess=True, verbose=False)
        elif config.VerboseLevel == "medium":
            m = PCG(r, self.Aprod, self.m0, self.Msolve, self.bilinearform,
                 atol=self.atol, rtol=self.rtol, iter_max=self.iter_max, 
                 initial_guess=True, verbose=True)
        elif config.VerboseLevel == "high":#    
            m = self.myPCG(self.m0, r)
        # magnetism
        cmag= self.compCmag(m[0])
        mf0=0.5*integrate((self.dmag_e)**2)/self.dd 
        mf=0.5*integrate((cmag-self.dmag_e)**2)/self.dd
        print('Initial misfit ', mf0 )
        print('  Final misfit ', mf)
        print("    difference ", mf0-mf)
        return m[0], cmag

class DataReader(object):
    EARTH_R=6378.0088*1000.
    toobig=1.e99
    def __init__(self, filename):
        f=netcdf_file(data_file, 'r', mmap=False)
        self.setReferencePoint()
        latitude=None
        longitude=None
        for latvar in "latitude", "lat", "y_range":
            try:
                latitude=f.variables[latvar]
                break
            except KeyError:
                pass
        for lonvar in "longitude", "lon", "x_range":
            try:
                longitude=f.variables[lonvar]
                break
            except KeyError:
                pass
        assert latitude is not None, "no latitude dimension found."
        assert longitude is not None, "no longitude dimension found."
        
        LAT_SOUTH=min(latitude)
        LAT_NORTH=max(latitude)
        print("Latitude range = %s, %s"%(LAT_NORTH, LAT_SOUTH))
        DataCenterLat=(LAT_SOUTH+LAT_NORTH)/2
        GRID_STEP_X=longitude[1]-longitude[0]
        GRID_STEP_Y=latitude[1]-latitude[0]

        # data region:
        DataSpacingX=GRID_STEP_X
        DataSpacingY=GRID_STEP_Y
        #DataSpacingX=GRID_STEP_X/180.*np.pi*np.cos(DataCenterLat/180.*np.pi)*self.EARTH_R
        #DataSpacingY=GRID_STEP_Y/180.*np.pi*self.EARTH_R
        self.DataSpacing=(DataSpacingX, DataSpacingY)
        self.table=None
        for datavar in f.variables.keys():
            if datavar not in [latvar, lonvar]:
                DATA=f.variables[datavar]
                print("Gravity data variable '%s' read from %s."%(datavar,data_file))
                # convert NaN (not-a-number) values to -1000 so the plotting works
                DATA=np.where(np.isnan(DATA[:]), -1000, DATA[:])

                # make sure data is in native format
                if DATA.dtype.byteorder != '=':

                    DATA=np.array(DATA, dtype=np.float32)
                self.table=DATA
        assert self.table is not None, "no gravity data found."
        print("data range : %s - %s"%(self.table.min(), self.table.max()))

    def setReferencePoint(self, Xref=(0., 0.)):
        self.Xref=Xref
        
    def extractData(self, domain):
        X=Function(domain).getX()
        result=interpolateTable(self.table, X[:2], self.Xref, self.DataSpacing, self.toobig)
        return result

########################################################################
### Input files and variables from file 
parser = argparse.ArgumentParser(description='Gravity inversion for plane data in netcdf format.', epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')
args = parser.parse_args()
config = importlib.import_module(args.config)
print("Configuration "+args.config+".py imported.")

mu       = config.mu
k_0      = config.k_0
rho_r    = config.rho_r
atol     = config.atol  
rtol     = config.rtol
pdetol   = config.pdetol
iter_max = config.iter_max 
data_file= config.data_file
data_scale = config.data_scale
Magmag   = config.Magmag
Incl     = config.Incl 
Decl     = config.Decl

bBx = Magmag*np.cos(Incl*np.pi/180.)*np.sin(Decl*np.pi/180.)
bBy = Magmag*np.cos(Incl*np.pi/180.)*np.cos(Decl*np.pi/180.)
bBz = Magmag*np.sin(-Incl*np.pi/180.)

# build domain 
filename, file_extension = os.path.splitext(config.mesh_name)
if file_extension == ".msh":
    dom=ReadGmsh(config.mesh_name, numDim = 3)
else:
    dom=ReadMesh(config.mesh_name, numDim = 3)
print("Mesh read from "+config.mesh_name)

w_e = escript.Scalar(0,escript.Function(dom))        
w_e.setTaggedValue("DataArea",1)
w_e.expand()

ground_e=escript.Scalar(0,escript.Function(dom))    
ground_e.setTaggedValue("Base",1)
ground_e.setTaggedValue("PaddingBase",1)
ground_e.expand()

# read in data and make g
dr=DataReader(config.data_file)
dmag_e = dr.extractData(dom)*wherePositive(w_e)*data_scale   #### convert to T? butdata in nT too 
del dr

m0 = Scalar(0., ContinuousFunction(dom)) 
k_0 = k_0*ground_e

#if config.VerboseLevel == "low":
cmag = ACmag(dom, w_e, bBx, bBy, bBz, dmag_e, k_0, mu, m0, atol, rtol, iter_max, pdetol, config.output_name, config.VerboseLevel)
m, cmb = cmag.solve()
#elif config.VerboseLevel == "medium":
#    cmag = ACmag(dom, dom, w_e, bBx, bBy, bBz, dmag_e, k_0, mu, m0, atol, rtol, iter_max, pdetol, config.output_name, True)
#    m, cmb = cmag.solve()
#elif config.VerboseLevel == "high":
#    cmag = ACmag(dom, w_e, bBx, bBy, bBz, dmag_e, k_0, mu, m0, atol, rtol, iter_max, pdetol, config.output_name)
#    m, cmb = cmag.solve2()

newk=m*k_0
saveSilo(config.output_name+"_final", c_mag = cmb , m=m, data=dmag_e, k=newk)
print('results silo saved to '+config.output_name+'_final.silo')
print("finished")



    
