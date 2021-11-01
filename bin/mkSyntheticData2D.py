#!/usr/bin/python3
#import esys.escript.unitsSI as U
from esys.escript import *
from esys.finley import ReadMesh, ReadGmsh
from esys.weipa import saveSilo
from esys.downunder.apps import GravityModel, MagneticModel3D
from converters import writeNetCDF, grepValuesByMask
from esys.escript.pdetools import Locator
import matplotlib.pyplot as plt
    
import numpy as np
import importlib, sys, os
sys.path.insert(0, os.getcwd())
import json
import argparse
parser = argparse.ArgumentParser(description='this generates gmsh geo file for synthetic application', epilog="version 11/2020 by l.gross@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')
parser.add_argument('--silo', '-s', dest='silo',  type=str, default="output", help='name silo file to show resuls')
parser.add_argument('--mag', '-m', dest='magnetic',  action='store_true', default=False, help='create synthetic magnetic data')
parser.add_argument('--grav', '-g', dest='gravity',  action='store_true', default=False, help='create synthetic gravity data')
parser.add_argument('--debug', '-d', dest='debug',  action='store_true', default=False, help='print a few more informations')
ToMgal=1e-5

print("** Generates synthetic gravity and/or magnetic data **")
args = parser.parse_args()

config = importlib.import_module(args.config)
print("configuration "+args.config+".py imported.")

filename, file_extension = os.path.splitext(config.meshfile)
if file_extension == ".msh":
    domain=ReadGmsh(config.meshfile, 3, optimize=True )
else:
    domain=ReadMesh(config.meshfile,  optimize=True)
print("mesh read from "+config.meshfile)

def createStructure(x, structure, show=False):
    """
    this create a structure  out a
    """
    out=Scalar(0., x.getFunctionSpace())
    for v, s in structure:
        m=Scalar(0., out.getFunctionSpace())
        for ssub in s:
            d=length(x-[ ssub['xc'],ssub['yc'], ssub['zc']  ])
            m+=whereNegative(d-ssub['r'])
            if show:
                print("spherical anomaly of %s at (%s, %s, %s) m with radius %s m added."%(v,ssub['xc'],ssub['yc'], ssub['zc'], ssub['r']))
        out+=v*wherePositive(m)
    if show:
        print("anomaly = ",out)
    return out

DataHeight=config.DataHeightAboveGround+config.DataMeshSizeVertical/2
#config.DataNumX
#config.DataNumY
DataX=np.linspace(config.DataRefX, config.DataRefX+config.LDataX, config.DataNumX, endpoint=True)
DataY=np.linspace(config.DataRefY, config.DataRefY+config.LDataY, config.DataNumY, endpoint=True)
#DataZ=np.linspace(config.DataMeshSizeVertical/2+config.DataHeightAboveGround, config.DataMeshSizeVertical/2+config.DataHeightAboveGround, 1, endpoint=True)
#xi=np.meshgrid(DataX, DataY, DataZ)
xi=np.meshgrid(DataX, DataY)


data_mask=Scalar(0.,ReducedFunction(domain))
data_mask.setTaggedValue("DataArea", 1)


#locator=Locator( x=[ [x_cord[ix,iy], y_cord[ix,iy],DataHeight] for ix,iy in np.ndindex(y_cord.shape)], where=ReducedFunction(domain))
#print(locator.getX())


output={"tags" : makeTagMap(Function(domain))}
if args.gravity:
    # Define model using class GravityModel
    # Assumes zero Dirichlet BC at top surface
    model=GravityModel(domain, fixBase=True)

    
    model.setDensity(createStructure(ReducedFunction(domain).getX(), config.true_density, show=True ))

    # Solve for gravity field anomaly
    potential=model.getGravityPotential()

    output['g_potential'] = potential
    print("gravity potential= ",potential )
    output['g'] = model.getGravityVector()/ToMgal
    output['gz'] = model.getzGravity()/ToMgal
    print("vertical gravity = ",output['gz'] )
    output['rho']=model.getDensity()
    
    G, xi =grepValuesByMask(xi, output['gz'] , data_mask)
    print("data region range = [%s, %s] x [%s, %s]"%(xi[0].min(), xi[0].max(), xi[1].min(), xi[1].max() ))
    print("gravity data range = [%s, %s] mgal"%(G.min(), G.max()))
    
    if config.noise>0:
        std=config.noise/100.
        pert=np.random.normal(0.0, scale=std, size=(config.DataNumY,config.DataNumX))
        #print(pert)
        #print(G)
        #print(pert.shape)
        #print(G.shape)
        G*=(1+pert)    
        print("%g %% noise added."%config.noise)
    else:
        print("no noise added.")
    n=writeNetCDF(filename=config.gravfile, 
                data=G,
                units='m',
                units_data='mgal',
                name='gravity',
                delta=(DataY[1]-DataY[0], DataX[1]-DataX[0]),
                origin=(DataY[0], DataX[0]),
                longname='vertical_gravity',
                error=None,
                title=config.project+' gravity data')
    print(f"gravity data written to file {n}")
    plt.figure()
    plt.imshow(G,  origin='lower', extent=(xi[0].min(), xi[0].max(), xi[1].min(), xi[1].max()))
    plt.colorbar()
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")    
    plt.title("Gravity anomaly at height %s "%(config.DataMeshSizeVertical/2+config.DataHeightAboveGround))
    plt.savefig("gravdataNoise.png")
    print("gravity image generated.")
    print("** gravity completed.")
    
if args.magnetic:
    model=MagneticModel3D(domain)

    model.setSusceptibility(createStructure(ReducedFunction(domain).getX(), config.true_magnetization, show=args.debug ))
    
    B_h=np.array([-config.B_hx, -config.B_hy, -config.B_hz])
    print("back ground magnetic field", B_h)
    model.setBackgroundMagneticField(B_h)
    
    # Solve for magnetic field anomaly
    B_a=model.getMagneticFieldAnomaly()
    print("magnetic field anomaly =", B_a)
    b_a=inner(B_a,B_h)/length(B_h)
    print("magnetic intensity anomaly =", b_a)

    output['B_a'] = B_a
    output['BI_a'] = b_a
    output['susceptibility']=model.getSusceptibility()

    B_A, xi =grepValuesByMask(xi, output['BI_a'] , data_mask)
    print("data region range = [%s, %s] x [%s, %s]"%(xi[0].min(), xi[0].max(), xi[1].min(), xi[1].max() ))
    print("magnetic field anomaly data range = [%s, %s] nT"%(B_A.min(), B_A.max()))
    
    if config.noise>0:
        std=config.noise/100.
        pert=np.random.normal(0.0, scale=std, size=(config.DataNumY,config.DataNumX))
        #print(pert)
        #print(B_A)
        #print(pert.shape)
        #print(B_A.shape)
        B_A*=(1+pert)    
        print("%g %% noise added."%config.noise)
        
    n=writeNetCDF(filename=config.magfile, 
                data=B_A,
                units='m',
                units_data='nT',
                name='magnetic',
                delta=(DataY[1]-DataY[0], DataX[1]-DataX[0]),
                origin=(DataY[0], DataX[0]),
                longname='magnetic_field_anomaly',
                error=None,
                title=config.project+' magnetic field anomaly data')
    print(f"magnetic data written to file {n}")
    plt.figure()
    plt.imshow(B_A,  origin='lower', extent=(xi[0].min(), xi[0].max(), xi[1].min(), xi[1].max()))
    plt.colorbar()
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")    
    plt.title("Magnetic field anomaly at height %s "%(config.DataMeshSizeVertical/2+config.DataHeightAboveGround))
    plt.savefig("magdataNoise.png")
    print("magnetic image generated.")
    print("** magnetic completed.")
    
if len(output)>1 and args.silo:
    saveSilo(args.silo, **output)
    print("finally results were written to %s.silo. All done :-)"%args.silo)
else:
    print("nothing to do.")
