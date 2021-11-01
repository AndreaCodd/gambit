#!/usr/bin/python3
#import esys.escript.unitsSI as U
import numpy as np
import importlib, sys, os
sys.path.insert(0, os.getcwd())
import json
import argparse
parser = argparse.ArgumentParser(description='this generates gmsh geo file to create for synthetic data', epilog="version 11/2020 by l.gross@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')
parser.add_argument('--geo', '-g', dest='geofile', metavar='GEO', type=str, default=None , help='name of the GMSH geo file. if not set project name is use.')
#parser.add_argument('--geo', '-g', dest='geo', metavar='geofile', type=str, default=None, help='gmsh geo file. if not given [project].geo is used.')
#parser.add_argument(dest='project', metavar='project', type=str, help='name of the project to be created')
#parser.add_argument('--novtx', '-nv', dest='ignorevtx', action='store_true', help="if set vtx points are ignored.")



print("** Generates a geo file for synthetic gravity and/or magnetic data file generation **")
args = parser.parse_args()

config = importlib.import_module(args.config)
print("configuration "+args.config+".py imported.")
if config.DataHeightAboveGround >0:
    GEOTEMPLATE=os.path.join(os.path.dirname(os.path.abspath(__file__)), "AboveGroundTemplate.geo")
    print('above ground template used.')
else:
    GEOTEMPLATE=os.path.join(os.path.dirname(os.path.abspath(__file__)), "OnGroundTemplate.geo")
    print('on ground template used.')    
DataSpacingX=config.LDataX/(config.DataNumX-1)
DataSpacingY=config.LDataY/(config.DataNumY-1)

mapping= { 
"DataRefX" : config.DataRefX-DataSpacingX, 
"DataRefY" : config.DataRefY-DataSpacingY,
"DataHeightAboveGround" : config.DataHeightAboveGround,
"DataSpacingX" : DataSpacingX,
"DataSpacingY" : DataSpacingY,
"DataNumX" : config.DataNumX,+2,
"DataNumY" : config.DataNumY,+2,
"DataMeshSizeVertical" : config.DataMeshSizeVertical,
"CoreThickness" : config.CoreThickness,
"AirLayerThickness" : config.AirLayerThickness,
"PaddingX" : config.PaddingX,
"PaddingY" : config.PaddingY,
"PaddingZ" : config.PaddingZ,
"PaddingAir" : config.PaddingAir,
"MeshSizeAirFactor" : config.MeshSizeAirFactor,
"MeshSizeCoreFactor" : config.MeshSizeCoreFactor,
"MeshSizePaddingFactor" : config.MeshSizePaddingFactor
}
print("Setting: ")
print(json.dumps(mapping, sort_keys=True, indent=4)[2:-2])

if args.geofile:
    GEOFN="syth_"+args.geofile
else:
    GEOFN="syth_"+config.project+".geo"
text=open(GEOTEMPLATE,'r').read().format(**mapping)
open(GEOFN, "w").write(text)
print("GMSH geofile has been written to ",GEOFN)
print("to generate mesh run:")
print("     gmsh -3 -format msh2 -o %s %s"%("syth_"+config.meshfile, GEOFN))


