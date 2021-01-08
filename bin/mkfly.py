#!/usr/bin/python3
from esys.finley import ReadGmsh
import sys
sys.path.append(os.getcwd())
'''
converts meshname.msh to meshname.fly

run-escript mkfly.py meshname
'''

fname=sys.argv[1]
mesh_file = fname+'.msh'
fly_file = fname+".fly"
 
print("converting "+mesh_file+" to "+fly_file)
dom = ReadGmsh(mesh_file, numDim=3) 
dom.write(fly_file)
print("finished")

