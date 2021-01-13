#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 16:41:15 2020

@author: andrea
"""

import numpy as np
import argparse
parser = argparse.ArgumentParser(description='this generates gmsh geo file for synthetic application', epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')

args = parser.parse_args()
config = importlib.import_module(args.config)

geoname = config.geo_name

# import measuring locations
pts = np.loadtxt(config.obsPts_file, delimiter=',')
# import nearest neighbour
nn = np.loadtxt(config.minDist_file, delimiter=',')

print('making geo for observatio points ',config.obsPts_file, ' and nearest neighbour file ', config.minDist_file)

# scale for element size based on nearest neighbour
nds=[]
for ind1 in range(len(pts)):
    nndist=nn[ind1]
    el = nndist / config.spacing0
    if nndist < config.mindist1:
        el = nndist/config.spacing1
    if nndist[ind1] < mindist2: 
        el= nndist/config.spacing2
    nds.append(el)
    print(nndist, el)

1/0 

numPts= len(pts)

# core
xmin = 1.01*min(pts[:,0])
xmax = 1.01*max(pts[:,0])
xspan = xmax-xmin

ymin = 1.05*min(pts[:,1])
ymax = 1.05*max(pts[:,1])
yspan = ymax-ymin
span=min(xspan,yspan)

zmin = -20000.
zmax = 10000.
grnd = 0.

# buffer
Bxmin = xmin - span/2
Bxmax = xmax + span/2
Bymin = ymin - span/2
Bymax = ymax + span/2
Bzmin = 2.*zmin
Bzmax = 2.*zmax

# meshsizes
MCtop = 3000
MCbase = 2000
MCground = 2000
MBtop = 20000
MBbase = 20000
MBground = 10000

#
out ="Cxmin=%e;\n"%xmin
out+="Cxmax=%e;\n"%xmax
out+="Cymin=%e;\n"%ymin
out+="Cymax=%e;\n"%ymax
out+="Czmin=%e;\n"%zmin
out+="Czmax=%e;\n"%zmax
out+="Bxmin=%e;\n"%Bxmin
out+="Bxmax=%e;\n"%Bxmax
out+="Bymin=%e;\n"%Bymin
out+="Bymax=%e;\n"%Bymax
out+="Bzmin=%e;\n"%Bzmin
out+="Bzmax=%e;\n"%Bzmax
out+="MCtop=%e;\n"%MCtop
out+="MCbase=%e;\n"%MCbase
out+="MCground=%e;\n"%MCground
out+="MBtop=%e;\n"%MBtop
out+="MBbase=%e;\n"%MBbase
out+="MBground=%e;\n"%MBground
out+="numPts=%e;\n"%numPts

Mpoints={}
text = ""


out+="""
// CORE
// core points
Point(1)  = {Cxmin,Cymin,Czmin,MCbase};
Point(2)  = {Cxmin,Cymin,0.,MCground};
Point(3)  = {Cxmin,Cymin,Czmax,MCtop};
Point(4)  = {Cxmax,Cymin,Czmin,MCbase};
Point(5)  = {Cxmax,Cymin,0.,MCground};
Point(6)  = {Cxmax,Cymin,Czmax,MCtop};
Point(7)  = {Cxmin,Cymax,Czmin,MCbase};
Point(8)  = {Cxmin,Cymax,0.,MCground};
Point(9)  = {Cxmin,Cymax,Czmax,MCtop};
Point(10) = {Cxmax,Cymax,Czmin,MCbase};
Point(11) = {Cxmax,Cymax,0.,MCground};
Point(12) = {Cxmax,Cymax,Czmax,MCtop};

// core lines, base, ground, top, sides
Line(1) = {1,4};
Line(2) = {7,10};
Line(3) = {1,7};
Line(4) = {4,10};
Line(5) = {2,5};
Line(6) = {8,11};
Line(7) = {2,8};
Line(8) = {5,11};
Line(9) = {3,6};
Line(10) = {9,12};
Line(11) = {3,9};
Line(12) = {6,12};
Line(13) = {1,2};
Line(14) = {2,3};
Line(15) = {4,5};
Line(16) = {5,6};
Line(17) = {7,8};
Line(18) = {8,9};
Line(19) = {10,11};
Line(20) = {11,12};

// core horizontal surfaces
Line Loop(1) = {-1,3,2,-4};
Plane Surface(1)={1};
Physical Surface(\"coreBase\") = {1};
Line Loop(2) = {-5,7,6,-8};
Plane Surface(2)={2};
Physical Surface(\"coreInterface\") = {2};
Line Loop(3) = {9,12,-10,-11};
Plane Surface(3)={3};
Physical Surface(\"coreTop\") = {3};
// core front surfaces
Line Loop(4) = {1,15,-5,-13};
Plane Surface(4)={4};
Line Loop(5) = {5,16,-9,-14};
Plane Surface(5)={5};
// core back surfaces
Line Loop(6) = {-2,17,6,-19};
Plane Surface(6)={6};
Line Loop(7) = {-6,18,10,-20};
Plane Surface(7)={7};
// core left surfaces
Line Loop(8) = {-3,13,7,-17};
Plane Surface(8)={8};
Line Loop(9) = {-7,14,11,-18};
Plane Surface(9)={9};
// core right surfaces
Line Loop(10) = {4,19,-8,-15};
Plane Surface(10)={10};
Line Loop(11) = {8,20,-12,-16};
Plane Surface(11)={11};

// BUFFER
// buffer points
Point(21) = {Bxmin,Bymin,Bzmin,MBbase};
Point(22) = {Bxmin,Bymin,0.,MBground};
Point(23) = {Bxmin,Bymin,Bzmax,MBtop};
Point(24) = {Bxmax,Bymin,Bzmin,MBbase};
Point(25) = {Bxmax,Bymin,0.,MBground};
Point(26) = {Bxmax,Bymin,Bzmax,MBtop};
Point(27) = {Bxmin,Bymax,Bzmin,MBbase};
Point(28) = {Bxmin,Bymax,0.,MBground};
Point(29) = {Bxmin,Bymax,Bzmax,MBtop};
Point(30) = {Bxmax,Bymax,Bzmin,MBbase};
Point(31) = {Bxmax,Bymax,0.,MBground};
Point(32) = {Bxmax,Bymax,Bzmax,MBtop};

// buffer lines, base, ground, top, sides
Line(21) = {21,24};
Line(22) = {27,30};
Line(23) = {21,27};
Line(24) = {24,30};
Line(25) = {22,25};
Line(26) = {28,31};
Line(27) = {22,28};
Line(28) = {25,31};
Line(29) = {23,26};
Line(30) = {29,32};
Line(31) = {23,29};
Line(32) = {26,32};
Line(33) = {21,22};
Line(34) = {22,23};
Line(35) = {24,25};
Line(36) = {25,26};
Line(37) = {27,28};
Line(38) = {28,29};
Line(39) = {30,31};
Line(40) = {31,32};

// buffer horizontal surfaces
Line Loop(21) = {-21,23,22,-24};
Plane Surface(21)={21};
Physical Surface(\"bufferBase\") = {21};
Line Loop(22) = {-25,27,26,-28};
Plane Surface(22)={22,-2};
Physical Surface(\"bufferInterface\") = {22};
Line Loop(23) = {29,32,-30,-31};
Plane Surface(23)={23};
Physical Surface(\"bufferTop\") = {23};
// buffer front surfaces
Line Loop(24) = {21,35,-25,-33};
Plane Surface(24)={24};
Line Loop(25) = {25,36,-29,-34};
Plane Surface(25)={25};
// buffer back surfaces
Line Loop(26) = {-22,37,26,-39};
Plane Surface(26)={26};
Line Loop(27) = {-26,38,30,-40};
Plane Surface(27)={27};
// buffer left surfaces
Line Loop(28) = {-23,33,27,-37};
Plane Surface(28)={28};
Line Loop(29) = {-27,34,31,-38};
Plane Surface(29)={29};
// buffer right surfaces
Line Loop(30) = {24,39,-28,-35};
Plane Surface(30)={30};
Line Loop(31) = {28,40,-32,-36};
Plane Surface(31)={31};


//  SENSOR LOCATIONS
"""
for i in range(len(pts)):
    pt=pts[i]
    nd = nds[i]
    out+="k=newp;\n"
    out+="Point(k)={%e,%e,%e,%e};\n"%(pt[0],pt[1],pt[2],nd/4.)
    out+="Point{k} In Surface{2};\n"
out+="""
// VOLUMES
// core volumes
Surface Loop(12) = {1,-2,4,6,8,10};
Volume(1) = {12}; 
Physical Volume(\"coreGround\") = {1};
Surface Loop(13) = {2,5,7,9,11,3};
Volume(2) = {13}; 
Physical Volume(\"coreAir\") = {2};
// buffer volumes
Surface Loop(14) = {21,22,24,26,28,30,-1,-4,-6,-8,-10};
Volume(3) = {14};
Physical Volume(\"bufferGround\") = {3};
Surface Loop(15) = {22,23,25,27,29,31,-5,-7,-9,-11,-3};
Volume(4) = {15}; 
Physical Volume(\"bufferAir\") = {4};
"""
open(geoname,'w').write(out)
print("mesh made and written to file "+geoname)















