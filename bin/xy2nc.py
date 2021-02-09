#!/usr/bin/python3
#import esys.escript.unitsSI as U
import numpy as np
import importlib, sys, os
from converters import writeNetCDF
sys.path.insert(0, os.getcwd())
import pathlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import argparse
parser = argparse.ArgumentParser(description='this converts an xy file to a NetCDF file for gravity and magnetic anomaly data', epilog="version 01/2021 by l.gross@uq.edu.au")
parser.add_argument(dest='xyfile', metavar='XYFILE', type=str, help='xy file.')
parser.add_argument('--out', '-o', dest='outfile', metavar='OUTFILE', type=str, default=None , help='name of the .')
parser.add_argument('--units', '-u', dest='units', metavar="UNITS", type=str, default=None, help = "units of the grav or mag data")
parser.add_argument('--mag', '-m', dest='magnetic',  action='store_true', default=False, help='these are magnetic data')
parser.add_argument('--name', '-n', dest='name',  metavar="NAME", default=None, help='name of the data field')
parser.add_argument('--plot', '-p', dest='plot',  type=str, default=None, help='name of plot file. if not set no plot is generated.')


args = parser.parse_args()

print("reading from "+args.xyfile)
title=pathlib.Path(args.xyfile).stem
basename=os.path.splitext(args.xyfile)[0]
if not args.outfile is None: 
    outfile=args.outfile
else:
    outfile=basename+".nc"

x,y,data=np.loadtxt(args.xyfile, comments='#', unpack=True)
print(len(x)," data point found.")
xN=np.nonzero( abs(x[1:] -x[0]) < 1e-8 * np.linalg.norm(x, ord=np.inf))[0][0]
yN=np.nonzero( abs(y[1:] -y[0]) < 1e-8 * np.linalg.norm(y, ord=np.inf))[0][0]
print(xN, yN)
print(x[0], x[xN-1], x[xN],x[xN+1],)
print(y[0], y[xN-1], y[xN],y[xN+1],)
if yN > 0:
    yy=x[:yN+1]
    xx=y[::yN+1]
    dd=data.reshape((len(xx), len(yy))).T
else:
    xx=x[:xN+1]
    yy=y[::xN+1]
    dd=data.reshape((len(yy), len(xx)))

print("x is [%s, %s]"%(xx[0], xx[-1]))
print("y is [%s, %s]"%(yy[0], yy[-1]))
print("grid shape is %s x %s"%(dd.shape))
print("data range %s - %s (mean =%s)"%(dd.min(), dd.max(), dd.mean()))
if args.units:
    dataunits=args.units
else:
    if args.magnetic:
        dataunits="nT"
    else:
        dataunits='mgal'


if args.name:
    name=args.name
    longname=args.name
else:
    if args.magnetic:
        name="magnetic"
        longname='magnetic_field_anomaly'
    else:
        name='gravity'
        longname='vertical_gravity',
        
n=writeNetCDF(filename=outfile, 
                data=dd,
                units='m',
                units_data=dataunits,
                name=name,
                delta=(yy[1]-yy[0], xx[1]-xx[0]),
                origin=(yy[0], xx[0]),
                longname=longname,
                error=None,
                title=title)

if args.plot:
    plt.figure()
    plt.imshow(dd,  origin='lower', extent=(xx[0]/1000, xx[-1]/1000, yy[0]/1000, yy[-1]/1000),  cmap=cm.rainbow)
    plt.colorbar()
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")    
    plt.title(title)
    plt.savefig(args.plot)
    print(f"{args.plot} image generated.")
    
print(f"NetCDF file is {outfile}")

