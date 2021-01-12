#!/usr/bin/python3
import esys.escript.unitsSI as U
import importlib, sys, os
sys.path.insert(0, os.getcwd())
import argparse
parser = argparse.ArgumentParser(description='Inputs needed to run gravity inversion for plane data.',epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')
