#!/usr/bin/python
__author__ = "morganlnance"

import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("dir", type=str, help="")
input_args = parser.parse_args()

import sys, os

cur_dir = os.getcwd() + '/'
if input_args.dir.endswith( '/' ):
    dir = input_args.dir
else:
    dir = input_args.dir + '/'

files = []
all_files = os.listdir( dir )
for f in all_files:
    if f.endswith( ".pdb" ) or f.endswith( ".fasc" ):
        files.append( f )

for f in files:
    orig = f
    new = f.replace( "_50_", "_25_" )
    
    cmd = "mv %s %s" %( orig, new )
    os.popen( cmd )
