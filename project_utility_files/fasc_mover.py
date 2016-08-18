#!/usr/bin/python

# This is to be ran in the directory where all of the fasc files live if you were to move into each pdb directory. It will move the .fasc files into the specified directory

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Use Python to move fasc files from all directories in a base directory into a given directory")
parser.add_argument("fasc_dir", type=str, help="the name of the directory where you want to move the fasc files")
parser.add_argument("-base_dir", type=str, default=None, help="the base directory where, in each directory contained here, is a fasc file. Default = current directory")
input_args = parser.parse_args()

# get current directory
cur_dir = os.getcwd()

# if nothing was passed for base_dir, use current working directory
if input_args.base_dir is None:
    input_args.base_dir = os.getcwd()

# check the validity of the directories
if not os.path.isdir( input_args.base_dir ):
    print "Something is wrong with your base directory. Check your input. Exiting."
    sys.exit()
else:
    base_dir = input_args.base_dir
if not os.path.isdir( input_args.fasc_dir ):
    print "Something is wrong with your fasc directory. I can't make directories. Exiting."
    sys.exit()
else:
    fasc_dir = input_args.fasc_dir

# get the absolute paths to everything
os.chdir( base_dir )
base_dir = os.path.abspath( os.getcwd() ) + '/'
os.chdir( cur_dir )

os.chdir( fasc_dir )
fasc_dir = os.path.abspath( os.getcwd() ) + '/'
os.chdir( cur_dir )

# get the single name of the fasc dir to ignore in case its in the base dir
fasc_dir_name = fasc_dir.split( '/' )[-2]

# move into the base dir and list the files and directories within
os.chdir( base_dir )
dirs = os.listdir( base_dir )

# move into each directory and copy the fasc file to its new location
for d in dirs:
    if d != fasc_dir_name:
        # if this name is a directory
        if os.path.isdir( d ):
        # make an appropriate path to this directory
            working_dir = base_dir + d + '/'
            os.chdir( working_dir )
            # list the files within
            files = os.listdir( working_dir )
            # if the file ends with .fasc
            for f in files:
                if f.endswith( ".fasc" ):
                    # move it to the new location
                    f_name = working_dir + f
                    cmd = "cp %s %s" %(f_name, fasc_dir)
                    os.popen( cmd )
    os.chdir( base_dir )
