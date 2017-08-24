#!/bin/bash

# This is an install script for WashUDAQ software. Run this script with the
# desired WashUDAQ version as a parameter, e.g.,
#
# ./install.sh washudaq-2.9-002
#
# The zipped software to install (e.g., washudaq-2.9-002.tar.gz) must be 
# present in the same directory as this script.

vers=$1
tar xzf $vers".tar.gz"
cd $vers
./configure --prefix=/home/wudaq/WashUDAQ/
make all install
