#!/bin/bash

# This is an install script for WashUDAQ software
# Run this script from the shell with the 
# desired WashUDAQ version as a parameter (i.e., X.X) 

vers=$1
tar xzf $vers".tar.gz"
cd $vers
./configure --prefix=/home/wudaq/WashUDAQ/
make all install
