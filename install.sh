#!/bin/bash
# Installs washudaqX.X (input as a parameter)
vers=$1
tar xzf $vers".tar.gz"
cd $vers
./configure --prefix=/home/wudaq/WashUDAQ/
make all install
