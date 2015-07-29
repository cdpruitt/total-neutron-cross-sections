#!/bin/bash
if [ -f ~/WashUDAQ/output/wutest.evt ]
  then
      rm ~/WashUDAQ/output/wutest.evt
fi

bin/readout testConfig/config.txt output/wutest.evt
