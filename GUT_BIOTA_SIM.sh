
#!/bin/bash

resC0=0
resC1=6
hBin=1
numsims=10

c++ -std=c++11 GUT_BIOTA_SIM.cpp -o exe_GUT_BIOTA_SIM -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas

foldername=SIM_RESULTS
mkdir $foldername
saveDirectory="./SIM_RESULTS/"

./exe_GUT_BIOTA_SIM $saveDirectory $numsims $resC0 $resC1 $hBin
