#!/bin/bash

#Este programa sirve para automatizar el analisis de las corridas sobre el studio de strain en NW de AuPd

mkdir scripts
mkdir data
for c in 0 1 2 3 4 5; do 
    sed 's|n=0|n='$c'| ; s|wf.dat|../data/wf'$c'.dat| ' MMB.f90 > scripts/E$c.f90
    #gfortran scripts/E$c.f90 -o scripts/E$c
    #./scripts/E$c.exe
done