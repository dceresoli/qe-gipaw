#!/bin/bash

PW=$HOME/Codes/qe-7.5/bin/pw.x
GIPAW=$HOME/Codes/qe-gipaw.git/bin/gipaw.x

tmpl=AlB2-scf
for magn in  0.002 0.004 0.006 0.008 0.010
do 
  mpirun -np 16 $PW -npool 16 <$tmpl-$magn.in >$tmpl-$magn.out
  mpirun -np 16 $GIPAW -npool 16 <hyperfine.in >hyperfine-$magn.out
done

