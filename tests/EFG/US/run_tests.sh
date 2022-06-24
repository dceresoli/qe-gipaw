#!/bin/bash

pw=$HOME/Codes/qe-7.1/bin/pw.x
gipaw=../../../bin/gipaw.x

set -x
mpirun -np 1 $pw <quartz-scf.in >quartz-scf.out_1cpu
mpirun -np 1 $gipaw <quartz-efg.in >quartz-efg.out_1cpu

mpirun -np 4 $pw <quartz-scf.in >quartz-scf.out_4cpu
mpirun -np 4 $gipaw <quartz-efg.in >quartz-efg.out_4cpu

mpirun -np 4 $pw -npool 2 <quartz-scf.in >quartz-scf.out_4cpu_2pool
mpirun -np 4 $gipaw -npool 2 <quartz-efg.in >quartz-efg.out_4cpu_2pool

rm -f *.magres
grep ! *scf.out* > RESULTS.txt
grep "Cq" *efg.out* >>RESULTS.txt
