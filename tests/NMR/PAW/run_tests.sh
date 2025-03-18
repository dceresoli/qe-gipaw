#!/bin/bash

pw=$HOME/Codes/qe-7.4.1/bin/pw.x
gipaw=../../../bin/gipaw.x

set -x
mpirun -np 1 $pw <quartz-scf.in >quartz-scf.out_1cpu
mpirun -np 1 $gipaw <quartz-nmr.in >quartz-nmr.out_1cpu

mpirun -np 4 $pw <quartz-scf.in >quartz-scf.out_4cpu
mpirun -np 4 $gipaw <quartz-nmr.in >quartz-nmr.out_4cpu

mpirun -np 4 $pw -npool 2 <quartz-scf.in >quartz-scf.out_4cpu_2pool
mpirun -np 4 $gipaw -npool 2 <quartz-nmr.in >quartz-nmr.out_4cpu_2pool

rm -f *.magres
grep ! *scf.out* > RESULTS.txt
grep "Total sigma:" *nmr.out* >>RESULTS.txt
