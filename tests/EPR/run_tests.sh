#!/bin/bash

pw=$HOME/Codes/qe-7.4.1/bin/pw.x
gipaw=../../bin/gipaw.x

set -x
mpirun -np 4 $pw <KCl_O2-scf.in >KCl_O2-scf.out_4cpu
mpirun -np 4 $gipaw <KCl_O2-gtensor.in >KCl_O2-gtensor.out_4cpu
mpirun -np 4 $gipaw <KCl_O2-hyperfine.in >KCl_O2-hyperfine.out_4cpu

rm -f *.magres
../../tools/EPR_diagonalize_g-tensor.py KCl_O2-gtensor.out_4cpu >GTENSOR.txt
grep -A18 spin-densities KCl_O2-hyperfine.out_4cpu >HYPERFINE.txt
