#!/bin/bash

. ../version
pw=$QEDIR/bin/pw.x
gipaw=../../bin/gipaw.x

#set -x
mpirun -np 4 $pw <CeO2-scf-U0.in >CeO2-scf-U0.out
mpirun -np 4 $gipaw <CeO2-nmr.in >CeO2-nmr-U0.out

mpirun -np 4 $pw <CeO2-scf-U5.in >CeO2-scf-U5.out
mpirun -np 4 $gipaw <CeO2-nmr.in >CeO2-nmr-U5.out

rm -f *.magres
grep ! *scf*.out* > RESULTS.txt
grep "Total sigma:" *nmr*.out* >>RESULTS.txt

mkdir $VERSION
cp RESULTS.txt $VERSION
diff reference/RESULTS.txt $VERSION/RESULTS.txt
