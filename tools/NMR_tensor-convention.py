#!/usr/bin/env python
# convert NMR shielding tensors from GIPAW output into different convetions
# Davide Ceresoli, Jan 2017

import numpy as np
import math
import sys

def read_tensor(l1, l2, l3):
    tens = []
    tens.append( list(map(float, l1.split())) )
    tens.append( list(map(float, l2.split())) )
    tens.append( list(map(float, l3.split())) )
    return np.array(tens)

def haerbelen(d):
     iso = (d[0] + d[1] + d[2]) / 3.0
     d = d[ np.argsort(np.abs(d - iso)) ]
     aniso = d[2] - iso
     if abs(aniso) > 1e-6:
	eta = (d[1]-d[0])/aniso
     else:
        eta = 0.0
     aniso *= 1.5
     return iso, aniso, eta

def herzfeld_berger(d):
     iso = (d[0] + d[1] + d[2]) / 3.0
     d = d[ np.argsort(d) ]
     span = d[2] - d[0]
     if span > 1e-6:
	skew = 3*(d[1]-iso)/span
     else:
        skew = 0.0
     return iso, span, skew


# parse command line arguments
if len(sys.argv) <= 1:
    sys.stderr.write("usage: %s nmr.out\n" % (sys.argv[0]))
    sys.exit(1)

# open GIPAW output
try:
    with open(sys.argv[1]) as fp:
        out = fp.readlines()
except:
    sys.stderr.write("error opening file: %s\n" % (sys.argv[1]))
    sys.exit(1)


# search for NMR shielding tensors
sigma = []
atoms = []
pos = []
for i in range(len(out)):
    if out[i].find("Total NMR chemical shifts in ppm:") >= 0:
        i += 3
        while out[i].find("Atom") >= 0:
            atoms.append(out[i].split()[2])
            out[i] = out[i].replace(")", " ")
            pos.append(list(map(float, out[i].split()[5:8])))
	    sigma.append(read_tensor(out[i+1], out[i+2], out[i+3]))
            i += 10

natoms = len(atoms)    

np.set_printoptions(precision = 6, suppress = True)

print("="*72)
print("Principal components and directions:")
print("="*72)
for i in range(natoms):
    sigma[i] = 0.5*(sigma[i] + sigma[i].T)  # symmetric part only
    d, v = np.linalg.eigh(sigma[i])
    print("%-2s %4i   sigma: %10.4f =>" % (atoms[i], i+1, d[0]), v[:,0])
    print("          sigma: %10.4f =>" % (d[1]), v[:,1])
    print("          sigma: %10.4f =>" % (d[2]), v[:,2])
print()


print("="*72)
print("HAERBELEN (SIMPSON) convention:")
print("="*72)
for i in range(natoms):
    d = np.linalg.eigvalsh(sigma[i])
    iso, aniso, eta = haerbelen(d)
    print("%-2s %4i   iso: %10.4f   aniso: %10.4f    eta:  %10.4f" % (atoms[i], i+1, iso, aniso, eta))
print()

print("="*72)
print("HERZFELD-BERGER convention:")
print("="*72)
for i in range(natoms):
    d = np.linalg.eigvalsh(sigma[i])
    iso, span, skew = herzfeld_berger(d)
    print("%-2s %4i   iso: %10.4f   span:  %10.4f    skew: %10.4f" % (atoms[i], i+1, iso, span, skew))

quit()
quit()
