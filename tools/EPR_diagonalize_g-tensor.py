#!/usr/bin/env python
# diagonalize the g-tensor from GIPAW output, according to the
# prescription of Weyl, Bolton "Electron Paramagnetic Resonance", Chap. 4.4
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

def diag_tensor(g):
    return np.linalg.eigh( np.dot(g, g.T) )

def print_tensor(g):
    np.savetxt(sys.stdout, g, fmt="%12.2f")


# parse command line arguments
if len(sys.argv) <= 1:
    sys.stderr.write("usage: %s g-tensor.out\n" % (sys.argv[0]))
    sys.exit(1)

# open GIPAW output
try:
    with open(sys.argv[1]) as fp:
        out = fp.readlines()
except:
    sys.stderr.write("error opening file: %s\n" % (sys.argv[1]))
    sys.exit(1)


# search for g-tensor
gparatec, qeq7 = None, None
for i in range(len(out)):
    if out[i].find("Delta_g total (SOO a la Paratec)") >= 0:
        gparatec = read_tensor(out[i+1], out[i+2], out[i+3])
    if out[i].find("Delta_g total (SOO as in Eq.(7))") >= 0:
        geq7 = read_tensor(out[i+1], out[i+2], out[i+3])
    
if gparatec is None or geq7 is None:
    sys.stderr.write("error reading g-tensors from file: %s\n" % (sys.argv[1]))
    sys.exit(1)

print("Delta g-tensors from GIPAW output:")
print("(SOO a la Paratec):")
print_tensor(gparatec)
print()
print("(SOO as in Eq.(7)):")
print_tensor(geq7)
print()

g_e = 2.00231930436182
print("using g_e =", g_e)
print()

np.set_printoptions(precision=6, suppress=True)

# diagonalize
print("Diagonal + principal components (cartesian):")
print("(SOO a la Paratec):")
gparatec = g_e*np.eye(3) + gparatec/1e6
g2, pc = diag_tensor(gparatec)
print("%12.6f =>" % (math.sqrt(g2[0])), pc[:,0])
print("%12.6f =>" % (math.sqrt(g2[1])), pc[:,1])
print("%12.6f =>" % (math.sqrt(g2[2])), pc[:,2])
print()

print("(SOO as in Eq.(7)):")
geq7 = g_e*np.eye(3) + geq7/1e6
g2, pc = diag_tensor(geq7)
print("%12.6f =>" % (math.sqrt(g2[0])), pc[:,0])
print("%12.6f =>" % (math.sqrt(g2[1])), pc[:,1])
print("%12.6f =>" % (math.sqrt(g2[2])), pc[:,2])
print()






