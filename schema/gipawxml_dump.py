#!/usr/bin/env python
# dump the XML GIPAW results in human-readable format

import sys
from xml.etree import ElementTree as ET
import numpy as np

np.set_printoptions(precision=4, suppress=True)
print_all = False   # set to True to debug all XML results

if len(sys.argv) <= 1:
    sys.stderr.write('usage: {0} prefix-gipaw.xml\n'.format(sys.argv[0]))
    sys.exit(1)

xml_file = sys.argv[1]
try:
    with open(xml_file) as f:
        xml_file_content = f.read()
    root = ET.fromstring(xml_file_content)
except Exception as ex:
    print("Unexpected error opening %s: %s" % (xml_file, ex))
    sys.exit(1)



job = root.find('input').find('job').text
print('GIPAW job:', job)
print()


output = root.find('output')

if job == 'nmr' or job == 'epr' or print_all:
    susc = output.find('susceptibility_low').text
    susc = np.array([float(x) for x in susc.split()], dtype=float).reshape(3,3)
    print('Susceptibility lower bound (10^{-6} cm^3/mol):')
    print(susc)
    print()
    
    susc = output.find('susceptibility_high').text
    susc = np.array([float(x) for x in susc.split()], dtype=float).reshape(3,3)
    print('Susceptibility upper bound (10^{-6} cm^3/mol):')
    print(susc)
    print()
    
if job == 'nmr' or print_all:
    sigma = output.find('shielding_tensors')
    print('NMR shielding tensors (ppm):')
    for atom in list(sigma):
        print('{0:4s} {1:4s} tau={2}'.format(atom.attrib['name'], atom.attrib['index'], atom.attrib['tau']))
        value = np.array([float(x) for x in atom.text.split()], dtype=float).reshape(3,3)
        print(value)
    print()

if job == 'efg' or print_all:
    efg = output.find('electric_field_gradients')
    print('NMR electric field gradients (MHz):')
    for atom in list(efg):
        print('{0:4s} {1:4s} tau={2}'.format(atom.attrib['name'], atom.attrib['index'], atom.attrib['tau']))
        value = np.array([float(x) for x in atom.text.split()], dtype=float).reshape(3,3)
        print(value)
    print()
        
if job == 'epr' or print_all:
    deltag = output.find('delta_g').text
    deltag = np.array([float(x) for x in deltag.split()], dtype=float).reshape(3,3)
    print('EPR delta_g (ppm):')
    print(deltag)
    print()

    deltag = output.find('delta_g_paratec').text
    deltag = np.array([float(x) for x in deltag.split()], dtype=float).reshape(3,3)
    print('EPR delta_g a la paratec (ppm):')
    print(deltag)
    print()

if job == 'hyperfine' or print_all:
    units = root.find('input').find('hfi_output_unit').text
    hf = output.find('hyperfine_dipolar')
    print('EPR hyperfine dipolar tensors ({0}):'.format(units))
    for atom in list(hf):
        print('{0:4s} {1:4s} tau={2}'.format(atom.attrib['name'], atom.attrib['index'], atom.attrib['tau']))
        value = np.array([float(x) for x in atom.text.split()], dtype=float).reshape(3,3)
        print(value)
    print()

    hf = output.find('hyperfine_fermi_contact')
    print('EPR hyperfine Fermi contact ({0}):'.format(units))
    for atom in list(hf):
        print('{0:4s} {1:4s} tau={2}'.format(atom.attrib['name'], atom.attrib['index'], atom.attrib['tau']))
        value = np.array([float(x) for x in atom.text.split()], dtype=float)
        print(value)
    print()


