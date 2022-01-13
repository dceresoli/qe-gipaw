#!/usr/bin/env python
# convert GIPAW XML to GIPAW input

import sys
from xml.etree import ElementTree as ET
import numpy as np

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

inp = root.find('input')
print('&inputgipaw')
for child in list(inp):
    tag, text = child.tag, child.text

    if text == None:
        text = ''

    if tag in ['job', 'prefix', 'tmp_dir', 'restart_mode', 'verbosity',
               'filcurr', 'filfield', 'filnics', 'diagonalization', 'hfi_output_unit']:
        print('    {0} = \'{1}\''.format(tag, text))

    elif tag in ['pawproj', 'use_nmr_macroscopic_shape', 'spline_ps', 'hfi_via_reconstruction_only']:
        text = text.replace('false', '.f.')
        text = text.replace('true', '.t.')
        print('    {0} = {1}'.format(tag, ', '.join(text.split())))

    elif tag in ['q_efg', 'hfi_nuclear_g_factor', 'nmr_macroscopic_shape']:
        print('    {0} = {1}'.format(tag, ', '.join(text.split())))

    else:
        print('    {0} = {1}'.format(tag, text))

print('/')

quit()

# job nmr
# prefix benzene
# tmp_dir ./scratch/
# conv_threshold 1.000000000000e-14
# restart_mode restart
# q_gipaw 3.887614121565e-2
# verbosity high
# filcurr None
# filfield None
# filnics None
# pawproj false false false false false false false false false false
# use_nmr_macroscopic_shape false
# nmr_macroscopic_shape 6.666666666667e-1 0.000000000000e0 0.000000000000e0 0.000000000000e0 6.666666666667e-1 0.000000000000e0 0.000000000000e0 0.000000000000e0 6.666666666667e-1
# spline_ps true
# diagonalization cg
# q_efg 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0
# max_seconds 1.000000000000e7
# r_rand 1.000000014901e-1
# hfi_output_unit MHz
# hfi_nuclear_g_factor 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0 1.000000000000e0
# core_relax_method 1
# hfi_via_reconstruction_only false
#

