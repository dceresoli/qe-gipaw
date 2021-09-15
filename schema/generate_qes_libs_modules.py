#! /usr/bin/env python3

import os, sys
sys.path.append('/home/ceresoli/Codes/XML_GIPAW/xmltool')
from helpers.file_generator import generate_file

input_xsd = "gipaw_xmlschema.xsd"
template_files = [
"templates/init/qes_init_module.f90.j2",
"templates/reset/qes_reset_module.f90.j2",
"templates/bcast/qes_bcast_module.f90.j2",
"templates/read/qes_read_module.f90.j2",
"templates/write/qes_write_module.f90.j2",
"templates/libs/qes_libs_module.f90.j2",
"templates/types/qes_types_module.f90.j2"
]

print("input xsd: "+input_xsd)
print("generate files:")
for template_file in template_files:
    output_file = os.path.splitext(os.path.basename(template_file))[0]
    output_file_path = os.path.join("..", "src", output_file )
    generate_file(input_xsd, template_file, output_file_path)
    print("  "+output_file_path)
