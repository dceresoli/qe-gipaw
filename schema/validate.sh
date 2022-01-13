#!/bin/bash

# validate prefix.xml and prefix-gipaw.xml
if [ -z $1 ]; then
   echo usage: $0 prefix
   exit 1
fi
prefix=$1

# validate PW XML file
xmllint --noout --schema qes_211101.xsd $ ${prefix}.xml

# validate GIPAW XML file
xmllint --noout --schema gipaw_xmlschema_import_local.xsd  ${prefix}-gipaw.xml
