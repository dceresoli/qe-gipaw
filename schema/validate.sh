#!/bin/bash

# validate prefix.xml and prefix-gipaw.xml
if [ -z $1 ]; then
   echo usage: $0 prefix
   exit 1
fi
prefix=$1

schemadir=$HOME/Codes/qe-gipaw/schema

# validate PW XML file
xmllint --noout --schema $schemadir/qes_211101.xsd ${prefix}.xml

# validate GIPAW XML file
xmllint --noout --schema $schemadir/gipaw_xmlschema_import_local.xsd  ${prefix}-gipaw.xml
