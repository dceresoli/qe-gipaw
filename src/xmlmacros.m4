dnl #define _NE(x) call XML_NewElement(xmlf, #x)
define(_NE,`call XML_NewElement(xmlf, "$1")')

dnl #define _ADDS(x) call XML_AddCharacters(xmlf, trim((x)))
define(_ADDS,`call XML_AddCharacters(xmlf, trim($1))')

dnl #define _ADDV(x) call XML_AddCharacters(xmlf, (x))
define(_ADDV,`call XML_AddCharacters(xmlf, $1)')

dnl #define _EE(x) call XML_EndElement(xmlf, #x)
define(_EE,`call XML_EndElement(xmlf, "$1")')

dnl #define _OUTS(x) call XML_NewElement(xmlf, #x); call XML_AddCharacters(xmlf, trim((x))); call XML_EndElement(xmlf, #x);
define(_OUTS,`call XML_NewElement(xmlf, "$1")
	call XML_AddCharacters(xmlf, trim($1))
	call XML_EndElement(xmlf, "$1")')

dnl #define _OUTV(x) call XML_NewElement(xmlf, #x); call XML_AddCharacters(xmlf, (x)); call XML_EndElement(xmlf, #x);
define(_OUTV,`call XML_NewElement(xmlf, "$1")
	call XML_AddCharacters(xmlf, $1)
	call XML_EndElement(xmlf, "$1")')

dnl #define _ATTR(x,y) call XML_AddAttribute(xmlf, #x, (y))
define(_ATTR,`call XML_AddAttribute(xmlf, "$1", $2)')
