# Autoconf template for ce-gipaw
AC_INIT(ce-gipaw, 4.3svn, dceresoli@gmail.com, ce-gipaw)

# The only option is the quantum-Espresso source directory (QE_SOURCE)
AC_ARG_WITH(qe-source,
   [AC_HELP_STRING([--with-qe-source=DIR],[quantum-Espresso source directory (QE_SOURCE)])],
   [qe_source=$withval],
   [qe_source=$QE_SOURCE])

# Look for make.sys
AC_CHECK_FILE(${qe_source}/make.sys, [], [AC_MSG_ERROR([Cannot find make.sys in ${qe_source}])])

# Look for Modules/version.f90
AC_CHECK_FILE(${qe_source}/Modules/version.f90, [], [AC_MSG_ERROR([Cannot find version.f90 in ${qe_source}/Modules])])

# Extract version
AC_MSG_CHECKING([quantum-Espresso version])
version=`grep version_number ${qe_source}/Modules/version.f90 | cut -d\' -f2`
AC_MSG_RESULT([$version])

case $version in
  "4.3a") ;;
  "4.3b") ;;
  "4.3") ;;
  *) AC_MSG_ERROR([Cannot compile against this version of quantum-Espresso]) ;;
esac

# export settings to generated files
AC_SUBST(qe_source)
AC_OUTPUT([makedeps.sh src/Makefile], [chmod +x makedeps.sh])

# default commands
./makedeps.sh

