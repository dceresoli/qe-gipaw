# Autoconf template for qe-gipaw
AC_INIT(qe-gipaw, 5.0svn, dceresoli@gmail.com, qe-gipaw)

# Option: Quantum-Espresso source directory (QE_SOURCE)
AC_ARG_WITH(qe-source,
   [AC_HELP_STRING([--with-qe-source=DIR],[Quantum-Espresso source directory (QE_SOURCE)])],
   [qe_source=$withval],
   [qe_source=$QE_SOURCE])

# Option: enable band parallelization
AC_ARG_ENABLE(band-parallel,
   [AC_HELP_STRING([--enable-band-parallel],[enable band parallelization (EXPERIMENTAL!)])],
   BAND_PARALLEL="$enableval",
   BAND_PARALLEL="no"
)

# Look for make.sys
AC_CHECK_FILE(${qe_source}/make.sys, [], [AC_MSG_ERROR([Cannot find make.sys in ${qe_source}])])

# Look for Modules/version.f90
AC_CHECK_FILE(${qe_source}/Modules/version.f90, [], [AC_MSG_ERROR([Cannot find version.f90 in ${qe_source}/Modules])])

# Extract version
AC_MSG_CHECKING([Quantum-Espresso version])
version=`grep version_number ${qe_source}/Modules/version.f90 | cut -d\' -f2`
AC_MSG_RESULT([$version])

case $version in
  4.99* ) ;;
  5.0* ) ;;
  *) AC_MSG_ERROR([Cannot compile against this version of quantum-Espresso]) ;;
esac

# process options
dflags_gipaw=''
if test $BAND_PARALLEL = yes ; then 
  dflags_gipaw="${dflags_gipaw} -D__BANDS"
fi

# export settings to generated files
AC_SUBST(qe_source)
AC_SUBST(dflags_gipaw)
AC_OUTPUT([makedeps.sh src/Makefile], [chmod +x makedeps.sh])

# default commands
./makedeps.sh

