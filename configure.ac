# Autoconf template for qe-gipaw
AC_INIT(qe-gipaw, 5.0svn, dceresoli@gmail.com, qe-gipaw)

# Option: Quantum-Espresso source directory (QE_SOURCE)
AC_ARG_WITH(qe-source,
   [AC_HELP_STRING([--with-qe-source=DIR],[Quantum-Espresso source directory (QE_SOURCE)])],
   [qe_source=$withval],
   [qe_source=$QE_SOURCE])

# make path absolute
abspath=`readlink -f ${qe_source} 2>/dev/null` && qe_source=$abspath

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
  5.0.2 ) ;;
  5.0.3 ) ;;
  5.1* ) ;;
  *) AC_MSG_ERROR([Cannot compile against this version of Quantum-Espresso]) ;;
esac

# Extract svn revision
AC_MSG_CHECKING([Quantum-Espresso SVN revision])
revision=`grep svn_revision ${qe_source}/Modules/version.f90 | cut -d\' -f2`
AC_MSG_RESULT([$revision])
good_revision=10590
if test $revision != "unknown" ; then
  if test $revision -lt ${good_revision} ; then
    AC_MSG_ERROR([Cannot compile against this SVN revision of Quantum-Espresso. Please checkout revision >$good_revision])
  fi
fi

# process options
dflags_gipaw=''
if test $BAND_PARALLEL = yes ; then 
  dflags_gipaw="${dflags_gipaw} -D__BANDS"
fi

# export settings to generated files
AC_SUBST(qe_source)
AC_SUBST(dflags_gipaw)
AC_OUTPUT([makedeps.sh src/Makefile src/Makefile-gpu], [chmod +x makedeps.sh])

# default commands
./makedeps.sh

