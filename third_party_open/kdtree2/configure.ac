#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(KDTREE 2, 2.13)
AC_CONFIG_SRCDIR([src-f90/kdtree2.f90])

AC_CONFIG_HEADER([config.h])
dnl Store generated files in config dir
AC_CONFIG_AUX_DIR([config])

dnl m4 macro's in m4 directory.
AC_CONFIG_MACRO_DIR([m4])

dnl Initialise automake
dnl set foreign option because we are missing README, INSTALL ...
AM_INIT_AUTOMAKE([1.9 tar-ustar foreign])


## ----------------------- ##
## Libtool initialisation. ##
## ----------------------- ##
LT_PREREQ(2.4)
LT_INIT
dnl LTDL_INIT dnl I don't think this is required anymore
dnl AC_PROG_LIBTOOL
AC_SUBST([LIBTOOL_DEPS])
AC_CONFIG_SUBDIRS(libltdl)
AC_SUBST(LIBLTDL)
AC_SUBST(INCLTDL)

AC_LANG_PUSH([Fortran])
AC_PROG_FC([ifort gfortran])
FCFLAGS=
AX_CHECK_COMPILE_FLAG([-warn all], [FCFLAGS="$FCFLAGS -warn all"], [], []) 
AX_CHECK_COMPILE_FLAG([-static-libcxa], [FCFLAGS="$FCFLAGS -static-libcxa"], [], []) 
AX_CHECK_COMPILE_FLAG([-ip], [FCFLAGS="$FCFLAGS -ip"], [], []) 
AX_CHECK_COMPILE_FLAG([-fpe0], [FCFLAGS="$FCFLAGS -fpe0"], [], []) 
AX_CHECK_COMPILE_FLAG([-fno-alias], [FCFLAGS="$FCFLAGS -fno-alias"], [], []) 
AX_CHECK_COMPILE_FLAG([-Wall], [FCFLAGS="$FCFLAGS -Wall"], [], []) 
AX_CHECK_COMPILE_FLAG([-fomit-framepointers], [FCFLAGS="$FCFLAGS -fomit-framepointers"], [], []) 
FCFLAGS="$FCFLAGS -O3"

AC_LANG_POP

AC_CONFIG_FILES(Makefile
    src-f90/Makefile
)

AC_OUTPUT

