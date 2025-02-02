#-------------------------------------------------------------------------------
#   Makefile kdtree_wrapper
#
#   30 nov 2017
#
#   Copyright (C)  Stichting Deltares, 2011-2017.
#-------------------------------------------------------------------------------

include $(top_srcdir)/common.am
noinst_LTLIBRARIES = libkdtreewrapper.la

.FORCE: 

.PHONY: .FORCE

#-------------------------------------------------------------------------------
#Includes

FMODULES += \
   -I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src \
	-I$(top_builddir)/third_party_open/kdtree2/src-f90

#-------------------------------------------------------------------------------
#Sources
libkdtreewrapper_la_SOURCES = \
         kdtreeWrapper.f90
         
#-------------------------------------------------------------------------------
#Libraries