## Makefile.am

# Set the module path.
AM_FCFLAGS	= -J$(top_builddir)/modules/ -I$(top_builddir)/modules/

lib_LTLIBRARIES	= libnested.la
libnested_la_SOURCES = mt95.F90 rand.F90 rand_slice.F90			\
	rand_nested.F90 treedef.F90 tree.F90 treeutil.F90 twomeans.F90	\
	nested.F90 pset.F90 nested_simple.F90 sampling.F90		\
	metropolis.F90 ellipsoidutil.F90 ellipsoid.F90			\
	nested_cluster.F90 distr.F90