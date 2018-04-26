## src/Makefile.am

# Set the path for libnested
AM_LDFLAGS	= -L$(top_builddir)/lib/
LDADD		= $(top_builddir)/lib/libnested.la

# Set the module path.
AM_FCFLAGS	= -I$(top_builddir)/modules/ -I$(top_builddir)/lib/

bin_PROGRAMS = test-sk test-sk-cl test-sk-skew test-sk-comp		\
test-sk-gen test-rb test-cl test-slice test-mph test-sk-sys test-sk-3	\
ex_min

test_sk_SOURCES	= test-sk.F90
test_sk_sys_SOURCES	= test-sk-sys.F90
test_sk_cl_SOURCES	= test-sk-cl.F90
test_sk_skew_SOURCES	= test-sk-skew.F90
test_sk_comp_SOURCES	= test-sk-comp.F90
test_sk_3_SOURCES	= test-sk-3.F90
test_sk_gen_SOURCES	= test-sk-gen.F90
test_mph_SOURCES	= test-mph.F90
test_rb_SOURCES	= test-rb.F90
test_cl_SOURCES	= test-cl.F90
ex_min_SOURCES	= ex_min.F90

test_slice_SOURCES = test-slice.F90
