lib_LIBRARIES = libMems.a
#lib_LTLIBRARIES = libMems.la

if DEBUG
CXXFLAGS = -Wall -g -D__WXGTK__
else
CXXFLAGS = -O3 -D__WXGTK__
endif

libMems_a_SOURCES = \
BigDiskSuffixArray.cpp       MimHash.cpp \
BigDnaSar.cpp                MimHashEntry.cpp \
MatchFinder.cpp \
MemHash.cpp             SmallDiskSuffixArray.cpp \
MemHashEntry.cpp        SmallDnaSar.cpp \
SuffixArray.cpp \
DiskSuffixArray.cpp          MemScorer.cpp

#check_PROGRAMS = TestRevComp
#TestRevComp_SOURCES = TestRevComp.cpp
#TestRevComp_LDADD = $(top_builddir)/src/libGenome.a

EXTRA_PROGRAMS = Mormems
Mormems_SOURCES = mormems.cpp
Mormems_INCLUDES = -I$(top_srcdir)/include/ -I@CONFIG_LIBGENOMEDIR@/include -I@CONFIG_LIBGENOMEDIR@/include/gn -I@CONFIG_WXWINDOWSDIR@/include/
Mormems_LDADD = $(top_builddir)/src/libMems.a -lGenome -lwx_gtk -L@CONFIG_LIBGENOMEDIR@/src

INCLUDES = -I$(top_srcdir)/include/ -I@CONFIG_LIBGENOMEDIR@/include -I@CONFIG_LIBGENOMEDIR@/include/gn -I@CONFIG_WXWINDOWSDIR@/include/
